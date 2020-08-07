import re
import math
import time
from queue import Queue

from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.Restriction import Analysis, RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqUtils import GC

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.string_manipulation import *

# fitness function score name in dict
score_host = 'eval_host'
score_restriction = 'eval_rest_sites'
score_repeats = 'eval_repeats'
score_homopolymers = 'eval_homopolymers'
score_hairpins = 'eval_hairpins'
score_gc = 'eval_gc'
score_names = [score_host, score_restriction, score_repeats, score_gc, score_homopolymers, score_hairpins]


def find_num_differences(sequence1: Seq, sequence2: Seq) -> int:
    return sum(1 for a, b in zip(sequence1, sequence2) if a != b)


def eval_host(params: dict, out_q: Queue) -> None:
    out = {score_host: {}}
    for sequence in params['population'].keys():
        score = find_num_differences(
            Seq(sequence, IUPAC.unambiguous_dna),
            params['codon_opt_seq']
        )
        out[score_host][sequence] = [score]
    out_q.put(out)
    return


def find_restriction_sites(restriction_batch: RestrictionBatch, sequence: Seq, linear=True):
    return Analysis(restrictionbatch=restriction_batch, sequence=sequence, linear=linear).full()


def eval_restriction_sites(params: dict, out_q: Queue) -> None:
    out = {score_restriction: {}}
    for sequence in params['population'].keys():
        rest_sites = find_restriction_sites(
            params['restriction_sites'],
            Seq(sequence, IUPAC.unambiguous_dna),
            params['linear']
        )
        score = 0
        for sites in rest_sites.values():
            score += len(sites)
        out[score_restriction][sequence] = [score]
    out_q.put(out)
    return


def eval_repeats(params: dict, out_q: Queue) -> None:
    out = {score_repeats: {}}
    for sequence in params['population'].keys():
        if params['locations']:
            locations = find_repeats(sequence, params['repeat_size'], params['overlapping'])
            score = [get_number_of_repeats_from_dict(locations), locations]
        else:
            if params['overlapping']:
                score = [find_number_of_overlapping_repeats(sequence, params['repeat_size'])]
            else:
                score = [find_number_of_non_overlapping_repeats(sequence, params['repeat_size'])]
        out[score_repeats][sequence] = score
    out_q.put(out)
    return


def eval_homopolymers(params: dict, out_q: Queue) -> None:
    out = {score_homopolymers: {}}
    for sequence in params['population'].keys():
        repeat_and_locations = find_repeats(sequence, params['homopolymer_size'], overlapping=True)
        # remove repeats with multiple letters
        locations = {}
        score = 0
        for k, v in repeat_and_locations.items():
            if len(set(k)) == 1:
                if params['locations']:
                    locations[k] = v
                score += len(v)
        out[score_homopolymers][sequence] = [score, locations]
    out_q.put(out)
    return


def eval_hairpins(params: dict, out_q: Queue) -> None:
    """ Current test does NOT have hairpins
    loops cannot be less than 3 bases long

    ideally 4-8 bases
    longer requires extra rna structures

    for each position
    look at current, current + separation
    then keep checking until stem length of 10result = {separation: set()}
    """
    out = {score_hairpins: {}}
    if params['shortest_loop_length'] < 3:
        raise AttributeError("Loop cannot be smaller than 3 (bio rules)")
    if params['longest_loop_length'] > 30:
        print("loop likely to be too unstable to exist without extra features inside")
    for sequence in params['population'].keys():
        palindrome_locations = find_separated_palindromes(
            sequence,
            params['shortest_loop_length'],
            params['longest_loop_length'],
            params['stem_length']
        )
        score = 0
        for x in palindrome_locations.values():
            score += len(x)
        locations = tuple(palindrome_locations.values()) if params['locations'] else tuple()
        out[score_hairpins][sequence] = [score, locations]
    out_q.put(out)
    return


def get_windowed_gc(sequence: str, window: int) -> list:
    seq = Seq(sequence, IUPAC.unambiguous_dna)
    gc_results = list()
    for idx in range(0, len(seq), window):
        sequence_window = seq[idx:idx + window] if idx + window < len(seq) else seq[idx:len(seq)]
        gc_percent = GC(sequence_window)
        gc_results.append(gc_percent)
    return gc_results


def eval_gc(params: dict, out_q: Queue) -> None:
    out = {score_gc: {}}
    if len(params['gc_params']) < 1:
        raise AttributeError("gc_params should not be empty")
    for gc_param in params['gc_params']:
        if len(gc_param) != 3 or \
                not isinstance(gc_param[0], float) or \
                not isinstance(gc_param[1], float) or \
                not isinstance(gc_param[2], int):
            raise AttributeError("gc_params should contain float low%, float high%, int window size")
    for sequence in params['population'].keys():
        score = 0
        for gc_param in params['gc_params']:
            gc_results = get_windowed_gc(sequence, gc_param[2])
            for result in gc_results:
                if result < gc_param[0] or result > gc_param[1]:
                    score += 1
        out[score_gc][sequence] = [score]
    out_q.put(out)
    return


"""
# are these necessary? ask the bio people
def eval_start_sites(individual, ribosome_binding_sites, table_name="Standard"):
    sequence = getattr(individual, "sequence")
    codon_table = CodonTable.unambiguous_dna_by_name[table_name]

    # find all start codon sites (xTG)
    start_codon_positions = [
        m.start()
        for start_codon in codon_table.start_codons
        for m in re.finditer(start_codon, str(sequence))
    ]

    # None found
    if not len(start_codon_positions):
        return 0

    # check each start site for RBS
    # 18 base pairs upstream of each xTG, ignore 3 bp closest to xTG
    _rbs_offset = 18
    rbs_positions = [
        pos - _rbs_offset for pos in start_codon_positions if pos >= _rbs_offset
    ]
    mutable_seq = sequence.tomutable()

    score = 0
    for rbs_start in rbs_positions:
        # ignore 3 bp closest to xTG
        rbs_stop = rbs_start + _rbs_offset - 3
        rbs_query_seq = str(mutable_seq[rbs_start:rbs_stop])

        # check each unwanted RBS in each potential fragment
        for rbs, site in ribosome_binding_sites.items():
            search = rbs_query_seq.find(site)

            count = 0  # counter to prevent infinite loop
            while search != -1 and count < 10:
                # mutate residues if site is found
                codon_pos = (search + rbs_start) // 3
                for ii in range(2):
                    codon_idx = slice((codon_pos + ii) * 3, (codon_pos + ii + 1) * 3)
                    score += 1

                # reset sequence and search again
                rbs_query_seq = str(mutable_seq[rbs_start: rbs_stop + 3])
                search = rbs_query_seq.find(site)
                count += 1
    return score

def eval_splice_sites(individual):
    sequence = getattr(individual, "sequence")

    def _pass_back_matches(list_of_sites, curr_dna):
        dna = str(curr_dna)
        sites = set(m for expr in list_of_sites for m in re.finditer(expr, dna))
        try:
            sites.remove(None)
        except KeyError:
            pass
        # remove redundancy
        sites = set((site.span(), site[0]) for site in sites)
        codon_bounds = [
            (s[0][0] // 3, -(-s[0][1] // 3)) for s in sorted(sites, key=lambda x: x[0])
        ]
        return codon_bounds

    def _get_splice_sites(curr_dna):
        # donor_sites = _pass_back_matches(splice_donors, curr_dna)
        # acceptor_sites = _pass_back_matches(splice_acceptors, curr_dna)
        # return set(donor_sites + acceptor_sites)
        pass

    return len(_get_splice_sites(sequence.tomutable))
"""

fitness_evals = [
    eval_host,
    eval_repeats,
    eval_restriction_sites,
    eval_homopolymers,
    eval_hairpins,
    eval_gc
]

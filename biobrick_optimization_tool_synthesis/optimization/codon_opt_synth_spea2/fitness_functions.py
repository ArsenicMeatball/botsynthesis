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


def find_num_differences(sequence1: Seq, sequence2: Seq) -> int:
    return sum(1 for a, b in zip(sequence1, sequence2) if a != b)


def eval_host(params: dict, out_q: Queue):
    out = {'eval_host': {}}
    for sequence in params['population'].keys():
        # convert to seq and send it
        out['eval_host'][sequence] = find_num_differences(
            Seq(sequence, IUPAC.unambiguous_dna),
            params['codon_opt_seq']
        )
    out_q.put(out)


def find_restriction_sites(restriction_batch: RestrictionBatch, sequence: Seq, linear=True):
    return Analysis(restrictionbatch=restriction_batch, sequence=sequence, linear=linear).full()


def eval_restriction_sites(params: dict, out_q: Queue):
    out = {'eval_rest_sites': {}}
    for sequence in params['population'].keys():
        rest_sites = find_restriction_sites(
            params['restriction_sites'],
            Seq(sequence, IUPAC.unambiguous_dna),
            params['linear']
        )
        score = 0
        for sites in rest_sites.values():
            score += len(sites)
        out['eval_rest_sites'][sequence] = score
    out_q.put(out)


def eval_repeats(params: dict, out_q: Queue):
    out = {'eval_repeats': {}}
    for sequence in params['population'].keys():
        if params['locations']:
            locations = find_repeats(sequence, params['repeat_size'], params['overlapping'])
            score = [get_number_of_repeats_from_dict(locations), locations]
        else:
            if params['overlapping']:
                score = [find_number_of_overlapping_repeats(sequence, params['repeat_size'])]
            else:
                score = [find_number_of_non_overlapping_repeats(sequence, params['repeat_size'])]
        out['eval_repeats'][sequence] = score
    out_q.put(out)


def eval_homopolymers(params: dict, out_q: Queue):
    out = {'eval_homopolymers': {}}
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
        out['eval_homopolymers'][sequence] = [score, locations]
    out_q.put(out)



"""

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
"""

"""
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


def eval_gc_content(individual, gc):
    sequence = getattr(individual, "sequence")
    window_size = gc.window_size  # tuples are immutable
    # some windows may be expressed as function of the sequence length
    if isinstance(window_size, str) and window_size.startswith("x"):
        window_size = int(float(window_size[1:]) * len(sequence))

    mutable_seq = sequence.tomutable()
    score = 0
    for i in range(len(mutable_seq)):
        window = slice(
            i,
            (i + window_size)
            if (i + window_size) < len(mutable_seq)
            else len(mutable_seq),
        )
        gc_percent = GC(mutable_seq[window]) / 100

        if gc_percent > gc.high:
            score += (gc_percent - gc.high) * 100
        if window.stop is len(mutable_seq):
            break
    return round(score)


def eval_hairpins(individual, stem_length=10):
    sequence = getattr(individual, "sequence")
    mutable_seq = sequence.tomutable()
    score = 0
    for i in range(0, len(mutable_seq) - stem_length, 3):
        stem_seq = mutable_seq[i: i + stem_length].toseq()
        # include wobble base pairing for G-[CT]
        hairpin_pattern = "".join(
            [nt if nt != "C" else "[CT]" for nt in stem_seq.reverse_complement()]
        )
        for hairpin in re.finditer(hairpin_pattern, str(mutable_seq)):
            if not (math.fabs(hairpin.start() - i) < stem_length or stem_length + 2 < hairpin.end() - hairpin.start()):
                score += 1
    return score
"""

fitness_evals = [
    eval_host,
    eval_repeats,
    eval_restriction_sites,
    eval_homopolymers
]
"""

eval_start_sites,

eval_homopolymers,
eval_splice_sites,
eval_gc_content,
eval_hairpins
"""

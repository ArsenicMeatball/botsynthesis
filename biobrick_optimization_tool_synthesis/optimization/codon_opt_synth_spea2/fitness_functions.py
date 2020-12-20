import logging

from queue import Queue

from Bio.Alphabet import IUPAC
from Bio.Restriction import Analysis, RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqUtils import GC

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.density import calculate_density, \
    __DENSITY_KEY__
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.domination import calculate_raw_fitness, \
    __RAW_FITNESS_KEY__
import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.mutations as mut
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.string_functions import \
    find_num_differences, get_number_of_repeats_from_repeats_dict, find_repeats, find_number_of_overlapping_repeats, \
    find_number_of_non_overlapping_repeats, find_separated_palindromes

# fitness function score name in dict
__SCORE_HOST__ = 'eval_host'
__SCORE_RESTRICTION__ = 'eval_rest_sites'
__SCORE_REPEATS__ = 'eval_repeats'
__SCORE_HOMOPOLYMERS__ = 'eval_homopolymers'
__SCORE_HAIRPINS__ = 'eval_hairpins'
__SCORE_GC__ = 'eval_gc'
__SCORE_NAMES__ = [
    __SCORE_HOST__,
    __SCORE_RESTRICTION__,
    __SCORE_REPEATS__,
    __SCORE_GC__,
    __SCORE_HOMOPOLYMERS__,
    __SCORE_HAIRPINS__,
]
__FITNESS_KEY__ = 'fitness'


def calculate_fitness(population: dict):
    """
    TODO: Multiprocessing
    Calculates final fitness
    Fitness = Raw fitness + density
    :param population:
    :return: None, update population
    """
    calculate_raw_fitness(population)
    calculate_density(population)
    for seq_id in population:
        population[seq_id][__FITNESS_KEY__] = population[seq_id][__RAW_FITNESS_KEY__] + \
                                                population[seq_id][__DENSITY_KEY__]


# Specific fitness

def eval_host(params: dict, out_q: Queue) -> None:
    out = {__SCORE_HOST__: {}}
    for seq_id in params['population'].keys():
        score = find_num_differences(
            params['population'][seq_id][mut.__SEQUENCE_KEY__],
            params['codon opt seq']
        )
        out[__SCORE_HOST__][seq_id] = [score]
    out_q.put(out)
    return


def find_restriction_sites(restriction_batch: RestrictionBatch, sequence: Seq, linear=True):
    return Analysis(restrictionbatch=restriction_batch, sequence=sequence, linear=linear).full()


def eval_restriction_sites(params: dict, out_q: Queue) -> None:
    out = {__SCORE_RESTRICTION__: {}}
    for seq_id in params['population'].keys():
        rest_sites = find_restriction_sites(
            params['restriction sites'],
            Seq(
                params['population'][seq_id][mut.__SEQUENCE_KEY__],
                IUPAC.unambiguous_dna
            ),
            params['linear']
        )
        score = [0, list(rest_sites.values())]
        for sites in rest_sites.values():
            score[0] += len(sites)
        out[__SCORE_RESTRICTION__][seq_id] = score
    out_q.put(out)
    return


def eval_repeats(params: dict, out_q: Queue) -> None:
    out = {__SCORE_REPEATS__: {}}
    for seq_id in params['population'].keys():
        if params['locations']:
            locations = find_repeats(
                params['population'][seq_id][mut.__SEQUENCE_KEY__],
                params['repeat size'],
                params['overlapping']
            )
            score = [get_number_of_repeats_from_repeats_dict(locations), locations]
        else:
            if params['overlapping']:
                score = [find_number_of_overlapping_repeats(
                    params['population'][seq_id][mut.__SEQUENCE_KEY__],
                    params['repeat size']
                )]
            else:
                score = [find_number_of_non_overlapping_repeats(
                    params['population'][seq_id][mut.__SEQUENCE_KEY__],
                    params['repeat size']
                )]
        out[__SCORE_REPEATS__][seq_id] = score
    out_q.put(out)
    return


def eval_homopolymers(params: dict, out_q: Queue) -> None:
    out = {__SCORE_HOMOPOLYMERS__: {}}
    for seq_id in params['population'].keys():
        repeat_and_locations = find_repeats(
            params['population'][seq_id][mut.__SEQUENCE_KEY__],
            params['homopolymer size'], overlapping=True
        )
        # remove repeats with multiple letters
        locations = {}
        score = 0
        for k, v in repeat_and_locations.items():
            if len(set(k)) == 1:
                if params['locations']:
                    locations[k] = v
                score += len(v)
        out[__SCORE_HOMOPOLYMERS__][seq_id] = [score, locations]
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
    out = {__SCORE_HAIRPINS__: {}}
    if params['shortest loop length'] < 3:
        raise AttributeError("Loop cannot be smaller than 3 (bio rules)")
    if params['longest loop length'] > 30:
        logging.warning("loop likely to be too unstable to exist without extra features inside")
    for seq_id in params['population'].keys():
        palindrome_locations = find_separated_palindromes(
            params['population'][seq_id][mut.__SEQUENCE_KEY__],
            params['shortest loop length'],
            params['longest loop length'],
            params['stem length']
        )
        logging.debug(palindrome_locations)
        score = 0
        for x in palindrome_locations.values():
            score += len(x)
        out[__SCORE_HAIRPINS__][seq_id] = [score, palindrome_locations]
    out_q.put(out)
    return


def get_windowed_gc(sequence: str, window: int) -> list:
    """
    Gets the
    :param sequence:
    :param window:
    :return: list of lists with idx 0 being the percent and idx 1-2 being the start-finish of the window
    """
    seq = Seq(sequence, IUPAC.unambiguous_dna)
    gc_results = list()
    for idx in range(0, len(seq)):
        sequence_window = seq[idx:idx + window] if idx + window < len(seq) else seq[idx:len(seq)]
        # this returns percent as 90.0 = 90% instead of 0.9
        gc_percent = GC(sequence_window)
        # convert to decimal
        gc_percent /= 100
        gc_results.append(gc_percent)

        if seq[idx:len(seq)] == sequence_window:
            break
    logging.debug(gc_results)
    return gc_results


def eval_gc(params: dict, out_q: Queue) -> None:
    out = {__SCORE_GC__: {}}
    if len(params['gc parameters']) != 3 or \
            not isinstance(params['gc parameters']['min'], float) or \
            not isinstance(params['gc parameters']['max'], float) or \
            not isinstance(params['gc parameters']['window size'], int):
        raise AttributeError("gc params should contain float low%, float high%, int window size")
    for seq_id in params['population'].keys():
        score = 0
        gc_results = get_windowed_gc(
            params['population'][seq_id][mut.__SEQUENCE_KEY__],
            params['gc parameters']['window size']
        )
        for result in gc_results:
            if result < params['gc parameters']['min']:
                score += 1
                logging.debug(
                    'min {0}, actual {1}, was found to be low'.format(
                        params['gc parameters']['min'],
                        result
                    )
                )
            elif result > params['gc parameters']['max']:
                score += 1
                logging.debug(
                    'max {0}, actual {1}, was found to be high'.format(
                        params['gc parameters']['max'],
                        result
                    )
                )
        out[__SCORE_GC__][seq_id] = [score, gc_results]
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

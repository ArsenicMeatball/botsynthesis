import logging

from queue import Queue

from Bio.Restriction import Analysis, RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqUtils import GC

from botsynthesis.spea2.density import calculate_density
from botsynthesis.spea2.domination import calculate_raw_fitness
from botsynthesis.utils.constants import FITNESS_KEY, DENSITY_KEY, \
    RAW_FITNESS_KEY, SCORE_HOST, SCORE_RESTRICTION, SCORE_REPEATS, \
    SCORE_HOMOPOLYMERS, SCORE_HAIRPINS, SCORE_GC, SEQUENCE_KEY
from botsynthesis.utils.string_functions import (
    find_num_differences,
    find_repeats,
    get_number_of_repeats_from_repeats_dict,
    find_number_of_overlapping_repeats,
    find_number_of_non_overlapping_repeats,
    find_separated_palindromes,
)


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
        population[seq_id][FITNESS_KEY] = (
            population[seq_id][RAW_FITNESS_KEY]
            + population[seq_id][DENSITY_KEY]
        )


# Specific fitness


def eval_host(params: dict, out_q: Queue) -> None:
    out = {SCORE_HOST: {}}
    for seq_id in params["population"].keys():
        score = find_num_differences(
            params["population"][seq_id][SEQUENCE_KEY],
            params["codon opt seq"],
        )
        out[SCORE_HOST][seq_id] = [score]
    out_q.put(out)
    return


def find_restriction_sites(
        restriction_batch: RestrictionBatch, sequence: Seq, linear=True
):
    return Analysis(
        restrictionbatch=restriction_batch, sequence=sequence, linear=linear
    ).full()


def eval_restriction_sites(params: dict, out_q: Queue) -> None:
    out = {SCORE_RESTRICTION: {}}
    for seq_id in params["population"].keys():
        rest_sites = find_restriction_sites(
            params["restriction sites"],
            Seq(
                params["population"][seq_id][SEQUENCE_KEY]
            ),
            params["linear"],
        )
        score = [0, list(rest_sites.values())]
        for sites in rest_sites.values():
            score[0] += len(sites)
        out[SCORE_RESTRICTION][seq_id] = score
    out_q.put(out)
    return


def eval_repeats(params: dict, out_q: Queue) -> None:
    out = {SCORE_REPEATS: {}}
    for seq_id in params["population"].keys():
        if params["locations"]:
            locations = find_repeats(
                params["population"][seq_id][SEQUENCE_KEY],
                params["repeat size"],
                params["overlapping"],
            )
            score = [
                get_number_of_repeats_from_repeats_dict(locations),
                locations,
            ]
        else:
            if params["overlapping"]:
                score = [
                    find_number_of_overlapping_repeats(
                        params["population"][seq_id][SEQUENCE_KEY],
                        params["repeat size"],
                    )
                ]
            else:
                score = [
                    find_number_of_non_overlapping_repeats(
                        params["population"][seq_id][SEQUENCE_KEY],
                        params["repeat size"],
                    )
                ]
        out[SCORE_REPEATS][seq_id] = score
    out_q.put(out)
    return


def eval_homopolymers(params: dict, out_q: Queue) -> None:
    out = {SCORE_HOMOPOLYMERS: {}}
    for seq_id in params["population"].keys():
        repeat_and_locations = find_repeats(
            params["population"][seq_id][SEQUENCE_KEY],
            params["homopolymer size"],
            overlapping=True,
        )
        # remove repeats with multiple letters
        locations = {}
        score = 0
        for k, v in repeat_and_locations.items():
            if len(set(k)) == 1:
                if params["locations"]:
                    locations[k] = v
                score += len(v)
        out[SCORE_HOMOPOLYMERS][seq_id] = [score, locations]
    out_q.put(out)
    return


def eval_hairpins(params: dict, out_q: Queue) -> None:
    """Current test does NOT have hairpins
    loops cannot be less than 3 bases long

    ideally 4-8 bases
    longer requires extra rna structures

    for each position
    look at current, current + separation
    then keep checking until stem length of 10result = {separation: set()}
    """
    out = {SCORE_HAIRPINS: {}}
    if params["shortest loop length"] < 3:
        raise AttributeError("Loop cannot be smaller than 3 (bio rules)")
    if params["longest loop length"] > 30:
        logging.warning(
            "loop likely to be too unstable to exist without extra features "
            "inside "
        )
    for seq_id in params["population"].keys():
        palindrome_locations = find_separated_palindromes(
            params["population"][seq_id][SEQUENCE_KEY],
            params["shortest loop length"],
            params["longest loop length"],
            params["stem length"],
        )
        logging.debug(palindrome_locations)
        score = 0
        for x in palindrome_locations.values():
            score += len(x)
        out[SCORE_HAIRPINS][seq_id] = [score, palindrome_locations]
    out_q.put(out)
    return


def get_windowed_gc(sequence: str, window: int) -> list:
    """
    Gets the
    :param sequence:
    :param window:
    :return: list of lists with idx 0 being the percent and idx 1-2 being the
    start-finish of the window
    """
    seq = Seq(sequence)
    gc_results = list()
    for idx in range(0, len(seq)):
        sequence_window = (
            seq[idx: idx + window]
            if idx + window < len(seq)
            else seq[idx: len(seq)]
        )
        # this returns percent as 90.0 = 90% instead of 0.9
        gc_percent = GC(sequence_window)
        # convert to decimal
        gc_percent /= 100
        gc_results.append(gc_percent)

        if seq[idx: len(seq)] == sequence_window:
            break
    logging.debug(gc_results)
    return gc_results


def eval_gc(params: dict, out_q: Queue) -> None:
    out = {SCORE_GC: {}}
    if (
            len(params["gc parameters"]) != 3
            or not isinstance(params["gc parameters"]["min"], float)
            or not isinstance(params["gc parameters"]["max"], float)
            or not isinstance(params["gc parameters"]["window size"], int)
    ):
        raise AttributeError(
            "gc params should contain float low%, float high%, int window size"
        )
    for seq_id in params["population"].keys():
        score = 0
        gc_results = get_windowed_gc(
            params["population"][seq_id][SEQUENCE_KEY],
            params["gc parameters"]["window size"],
        )
        for result in gc_results:
            if result < params["gc parameters"]["min"]:
                score += 1
                logging.debug(
                    "min {0}, actual {1}, was found to be low".format(
                        params["gc parameters"]["min"], result
                    )
                )
            elif result > params["gc parameters"]["max"]:
                score += 1
                logging.debug(
                    "max {0}, actual {1}, was found to be high".format(
                        params["gc parameters"]["max"], result
                    )
                )
        out[SCORE_GC][seq_id] = [score, gc_results]
    out_q.put(out)
    return


fitness_evals = [
    eval_host,
    eval_repeats,
    eval_restriction_sites,
    eval_homopolymers,
    eval_hairpins,
    eval_gc,
]

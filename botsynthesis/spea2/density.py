from math import sqrt

from botsynthesis.utils.constants import DENSITY_KEY, SCORE_HAIRPINS, \
    SCORE_HOMOPOLYMERS, SCORE_GC, SCORE_REPEATS, SCORE_RESTRICTION, \
    SEQUENCE_KEY
from botsynthesis.utils.dict_functions import (
    get_number_of_differences_dict_list,
    sort_dict_by_value_get_list_of_keys,
)
from botsynthesis.utils.list_functions import (
    number_of_differences_between_two_lists,
)
from botsynthesis.utils.string_functions import (
    find_num_differences,
)


def calculate_density(population: dict):
    """
    calculates the density of a sequence based on the kth nearest neighbour
    D = (distance to kth nearest neighbour + 2)^-1,
    where k is sqrt(population + archive size)
    :param population: dict containing all sequences in the population
        and the archive
    :return: None, updates population
    """
    k = round(sqrt(len(population)))
    all_distances = calculate_distances_for_all_sequences(population)
    for seq_id, neighbours in all_distances.items():
        k_nearest_neighbour = get_kth_neirest_neighbour_dict(neighbours, k)
        population[seq_id][DENSITY_KEY] = pow(
            (all_distances[seq_id][k_nearest_neighbour] + 2), -1
        )


def get_kth_neirest_neighbour_dict(neighbours: dict, k: int) -> str:
    """
    gets the kth minimum distance in a dictionary
    :param k: whicvh neighbour to look for, eg k=4, 4th neighbour
    :param neighbours: dict of neighbours {sequence: distance}
    :return: the kth neirest neighbour
    """
    in_order = sort_dict_by_value_get_list_of_keys(neighbours)
    return in_order[k]


def calculate_distances_for_all_sequences(population: dict) -> dict:
    """
    Calculates the distances between all sequences and puts it into a graph
    :param population:
    :return:
    """
    distances = {}
    for seq_id1, seq1_vals in population.items():
        if seq_id1 not in distances:
            distances[seq_id1] = {}
        for seq_id2, seq2_vals in population.items():
            if seq_id2 not in distances:
                distances[seq_id2] = {}
            if seq_id2 not in distances[seq_id1]:
                dist = get_total_distance(
                    seq_id1, seq1_vals, seq_id2, seq2_vals
                )
                distances[seq_id1][seq_id2] = dist
                distances[seq_id2][seq_id1] = dist
    return distances


def get_total_distance(
        seq_id1: str, values_1: dict, seq_id2: str, values_2: dict
) -> float:
    """
    just the sum of all the distances between 2 sequences
    :param values_2:
    :param values_1:
    :param seq_id1:
    :param seq_id2:
    :return:
    """
    total_distance = (
        get_host_distance(
            values_1[SEQUENCE_KEY], values_2[SEQUENCE_KEY]
        ) + get_restriction_distance(
            values_1[SCORE_RESTRICTION][1],
            values_2[SCORE_RESTRICTION][1],
        ) + get_repeat_distance(
            values_1[SCORE_REPEATS][1],
            values_2[SCORE_REPEATS][1],
        ) + get_gc_distance(
            values_1[SCORE_GC][1],
            values_2[SCORE_GC][1],
        ) + get_homopolymers_distance(
            values_1[SCORE_HOMOPOLYMERS][1],
            values_2[SCORE_HOMOPOLYMERS][1],
        ) + get_hairpin_distance(
            values_1[SCORE_HAIRPINS][1],
            values_2[SCORE_HAIRPINS][1],
        )
    )
    return total_distance


def get_host_distance(str1: str, str2: str) -> int:
    """
    Hamming distance
    determine the differences between the two sequences instead of the
    difference with the host
    :argument str1 (str) the first string to compare
    :argument str2 (str) the second string to compare
    :return (int) the number of differences between strings
    """
    return find_num_differences(str1, str2)


def get_restriction_distance(locations1: list, locations2: list) -> int:
    """
    determine the differences in number of restriction sites and their
    locations
    :param locations1: (list) list containing lists of locations where a
    restriction site occurs
    :param locations2: (list) list containing lists of locations where a
    restriction site occurs
    expects [[idx],[],[idx, idx, idx] ... ]
    :return: (int) the number of differences between the two lists
    """
    differences = 0
    # for each restriction enzyme
    for enz_left, enz_right in zip(locations1, locations2):
        differences = number_of_differences_between_two_lists(
            enz_left, enz_right
        )
    return differences


def get_repeat_distance(
        repeat_n_locations1: dict, repeat_n_locations2: dict
) -> int:
    """
    determine the differences in the repeats and their locations between the
    two dictionaries
    expects {'repeat_sequence':[idx, idx], 'seq':[idx] ... }
    :param repeat_n_locations1: the sequences and the locations and
    the indices they are in
    :param repeat_n_locations2: the sequences and the locations and
    the indices they are in
    :return: number of differences between the dictionaries
    """
    return get_number_of_differences_dict_list(
        repeat_n_locations1, repeat_n_locations2
    )


def get_gc_distance(percents1: list, percents2: list) -> float:
    """
    determine the differences in the percentages between the lists
    :param percents1: first list containing percentages
    :param percents2: second list containign percentages
    # must be of the same length
    :return: the sum of the differences
    """
    sum_difference = 0
    for percent1, percent2 in zip(percents1, percents2):
        sum_difference += abs(percent1 - percent2)
    return sum_difference


def get_homopolymers_distance(
        homo_n_locations1: dict, homo_n_locations2: dict
) -> int:
    """
    determine the differences between the two dictionaries
    expects {'repeat_sequence':[idx, idx], 'seq':[idx] ... }
    :return: number of differences between the dictionaries
    """
    return get_number_of_differences_dict_list(
        homo_n_locations1, homo_n_locations2
    )


def get_hairpin_distance(
        hairpin_lengths_n_locations1: dict, hairpin_lengths_n_locations2: dict
) -> int:
    """

    :param hairpin_lengths_n_locations1:
    :param hairpin_lengths_n_locations2:
    :return:
    """
    return get_number_of_differences_dict_list(
        hairpin_lengths_n_locations1, hairpin_lengths_n_locations2
    )

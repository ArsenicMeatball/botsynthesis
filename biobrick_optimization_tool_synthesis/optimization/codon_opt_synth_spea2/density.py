# k = sqrt(sample size)
# D = 1 / (kth distance + 2)
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.dict_functions import \
    get_number_of_differences_dict_list
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions import score_names
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.list_functions import \
    number_of_differences_between_two_lists
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.string_functions import \
    find_num_differences

score_names = [score_host, score_restriction, score_repeats, score_gc, score_homopolymers, score_hairpins]


def get_host_distance(str1: str, str2: str) -> int:
    """
    determine the differences between the two sequences instead of the difference with the host
    :argument str1 (str) the first string to compare
    :argument str2 (str) the second string to compare
    :return (int) the number of differences between strings
    """
    return find_num_differences(str1, str2)


def get_restriction_distance(locations1: list, locations2: list) -> int:
    """
    determine the differences in restriction site locations
    :param locations1: (list) list containing lists of locations where a restriction site occurs
    :param locations2: (list) list containing lists of locations where a restriction site occurs (same order as the above)
    expects [[idx],[],[idx, idx, idx] ... ]
    :return: (int) the number of differences between the two lists
    """
    differences = 0
    # for each restriction enzyme
    for idx in range(1, len(locations1)):
        differences = number_of_differences_between_two_lists(locations1[idx], locations2[idx])
    return differences


def get_repeat_distance(repeat_n_locations1: dict, repeat_n_locations2: dict) -> int:
    """
    determine the differences between the two dictionaries
    expects {'repeat_sequence':[idx, idx], 'seq':[idx] ... }
    :param repeat_n_locations1: the sequences and the locations and the indices they are in
    :param repeat_n_locations2: the sequences and the locations and the indices they are in
    :return: number of differences between the dictionaries
    """
    return get_number_of_differences_dict_list(repeat_n_locations1, repeat_n_locations2)


def get_gc_distance(percents1: list, percents2: list) -> float:
    """
    determine the differences in the percentages between the lists
    :param percents1: first list containing percentages
    :param percents2: second list containign percentages
    # must be of the same length
    :return: the sum of the differences
    """
    sum_difference = 0
    for idx in range(len(percents1)):
        sum_difference += abs(percents1[idx] - percents2[idx])
    return sum_difference


def get_homopolymers_distance(homo_n_locations1: dict, homo_n_locations2: dict) -> int:
    """
    determine the differences between the two dictionaries
    expects {'repeat_sequence':[idx, idx], 'seq':[idx] ... }
    :return: number of differences between the dictionaries
    """
    return get_number_of_differences_dict_list(homo_n_locations1, homo_n_locations2)


def get_hairpin_distance(hairpin_lengths_n_locations1:dict, hairpin_lengths_n_locations2:dict) -> int:
    """

    :param hairpin_lengths_n_locations1:
    :param hairpin_lengths_n_locations2:
    :return:
    """
    return get_number_of_differences_dict_list(hairpin_lengths_n_locations1, hairpin_lengths_n_locations2)

if __name__ == '__main__':
    x = ()
    print(x)
    print(len(x))

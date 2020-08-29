from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.list_functions import \
    differences_between_two_lists


def get_differences_dict_list(dict1: dict, dict2: dict) -> dict:
    """
    grabs the differences between two dicts
    :param dict1:
    :param dict2:
    :return: a dict containing only the uncommon values
    """
    result_dict = {}
    for key in dict1.keys():
        if key not in dict2:
            result_dict[key] = dict1[key]
        else:
            result_dict[key] = differences_between_two_lists(dict1[key], dict2[key])
    for key in dict2.keys():
        if key not in dict1:
            result_dict[key] = dict2[key]
    return result_dict


def get_number_of_differences_dict_list(dict1: dict, dict2: dict) -> int:
    """
    Gets number of differences for a dict where the values are lists
    {'k':[list], 'k':[list]...}
    :param dict1:
    :param dict2:
    :return: num differences
    """
    differences = 0
    only_diffs = get_differences_dict_list(dict1, dict2)
    for value in only_diffs.values():
        differences += len(value)
    return differences


def find_key_n_minimum_value(d: dict) -> list:
    """
    gets minimum value key(s)
    :param d: dict with values
    :return: key(s) of minimum
    """
    temp = min(d.values())
    res = [key for key in d if d[key] == temp]
    return res


def sort_dict_by_value_get_list_of_keys(d: dict) -> list:
    """
    sorts dict by value and returns list of keys
    :param d: dict with values
    :return: keys in order
    """
    return sorted(d.keys(), key=d.get)


def sort_dict_by_value(d: dict) -> list:
    """
    sorts dict by value and returns list of values
    :param d: dict with values
    :return: values in order
    """
    return sorted(d.values())

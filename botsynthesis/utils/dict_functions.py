from botsynthesis.utils.list_functions import (
    differences_between_two_lists,
)


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
            only_differences = differences_between_two_lists(
                dict1[key], dict2[key]
            )
            if len(only_differences) > 0:
                result_dict[key] = only_differences
    for key in dict2.keys():
        if key not in dict1:
            result_dict[key] = dict2[key]
    return result_dict


def get_number_of_differences_dict_list(dict1: dict, dict2: dict) -> int:
    """
    Gets number of differences for a dict where the values are lists
    {'k':[list], 'k':[list]...}
    empty dicts that exixst only on 1 side do not count
    :param dict1:
    :param dict2:
    :return: num differences
    """
    differences = 0
    only_diffs = get_differences_dict_list(dict1, dict2)
    for value in only_diffs.values():
        differences += len(value)
    return differences


def find_keys_of_minimum_value(d: dict) -> list:
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


def invert_dict(d: dict, is_naive=False) -> dict:
    """
    reverses a dict
    :param is_naive: True for pretending like it is a 1-1 mapping
    (may lose entries), False to have sets for keys with
        multiple mappings
    :param d: the dict to invert
    :return: inverted dict
    """
    if is_naive:
        return {v: k for k, v in d.items()}
    inv_dict = {}
    for k, v in d.items():
        if v not in inv_dict:
            inv_dict[v] = set()
        inv_dict[v].add(k)
    return inv_dict

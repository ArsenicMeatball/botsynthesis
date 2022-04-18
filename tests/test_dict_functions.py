import pytest

from Bio.Data import CodonTable


# TODO: Test fail paths
from botsynthesis.dict_functions import get_differences_dict_list, \
    get_number_of_differences_dict_list, find_keys_of_minimum_value, \
    sort_dict_by_value_get_list_of_keys, sort_dict_by_value, invert_dict


@pytest.fixture
def d1():
    return {
        1: ["123", "23"],
        2: ["mm", 3, 5, 4, 6, 3],
        3: [],
        5: [1, 2, 3, 4, 5],
    }


@pytest.fixture
def d2():
    return {1: ["23"], 2: ["mm", 3, 5, 4, 6, 3], 4: [], 5: [1, 2]}


def test_get_differences_dict_list(d1, d2):
    expected = {1: ["123"], 3: [], 4: [], 5: [3, 4, 5]}
    actual = get_differences_dict_list(d1, d2)
    assert actual == expected


def test_get_number_of_differences_dict_list(d1, d2):
    expected = 4
    actual = get_number_of_differences_dict_list(d1, d2)
    assert actual == expected


def test_find_keys_of_minimum_value():
    d = {
        1: 78787787787,
        3: 3,
        4: 7,
        "aohd": 92000,
        "ajfp": 4,
        "we": 5,
        "q": 6,
        "eoqafi": 3,
        "o0aihf": 4,
        "ofajiw": 3,
    }
    expected = [3, "eoqafi", "ofajiw"]
    actual = find_keys_of_minimum_value(d)
    assert actual == expected


def test_sort_dict_by_value_get_list_of_keys():
    d = {
        1: 78787787787,
        3: 3,
        4: 7,
        "aohd": 92000,
        "ajfp": 4,
        "we": 5,
        "q": 6,
    }
    expected = [3, "ajfp", "we", "q", 4, "aohd", 1]
    actual = sort_dict_by_value_get_list_of_keys(d)
    assert actual == expected


def test_sort_dict_by_value():
    d = {
        1: 78787787787,
        3: 3,
        4: 7,
        "aohd": 92000,
        "ajfp": 4,
        "we": 5,
        "q": 6,
    }
    expected = [3, 4, 5, 6, 7, 92000, 78787787787]
    actual = sort_dict_by_value(d)
    assert actual == expected


def test_invert_dict():
    table = CodonTable.generic_by_id[1].forward_table
    print(table)
    actual = invert_dict(table)
    expected = {}
    for codon, amino in table.items():
        if amino not in expected:
            expected[amino] = set()
        expected[amino].add(codon)
    assert actual == expected

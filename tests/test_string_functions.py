import pytest

from botsynthesis.all_algorithm_parameters import algorithm_params
from botsynthesis.string_functions import find_number_of_overlapping_repeats, \
    find_repeats, get_number_of_repeats_from_repeats_dict, \
    find_number_of_non_overlapping_repeats, find_separated_palindromes, \
    find_num_differences


@pytest.fixture
def sequences():
    return [
        "123123123456456456789789789",
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        str(algorithm_params["codon opt seq"]),
    ]


@pytest.fixture
def low_repeat_sizes():
    return [n for n in range(-1, 1)]


@pytest.fixture
def repeat_sizes():
    return [n for n in range(2, 11)]


def test_overlapping_repeats(sequences, repeat_sizes):
    for sequence in sequences:
        for n in repeat_sizes:
            result1 = find_number_of_overlapping_repeats(
                string=sequence, min_repeat_size=n
            )
            repeats = find_repeats(
                string=sequence, min_repeat_size=n, overlapping=True
            )
            result2 = get_number_of_repeats_from_repeats_dict(repeats)
            assert result1 == result2


def test_non_overlapping_repeats(sequences, repeat_sizes):
    for sequence in sequences:
        for n in repeat_sizes:
            result1 = find_number_of_non_overlapping_repeats(
                string=sequence, min_repeat_size=n
            )
            repeats = find_repeats(
                string=sequence, min_repeat_size=n, overlapping=False
            )
            result2 = get_number_of_repeats_from_repeats_dict(repeats)
            assert result1 == result2


def test_small_repeat_size(sequences, low_repeat_sizes):
    for sequence in sequences:
        for n in low_repeat_sizes:
            with pytest.raises(ArithmeticError):
                find_number_of_non_overlapping_repeats(sequence, n)
            with pytest.raises(ArithmeticError):
                find_number_of_overlapping_repeats(sequence, n)
            with pytest.raises(ArithmeticError):
                find_repeats(sequence, n)


def test_large_min_repeat_size(sequences):
    increase = 5
    for sequence in sequences:
        # a repeat size larger than sequence should yield 0 or empty dict
        result = find_number_of_non_overlapping_repeats(
            string=sequence, min_repeat_size=len(sequence) + increase
        )
        assert result == 0
        result = find_number_of_overlapping_repeats(
            string=sequence, min_repeat_size=len(sequence) + increase
        )
        assert result == 0
        result = find_repeats(
            string=sequence, min_repeat_size=len(sequence) + increase
        )
        assert result == {}


def test_find_separated_palindromes_pass():
    palindrome = "1QWERTYUIOP788978POIUYTREWQ1ASDFGHJKL:URT:LKJHGFDSA"
    expected = {
        6: [1, "788978", 26],
        3: [28, "URT", 50],
        8: [0, "P788978P", 27],
    }
    actual = find_separated_palindromes(palindrome)
    assert actual == expected
    palindrome = ""
    expected = {}
    actual = find_separated_palindromes(palindrome)
    assert actual == expected


def test_find_num_differences_pass():
    string1 = "asdfghjklqwertyuiop"
    string2 = "aadfhhjklewervbuiop"
    expected = 5
    actual = find_num_differences(string1, string2)
    assert actual == expected


def test_get_number_of_repeats_from_dict_pass():
    # TODO
    assert True

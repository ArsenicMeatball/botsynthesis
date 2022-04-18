import pytest
from Bio import Restriction

import botsynthesis.cleaning as cleaning


def test_clean_sequence_pass():
    expected = "ACTGGCGTACGTCATGCATCTACACGATCG"
    sequence = "act ggc gt acg tca tg ca tc tac ac g at cg"
    actual = cleaning.clean_sequence(sequence)
    assert expected == actual


def test_clean_sequence_fail():
    sequence = {"aaaaa": "aaaaaa", "bbbbbbb": "bbbbb"}
    with pytest.raises(AttributeError):
        cleaning.clean_sequence(sequence)


def test_get_rest_enzymes_from_list_pass():
    input_list = ["NdeI", "XhoI", "HpaI", "PstI", "EcoRV", "NcoI", "BamHI"]
    expected = Restriction.RestrictionBatch(
        [Restriction.AllEnzymes.get(enz) for enz in input_list]
    )
    actual = cleaning.get_rest_enzymes_from_list(input_list)
    assert expected == actual


def test_get_rest_enzymes_from_list_fail():
    input_list = [
        "REEEEEEEEEEEEEEEEE",
        "XhoI",
        "HpaI",
        "PstI",
        "EcoRV",
        "NcoI",
        "BamHI",
    ]
    with pytest.raises(ValueError):
        cleaning.get_rest_enzymes_from_list(input_list)


def test_get_rest_enzymes_from_string_pass():
    input_list = "XhoI HpaI PstI EcoRV NcoI BamHI"
    expected = Restriction.RestrictionBatch(
        [Restriction.AllEnzymes.get(enz) for enz in input_list.split()]
    )
    actual = cleaning.get_rest_enzymes_from_string(input_list)
    assert expected == actual


def test_get_rest_enzymes_from_string_fail():
    input_list = "REEEEEEEEEEEEEEEEE XhoI HpaI PstI EcoRV NcoI BamHI"
    with pytest.raises(ValueError):
        cleaning.get_rest_enzymes_from_string(input_list)

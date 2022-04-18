import pytest

from botsynthesis.codon_optimization_host import *


def test_determine_ideal_codon_optimized_sequence_pass():
    sequence = Seq("ATRPKKTT")
    usage = {
        "A": {"aaa": 0.9, "bbb": 0.1},
        "T": {"ttt": 0.2, "bbb": 0.1},
        "R": {"rrr": 0.5, "bbb": 0.5},
        "P": {"ppp": 0.9, "bbb": 0.1},
        "K": {"kkk": 0.9, "bbb": 0.1},
    }
    expected = "aaatttrrrpppkkkkkktttttt"
    actual = determine_ideal_codon_optimized_sequence(sequence, usage)
    assert actual == expected


def test_determine_ideal_codon_optimized_sequence_fail():
    sequence = Seq("ATRPKKTT")
    # wrong keys
    usage = {}
    with pytest.raises(KeyError):
        determine_ideal_codon_optimized_sequence(sequence, usage)
    # no value
    usage = {
        "A": {},
        "T": {"ttt": 0.2, "bbb": 0.1},
        "R": {"rrr": 0.5, "bbb": 0.5},
        "P": {"ppp": 0.9, "bbb": 0.1},
        "K": {"kkk": 0.9, "bbb": 0.1},
    }
    with pytest.raises(ValueError):
        determine_ideal_codon_optimized_sequence(sequence, usage)
    # bad numerical value
    usage = {'A': {'aaa': 'REEE', 'bbb': 'ARRR'},
             'T': {'ttt': 0.2, 'bbb': 0.1}, 'R': {'rrr': 0.5, 'bbb': 0.5},
             'P': {'ppp': 0.9, 'bbb': 0.1},
             'K': {'kkk': 0.9, 'bbb': 0.1}}
    # todo check this at start
    # with pytest.raises(ValueError):
    #     determine_ideal_codon_optimized_sequence(sequence, usage)

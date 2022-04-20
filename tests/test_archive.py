import copy
import botsynthesis.spea2.archive as arch
import pytest


@pytest.fixture
def seqs():
    return [{"1": 5}, {"1": 4}, {"1": 3}, {"1": 3}, {"1": 2.5}, {"1": 1}]


@pytest.fixture
def pops(seqs):
    return [
        {"1": seqs[0]},
        {
            "1": seqs[0],
            "2": seqs[1],
            "3": seqs[2],
            "4": seqs[3],
            "5": seqs[4],
        }, {
            "1": seqs[0],
            "2": seqs[1],
            "3": seqs[2],
            "4": seqs[3],
            "5": seqs[4],
            "6": seqs[5],
        },
        {"6": seqs[5]},
        {
            "6": seqs[4],
            "5": seqs[3],
            "4": seqs[2],
            "3": seqs[1],
            "2": seqs[0],
        }, {
            "6": seqs[5],
            "5": seqs[4],
            "4": seqs[3],
            "3": seqs[2],
            "2": seqs[1],
            "1": seqs[0],
        }
    ]


def test_build_archive_not_existing():
    pop = {}
    with pytest.raises(ValueError):
        arch.build_archive(pop, len(pop) - 1, "1")
    with pytest.raises(ValueError):
        arch.build_archive(pop, 0, "1")
    with pytest.raises(ValueError):
        arch.build_archive(pop, -1, "1")


def test_build_archive(pops):
    for pop in pops:
        archive = arch.build_archive(
            pop, len(pop) + 1, "1"
        )
        assert len(archive) == len(pop)

        archive = arch.build_archive(
            pop, len(pop), "1"
        )
        assert len(archive) == len(pop)

        archive = arch.build_archive(
            pop, len(pop) - 1, "1"
        )
        assert len(archive) == len(pop) - 1

        archive = arch.build_archive(
            pop, 1, "1"
        )
        assert len(archive) == 1

        archive = arch.build_archive(
            pop, 0, "1"
        )
        assert len(archive) == 0

        with pytest.raises(IndexError):
            arch.build_archive(
                pop, -1, "1"
            )

        with pytest.raises(KeyError):
            arch.build_archive(pop, len(pop) - 1, "not a key")


def test_add_dominated_to_archive_pass(seqs, pops):
    # should add the best dominated members from the population
    # into the archive
    expected = {
        "5": seqs[4],
        "6": seqs[5],
        "4": seqs[3],
        "3": seqs[2],
    }
    actual = {"5": seqs[4], "6": seqs[5]}
    pop = copy.deepcopy(pops[5])
    arch.add_dominated_to_archive(
        actual, pop, 4, "1"
    )
    assert actual == expected


def test_add_dominated_to_archive_pass_empty_archive(seqs, pops):
    # same thing but with an empty dict
    expected = {
        "5": seqs[4],
        "6": seqs[5],
        "4": seqs[3],
        "3": seqs[2],
    }
    actual = {}
    pop = copy.deepcopy(pops[5])
    arch.add_dominated_to_archive(
        actual, pop, 4, "1"
    )
    assert actual == expected


def test_add_dominated_to_archive_fail(pops):
    # population too small
    archive = population = {}
    desired_size = 10000000
    with pytest.raises(ValueError):
        arch.add_dominated_to_archive(archive, population, desired_size, "1")
    # wrong key
    with pytest.raises(KeyError):
        arch.add_dominated_to_archive(archive, pops[5], 1, "not a key")


def test_truncate_archive_pass(seqs, pops):
    # should remove worst (highest scores) until correct size
    actual = copy.deepcopy(pops[5])
    desired_size = 0
    expected = {}
    arch.truncate_archive(actual, desired_size, "1")
    assert actual == expected

    actual = copy.deepcopy(pops[5])
    desired_size = 4
    expected = {
        "5": seqs[4],
        "6": seqs[5],
        "4": seqs[3],
        "3": seqs[2],
    }
    arch.truncate_archive(actual, desired_size, "1")
    assert actual == expected


def test_truncate_archive_fail(pops):
    # archive smaller than desired size
    with pytest.raises(ValueError):
        arch.truncate_archive({}, 10000, "1")
    # wrong key
    with pytest.raises(KeyError):
        arch.truncate_archive(pops[5], 0, "not a key")

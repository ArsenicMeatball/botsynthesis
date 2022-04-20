import pytest
from Bio.Seq import Seq
from Bio.Data import CodonTable

import botsynthesis.spea2.mutations as mut
import botsynthesis.utils.dict_functions as dictf
import botsynthesis.utils.string_functions as strf
import botsynthesis.spea2.fitness_functions as fitf

# TODO finish tests
from botsynthesis.utils.constants import FITNESS_KEY, SEQUENCE_KEY


@pytest.fixture()
def codon_table():
    return CodonTable.unambiguous_dna_by_id[1]


def test_mutate_codon_pass():
    # regular case
    initial_codon = "aac".upper()
    codon_table = CodonTable.unambiguous_dna_by_id[1]
    new_codon = mut.mutate_codon(initial_codon, codon_table)
    assert initial_codon != new_codon
    possible_codons = dictf.invert_dict(codon_table.forward_table)[
        codon_table.forward_table[initial_codon]
    ]
    assert new_codon in possible_codons
    # works with another translation dict
    codon_table = CodonTable.unambiguous_dna_by_id[2]
    new_codon = mut.mutate_codon(initial_codon, codon_table)
    assert new_codon != initial_codon
    possible_codons = dictf.invert_dict(codon_table.forward_table)[
        codon_table.forward_table[initial_codon]
    ]
    assert new_codon in possible_codons
    # no other codons case
    initial_codon = "tgg".upper()
    codon_table = CodonTable.unambiguous_dna_by_id[1]
    new_codon = mut.mutate_codon(initial_codon, codon_table)
    assert initial_codon == new_codon


def test_mutate_codon_fail():
    # smol codon
    initial_codon = "ac"
    codon_table = "not a codon tbale object"
    with pytest.raises(ValueError):
        mut.mutate_codon(initial_codon, codon_table)
    # bad codon usage
    initial_codon = "acg"
    with pytest.raises(AttributeError):
        mut.mutate_codon(initial_codon, codon_table)


def test_mutate_seq_pass():
    seq1 = "actagctacggaattagcagaagcaatgctagc".upper()
    # regular, make sure length doesnt change
    mutation_chance = 0.5
    codon_table = CodonTable.unambiguous_dna_by_id[1]
    mutant = mut.mutate_seq(seq1, mutation_chance, codon_table)
    assert len(mutant) == len(seq1)
    # make sure the number of differences makes sense at 100% mutation
    mutation_chance = 1
    mutant = mut.mutate_seq(seq1, mutation_chance, codon_table)
    minimum_number_of_differences = 11
    maximum_number_of_differences = 20
    actual_number_of_differences = strf.find_num_differences(seq1, mutant)
    assert actual_number_of_differences >= minimum_number_of_differences
    assert actual_number_of_differences <= maximum_number_of_differences
    # make sure the number of differences makes sense at 0% mutation
    mutation_chance = 0
    mutant = mut.mutate_seq(seq1, mutation_chance, codon_table)
    actual_number_of_differences = strf.find_num_differences(seq1, mutant)
    assert actual_number_of_differences == 0


def test_mutate_seq_fail(codon_table):
    # sequence too small
    seq = ""
    mutation_chance = -1
    with pytest.raises(ValueError):
        mut.mutate_seq(seq, mutation_chance, codon_table)
    # sequence not multiple of 3
    seq = "reeeeeeeeeeeeeeeee"
    with pytest.raises(ValueError):
        mut.mutate_seq(seq, mutation_chance, codon_table)
    # mutation chance too low
    seq = "ACTGTA"
    with pytest.raises(ValueError):
        mut.mutate_seq(seq, mutation_chance, codon_table)
    # mutation chance too high
    mutation_chance = 2
    with pytest.raises(ValueError):
        mut.mutate_seq(seq, mutation_chance, codon_table)
    # sequence not DNA
    seq = "reeeeeeee"
    mutation_chance = 0.5
    with pytest.raises(KeyError):
        mut.mutate_seq(seq, mutation_chance, codon_table)


def test_initialize_population_pass(codon_table):
    # can it create a good population with default parameters
    seq = "actagctacggaattagcagaagcaatgctagc".upper()
    desired_size = 10
    population = mut.initialize_population(
        desired_population_size=desired_size,
        parent_sequence=seq,
        mutation_chance=0.7, codon_table=codon_table
    )
    assert len(population) == desired_size
    # can it create a good population with tougher parameters?
    population = mut.initialize_population(
        desired_population_size=desired_size,
        parent_sequence=seq,
        mutation_chance=0.2, codon_table=codon_table
    )
    assert len(population) == desired_size


def test_initialize_population_fail(codon_table):
    # bad population
    seq = "actagctacggaattagcagaagcaatgctagc".upper()
    with pytest.raises(ValueError):
        mut.initialize_population(0, seq, 0.7, codon_table)


def test_tournament_selection_without_replacement_pass():
    key = "key"
    population = {
        1: {key: 9099990.393},
        2: {key: 9099990.392},
        3: {key: 9099990.394},
    }
    winner = mut.tournament_selection_without_replacement(
        population, fitness_key_name=key, minima=True, n_ary=2
    )
    assert winner in population.keys()
    winner = mut.tournament_selection_without_replacement(
        population, n_ary=100, fitness_key_name=key, minima=True)
    assert winner == 2
    winner = mut.tournament_selection_without_replacement(
        population, n_ary=100, minima=False, fitness_key_name=key
    )
    assert winner == 3


def test_tournament_selection_without_replacement_fail():
    key = "key"
    population = {
        1: {key: 9099990.393},
        2: {key: 9099990.392},
        3: {key: 9099990.394},
    }
    with pytest.raises(ValueError):
        mut.tournament_selection_without_replacement(population, 1, True,
                                                     FITNESS_KEY)


def test_generate_mating_pool_from_archive_pass():
    # can generate a mating pool from an archive
    key = "key"
    archive = {
        1: {key: 1},
        2: {key: 2},
        3: {key: 3},
        4: {key: 4},
        5: {key: 5},
    }
    desired_size = 3
    mating_pool = mut.generate_mating_pool_from_archive(
        archive, desired_size, key
    )
    assert len(mating_pool) == desired_size
    assert mating_pool.items() <= archive.items()
    # default key works
    key = FITNESS_KEY
    archive = {
        1: {key: 1},
        2: {key: 2},
        3: {key: 3},
        4: {key: 4},
        5: {key: 5},
    }
    mating_pool = mut.generate_mating_pool_from_archive(
        archive, desired_size, key
    )
    assert len(mating_pool) == desired_size
    assert mating_pool.items() <= archive.items()


def test_generate_mating_pool_from_archive_fail():
    # smol archive
    with pytest.raises(ValueError):
        mut.generate_mating_pool_from_archive({}, 5, FITNESS_KEY)
    # smol mating pool
    archive = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
    with pytest.raises(ValueError):
        mut.generate_mating_pool_from_archive(archive, 0, FITNESS_KEY)
    # big mating pool
    with pytest.raises(ValueError):
        mut.generate_mating_pool_from_archive(archive, 100,
                                              FITNESS_KEY)


def test_recombine_pass():
    # test that from 2 strings, the new one is created
    # new one is same aa
    s1 = "actagctacggaattagcagaagcaatgctagc"
    s2 = "acaagttatgggatctcgcgatcgaacgcaagt"
    aa1 = Seq(s1).translate()
    aa2 = Seq(s2).translate()
    assert aa2 == aa1
    for num in range(1, 20):
        r = mut.recombine_dna_sequence(s1, s2, num, 3)
        raa = Seq(r).translate()
        assert raa == aa1
        for pos in range(len(r)):
            assert r[pos] == s1[pos] or r[pos] == s2[pos]
        if num > 10:
            assert r != s1
    new_seq = mut.recombine_dna_sequence(s1, s1, 10, 3)
    assert s1 == new_seq


def test_recombine_fail():
    s1 = "actagctacggaattagcagaagcaatgctagc"
    s2 = "acaagttatgggatctcgcgatcgaacgcaagt"
    # detect invalid number of sites
    with pytest.raises(ValueError):
        mut.recombine_dna_sequence(s1, s2, 0, 3)
    # detect unequal sequences
    with pytest.raises(ValueError):
        mut.recombine_dna_sequence(s1, s2[:-2], 5, 3)


def test_generate_population_from_archive_pass(codon_table):
    # can it run regularly with default vars
    s1 = "actagctacggaattagcCGCagcaatgctagc".upper()
    s2 = "ACAAGCTACGGTATTTCTCGTAGTAATGCTTCA"
    s3 = "ACTAGCTACGGCATCAGCCGCAGCAATGCTAGC"
    s4 = "ACTAGTTACGGGATTTCTCGGAGCAACGCAAGC"
    s5 = "ACGTCATACGGAATTAGCCGAAGTAATGCTAGC"
    s6 = "ACCAGCTATGGAATTAGTCGAAGCAACGCGAGC"
    expected_aa = Seq(s1).translate()
    desired_population_size = 15
    recomb_chance = mutation_chance = 0.2
    num_recombination_sites = mut.get_rec_sites_for_len(
        len(s1), recomb_chance
    )
    fit_key = FITNESS_KEY
    seq_key = SEQUENCE_KEY
    archive = {
        1: {seq_key: s1, fit_key: 0.2},
        2: {seq_key: s2, fit_key: 1.9},
        3: {seq_key: s3, fit_key: 0.1},
        4: {seq_key: s4, fit_key: 0.7},
        5: {seq_key: s5, fit_key: 1.2},
        6: {seq_key: s6, fit_key: 0.9},
    }

    desired_mating_pool_size = 4
    population = mut.generate_population_from_archive(
        archive,
        desired_mating_pool_size,
        num_recombination_sites,
        mutation_chance,
        desired_population_size, FITNESS_KEY, codon_table
    )
    assert len(population) == desired_population_size
    for ind in population.values():
        if ind[seq_key] in [s1, s2, s3, s4, s5, s6]:
            assert expected_aa == Seq(ind[seq_key]).translate()
    # can it run regularly with custom vars
    codon_table = CodonTable.unambiguous_dna_by_id[2]
    expected_aa = Seq(s1).translate(
        codon_table)
    fit_key = "reeeee"
    new_archive = {
        1: {seq_key: s1, fit_key: 0.2},
        2: {seq_key: s2, fit_key: 1.9},
        3: {seq_key: s3, fit_key: 0.1},
        4: {seq_key: s4, fit_key: 0.7},
        5: {seq_key: s5, fit_key: 1.2},
        6: {seq_key: s6, fit_key: 0.9},
    }
    population = mut.generate_population_from_archive(
        new_archive,
        desired_mating_pool_size,
        num_recombination_sites,
        mutation_chance,
        desired_population_size,
        fit_key,
        codon_table,
    )
    assert len(population) == desired_population_size
    for ind in population.values():
        assert Seq(ind[seq_key]).translate(
            codon_table
        ) == expected_aa


def test_get_rec_sites_for_len_pass():
    # TODO
    assert True


def test_get_rec_sites_for_len_fail():
    # TODO
    assert True

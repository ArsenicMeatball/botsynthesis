import logging
import unittest
from typing import Dict, Union

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Data import CodonTable

import botsynthesis.mutations as mut
import botsynthesis.dict_functions as dictf
import botsynthesis.string_functions as strf
import botsynthesis.fitness_functions as fitf


# TODO finish tests
class TestMutations(unittest.TestCase):
    def test_mutate_codon_pass(self):
        # regular case
        initial_codon = 'aac'.upper()
        codon_table = CodonTable.unambiguous_dna_by_id[1]
        new_codon = mut.mutate_codon(initial_codon, codon_table)
        self.assertNotEqual(initial_codon, new_codon, 'mutate codon failed to change the codon when given the '
                                                      'opportunity')
        possible_codons = dictf.invert_dict(codon_table.forward_table)[codon_table.forward_table[initial_codon]]
        self.assertIn(new_codon, possible_codons, 'codon did not become a possible codon - wtf happened')
        # works with another translation dict
        codon_table = CodonTable.unambiguous_dna_by_id[2]
        new_codon = mut.mutate_codon(initial_codon, codon_table)
        self.assertNotEqual(initial_codon, new_codon, 'mutate codon failed to change the codon when given the '
                                                      'opportunity - using different table')
        possible_codons = dictf.invert_dict(codon_table.forward_table)[codon_table.forward_table[initial_codon]]
        self.assertIn(new_codon, possible_codons, 'codon did not become a possible codon - wtf happened - '
                                                  'using different table')
        # no other codons case
        initial_codon = 'tgg'.upper()
        codon_table = CodonTable.unambiguous_dna_by_id[1]
        new_codon = mut.mutate_codon(initial_codon, codon_table)
        self.assertEqual(initial_codon, new_codon, 'mutate codon failed to keep the codon when given only 1 option')

    def test_mutate_codon_fail(self):
        # smol codon
        initial_codon = 'ac'
        codon_table = 'not a codon tbale object'
        self.assertRaises(ValueError, mut.mutate_codon, initial_codon, codon_table)
        # bad codon usage
        initial_codon = 'acg'
        self.assertRaises(AttributeError, mut.mutate_codon, initial_codon, codon_table)

    def test_mutate_seq_pass(self):
        seq1 = 'actagctacggaattagcagaagcaatgctagc'.upper()
        # regular, make sure length doesnt change
        mutation_chance = 0.5
        codon_table = CodonTable.unambiguous_dna_by_id[1]
        mutant = mut.mutate_seq(seq1, mutation_chance, codon_table)
        self.assertEqual(len(mutant), len(seq1), 'lengths changed unexpectedly')
        # make sure the number of differences makes sense at 100% mutation
        mutation_chance = 1
        mutant = mut.mutate_seq(seq1, mutation_chance, codon_table)
        minimum_number_of_differences = 11
        maximum_number_of_differences = 20
        actual_number_of_differences = strf.find_num_differences(seq1, mutant)
        self.assertGreaterEqual(actual_number_of_differences, minimum_number_of_differences,
                                'too few mutations when 100%')
        self.assertLessEqual(actual_number_of_differences, maximum_number_of_differences,
                             'too many mutations when 100%')
        # make sure the number of differences makes sense at 0% mutation
        mutation_chance = 0
        mutant = mut.mutate_seq(seq1, mutation_chance, codon_table)
        actual_number_of_differences = strf.find_num_differences(seq1, mutant)
        self.assertEqual(actual_number_of_differences, 0, '0% mutations should result in the same string')

    def test_mutate_seq_fail(self):
        # sequence too small
        seq = ''
        mutation_chance = -1
        self.assertRaises(ValueError, mut.mutate_seq, seq, mutation_chance)
        # sequence not multiple of 3
        seq = 'reeeeeeeeeeeeeeeee'
        self.assertRaises(ValueError, mut.mutate_seq, seq, mutation_chance)
        # mutation chance too low
        seq = 'ACTGTA'
        self.assertRaises(ValueError, mut.mutate_seq, seq, mutation_chance)
        # mutation chance too high
        mutation_chance = 2
        self.assertRaises(ValueError, mut.mutate_seq, seq, mutation_chance)
        # sequence not DNA
        seq = 'reeeeeeee'
        mutation_chance = 0.5
        mut.mutate_seq(seq, mutation_chance)
        self.assertRaises(KeyError, mut.mutate_seq, seq, mutation_chance)

    def test_initialize_population_pass(self):
        # can it create a good population with default parameters
        seq = 'actagctacggaattagcagaagcaatgctagc'.upper()
        desired_size = 10
        population = mut.initialize_population(desired_population_size=desired_size, parent_sequence=seq,
                                               mutation_chance=0.7)
        self.assertEqual(len(population), desired_size, 'Failed to create a population diverse enough')
        # can it create a good population with tougher parameters?
        population = mut.initialize_population(desired_population_size=desired_size, parent_sequence=seq,
                                               mutation_chance=0.2)
        self.assertEqual(
            len(population), desired_size,
            'Failed to create a population diverse enough, using a low mutation chance'
        )

    def test_initialize_population_fail(self):
        # bad population
        seq = 'actagctacggaattagcagaagcaatgctagc'.upper()
        desired_size = 0
        self.assertRaises(ValueError, mut.initialize_population, desired_size, seq, 0.7)

    def test_tournament_selection_without_replacement_pass(self):
        key = 'key'
        population = {
            1: {key: 9099990.393}, 2: {key: 9099990.392}, 3: {key: 9099990.394}
        }
        winner = mut.tournament_selection_without_replacement(population, fitness_key_name=key)
        self.assertIn(winner, population.keys(), 'winner should be in the population')
        winner = mut.tournament_selection_without_replacement(population, n_ary=100, fitness_key_name=key)
        expected = 2
        self.assertEqual(expected, winner, 'winner should be 2, (2/3)^200 chance of not occurring')
        winner = mut.tournament_selection_without_replacement(population, n_ary=100, minima=False, fitness_key_name=key)
        expected = 3
        self.assertEqual(expected, winner, 'winner should be 3, (2/3)^200 chance of not occurring')

    def test_tournament_selection_without_replacement_fail(self):
        key = 'key'
        population = {
            1: {key: 9099990.393}, 2: {key: 9099990.392}, 3: {key: 9099990.394}
        }
        self.assertRaises(ValueError, mut.tournament_selection_without_replacement, population, 1)

    def test_generate_mating_pool_from_archive_pass(self):
        # can generate a mating pool from an archive
        key = 'key'
        archive = {1: {key: 1}, 2: {key: 2}, 3: {key: 3}, 4: {key: 4}, 5: {key: 5}}
        desired_size = 3
        mating_pool = mut.generate_mating_pool_from_archive(archive, desired_size, key)
        self.assertEqual(len(mating_pool), desired_size, 'mating pool did not reach desired size')
        self.assertTrue(mating_pool.items() <= archive.items(), 'mating pool is not a subset of the archive')
        # default key works
        key = fitf.__FITNESS_KEY__
        archive = {1: {key: 1}, 2: {key: 2}, 3: {key: 3}, 4: {key: 4}, 5: {key: 5}}
        mating_pool = mut.generate_mating_pool_from_archive(archive, desired_size)
        self.assertEqual(len(mating_pool), desired_size, 'mating pool did not reach desired size with default key')
        self.assertTrue(mating_pool.items() <= archive.items(),
                        'mating pool is not a subset of the archive with default key')

    def test_generate_mating_pool_from_archive_fail(self):
        # smol archive
        self.assertRaises(ValueError, mut.generate_mating_pool_from_archive, {}, 5)
        # smol mating pool
        archive = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
        self.assertRaises(ValueError, mut.generate_mating_pool_from_archive, archive, 0)
        # big mating pool
        self.assertRaises(ValueError, mut.generate_mating_pool_from_archive, archive, 100)

    def test_recombine_pass(self):
        # test that from 2 strings, the new one is created
        # new one is same aa
        s1 = 'actagctacggaattagcagaagcaatgctagc'
        s2 = 'acaagttatgggatctcgcgatcgaacgcaagt'
        aa1 = Seq(s1, IUPAC.unambiguous_dna).translate()
        aa2 = Seq(s2, IUPAC.unambiguous_dna).translate()
        self.assertEqual(aa2, aa1, 'The setup for recombine is bad, the test cannot be ran')
        for num in range(1, 20):
            r = mut.recombine_dna_sequence(s1, s2, num)
            raa = Seq(r, IUPAC.unambiguous_dna).translate()
            self.assertEqual(raa, aa1, 'Recombined sequence as AA is not the same AA sequence')
            for pos in range(len(r)):
                self.assertTrue(r[pos] == s1[pos] or r[pos] == s2[pos], 'the sequence is not composed of its parents')
            if num > 10:
                self.assertTrue(r != s1, 'at large numbers of sites, the string should theoretically nearly never be '
                                         'the same as the parent chance=1/11^10')
        # sending in the same sequence twice should result in the sequence being returned
        new_seq = mut.recombine_dna_sequence(s1, s1, 10)
        self.assertEqual(s1, new_seq, 'Failed to return parent sequence')

    def test_recombine_fail(self):
        s1 = 'actagctacggaattagcagaagcaatgctagc'
        s2 = 'acaagttatgggatctcgcgatcgaacgcaagt'
        # detect invalid number of sites
        self.assertRaises(ValueError, mut.recombine_dna_sequence, s1, s2, 0)
        # detect unequal sequences
        self.assertRaises(ValueError, mut.recombine_dna_sequence, s1, s2[:-2], 5)

    def test_generate_population_from_archive_pass(self):
        # can it run regularly with default vars
        s1 = 'actagctacggaattagcCGCagcaatgctagc'.upper()
        s2 = 'ACAAGCTACGGTATTTCTCGTAGTAATGCTTCA'
        s3 = 'ACTAGCTACGGCATCAGCCGCAGCAATGCTAGC'
        s4 = 'ACTAGTTACGGGATTTCTCGGAGCAACGCAAGC'
        s5 = 'ACGTCATACGGAATTAGCCGAAGTAATGCTAGC'
        s6 = 'ACCAGCTATGGAATTAGTCGAAGCAACGCGAGC'
        expected_aa = Seq(s1, IUPAC.unambiguous_dna).translate()
        desired_population_size = 15
        recomb_chance = mutation_chance = 0.2
        num_recombination_sites = mut.get_rec_sites_for_len(len(s1), recomb_chance)
        fit_key = fitf.__FITNESS_KEY__
        seq_key = mut.__SEQUENCE_KEY__
        archive = {1: {seq_key: s1, fit_key: 0.2},
                   2: {seq_key: s2, fit_key: 1.9},
                   3: {seq_key: s3, fit_key: 0.1},
                   4: {seq_key: s4, fit_key: 0.7},
                   5: {seq_key: s5, fit_key: 1.2},
                   6: {seq_key: s6, fit_key: 0.9}}

        desired_mating_pool_size = 4
        population = mut.generate_population_from_archive(
            archive, desired_mating_pool_size, num_recombination_sites,
            mutation_chance, desired_population_size)
        self.assertEqual(len(population), desired_population_size,
                         'Failed to create a population of correct size')
        for ind in population.values():
            if ind[seq_key] in [s1, s2, s3, s4, s5, s6]:
                logging.warning('new population member the same as an old member - undesirable behaviour')
            self.assertEqual(expected_aa, Seq(ind[seq_key], IUPAC.unambiguous_dna).translate(),
                             'Did not maintain Amino Acid integrity')
        # can it run regularly with custom vars
        codon_table = CodonTable.unambiguous_dna_by_id[2]
        expected_aa = Seq(s1, IUPAC.unambiguous_dna).translate(codon_table)
        fit_key = 'reeeee'
        new_archive = {1: {seq_key: s1, fit_key: 0.2},
                   2: {seq_key: s2, fit_key: 1.9},
                   3: {seq_key: s3, fit_key: 0.1},
                   4: {seq_key: s4, fit_key: 0.7},
                   5: {seq_key: s5, fit_key: 1.2},
                   6: {seq_key: s6, fit_key: 0.9}}
        population = mut.generate_population_from_archive(
            new_archive, desired_mating_pool_size, num_recombination_sites,
            mutation_chance, desired_population_size, fit_key, codon_table)
        self.assertEqual(len(population), desired_population_size,
                         'Failed to create a population of correct size - custom')
        for ind in population.values():
            self.assertEqual(expected_aa, Seq(ind[seq_key], IUPAC.unambiguous_dna).translate(codon_table),
                             'Did not maintain Amino Acid integrity - custom')

    def test_get_rec_sites_for_len_pass(self):
        pass

    def test_get_rec_sites_for_len_fail(self):
        pass


if __name__ == '__main__':
    unittest.main()

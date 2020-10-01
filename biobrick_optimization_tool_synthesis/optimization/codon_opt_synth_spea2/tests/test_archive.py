import unittest

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.test_parameters import algorithm_params
import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions as fit_func
import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.archive as arch


class TestArchive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.correct_sort_key = '1'
        cls.not_exist_sort_key = '3'
        seq1 = {cls.correct_sort_key: 5}
        seq2 = {cls.correct_sort_key: 4}
        seq3 = {cls.correct_sort_key: 3}
        seq4 = {cls.correct_sort_key: 3}
        seq5 = {cls.correct_sort_key: 2.5}
        seq6 = {cls.correct_sort_key: 1}

        cls.empty = {}
        cls.pop_1 = {'1': seq1}
        cls.pop_5 = {'1': seq1, '2': seq2, '3': seq3, '4': seq4, '5': seq5}
        cls.pop_6 = {'1': seq1, '2': seq2, '3': seq3, '4': seq4, '5': seq5, '6': seq6}
        cls.pop_r1 = {'6': seq6}
        cls.pop_r5 = {'6': seq6, '5': seq5, '4': seq4, '3': seq3, '2': seq2}
        cls.pop_r6 = {'6': seq6, '5': seq5, '4': seq4, '3': seq3, '2': seq2, '1': seq1}

    # build archive black box
    # empty population, population of 1, population of many (5)
    # archive size of -1, 0, 1, between 1 and size of population,(size of population or greater than size of population)
    # sort key correct, not existing
    # if nothing goes badly, an archive (subset of the population) with size equal to archive size will be returned
    # else descriptive error message
    def test_build_archive_not_existing(self):
        pops = [self.pop_5, self.pop_1]
        for pop in pops:
            self.assertRaises(KeyError, arch.build_archive(pop, len(pop) + 1, self.not_exist_sort_key))
            self.assertRaises(KeyError, arch.build_archive(pop, len(pop), self.not_exist_sort_key))
            self.assertRaises(KeyError, arch.build_archive(pop, len(pop) - 1, self.not_exist_sort_key))
            self.assertRaises(KeyError, arch.build_archive(pop, 1, self.not_exist_sort_key))
            self.assertRaises(KeyError, arch.build_archive(pop, 0, self.not_exist_sort_key))
            self.assertRaises(KeyError, arch.build_archive(pop, -1, self.not_exist_sort_key))

    def test_build_archive_pop5_arch_size_gt_pop(self):
        archive = arch.build_archive(self.pop_5, len(self.pop_5) + 1, self.correct_sort_key)
        if len(archive) != len(self.pop_5):
            self.fail('Incorrect archive size')
        for x in range(1, 6):
            if str(x) not in archive:
                self.fail('Failed to find expected key in result archive')

    def test_build_archive_pop5_arch_size_eq_pop(self):
        archive = arch.build_archive(self.pop_5, len(self.pop_5), self.correct_sort_key)
        if len(archive) != len(self.pop_5):
            self.fail('Incorrect archive size')
        for x in range(1, 6):
            if str(x) not in archive:
                self.fail('Failed to find expected key in result archive')

    def test_build_archive_pop5_arch_size_lt_pop(self):
        archive = arch.build_archive(self.pop_5, len(self.pop_5) - 1, self.correct_sort_key)
        if len(archive) != len(self.pop_5) - 1:
            self.fail('Incorrect archive size')
        for x in range(2, 6):
            if str(x) not in archive:
                self.fail('Failed to find expected key in result archive')

    def test_build_archive_pop5_arch_size_1(self):
        archive = arch.build_archive(self.pop_5, 1, self.correct_sort_key)
        if len(archive) != 1:
            self.fail('Incorrect archive size')
        if '5' not in archive:
            self.fail('Failed to find expected key in result archive')

    def test_build_archive_pop5_arch_size_0(self):
        archive = arch.build_archive(self.pop_5, 0, self.correct_sort_key)
        if len(archive) != 0:
            self.fail('Incorrect archive size')

    def test_build_archive_pop5_arch_size_neg_1(self):
        self.assertRaises(ValueError, arch.build_archive(self.pop_5, -1, self.correct_sort_key))

    def test_build_archive_pop1_arch_size_gt_pop(self):
        archive = arch.build_archive(self.pop_1, len(self.pop_1) + 1, self.correct_sort_key)
        if len(archive) != len(self.pop_1):
            self.fail('Incorrect archive size')
        if '1' not in archive:
            self.fail('Failed to find expected key in result archive')

    def test_build_archive_pop1_arch_size_eq_pop(self):
        archive = arch.build_archive(self.pop_1, len(self.pop_1), self.correct_sort_key)
        if len(archive) != len(self.pop_1):
            self.fail('Incorrect archive size')
        if '1' not in archive:
            self.fail('Failed to find expected key in result archive')

    def test_build_archive_pop1_arch_size_0(self):
        archive = arch.build_archive(self.pop_1, 0, self.correct_sort_key)
        if len(archive) != 0:
            self.fail('Incorrect archive size')

    def test_build_archive_pop1_arch_size_neg_1(self):
        self.assertRaises(ValueError, arch.build_archive(self.pop_1, -1, self.correct_sort_key))

    def test_build_archive_empty_arch_size_1(self):
        archive = arch.build_archive(self.empty, 1, self.correct_sort_key)
        if len(archive) != 0:
            self.fail('Expected empty archive')

    def test_build_archive_empty_arch_size_0(self):
        archive = arch.build_archive(self.empty, 0, self.correct_sort_key)
        if len(archive) != 0:
            self.fail('Expected empty archive')

    def test_build_archive_empty_arch_size_neg_1(self):
        self.assertRaises(ValueError, arch.build_archive(self.empty, -1, self.correct_sort_key))

    # add dominated to archive black box
    # empty archive, archive size 1, archive smaller than population, archive == population, archive larger than population
    # empty population, population of 1, population smaller than archive, archive == population, population larger than archive
    # archive size of -1, 0, 1, 1<size<pop size <= arch size, 1<size<arch size <= pop size, size = archive, size = population, size > archive, size > population
    # sort key correct, not existing
    # if nothing goes badly, the size of the archive will INCREASE to the archive size
    # else an error should pop up
    # def test_add_dominated_to_archive_arch5_pop1_arch_size_gt_pop

    # truncate archive
    # empty archive, population of 1, population of many (5)
    # archive size of -1, 0, 1, between 1 and size of archive, size of archive, greater than size of archive
    # sort key correct, not existing
    # if nothing goes badly, an archive (subset of the population) with size equal to archive size will be returned


if __name__ == '__main__':
    unittest.main()

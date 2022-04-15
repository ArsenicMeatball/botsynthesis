import unittest
import copy
import src.archive as arch


class TestArchive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.correct_sort_key = '1'
        cls.not_exist_sort_key = 'tronbikesupreme'
        cls.seq1 = {cls.correct_sort_key: 5}
        cls.seq2 = {cls.correct_sort_key: 4}
        cls.seq3 = {cls.correct_sort_key: 3}
        cls.seq4 = {cls.correct_sort_key: 3}
        cls.seq5 = {cls.correct_sort_key: 2.5}
        cls.seq6 = {cls.correct_sort_key: 1}

        cls.empty = {}
        cls.pop_1 = {'1': cls.seq1}
        cls.pop_5 = {'1': cls.seq1, '2': cls.seq2, '3': cls.seq3, '4': cls.seq4, '5': cls.seq5}
        cls.pop_6 = {'1': cls.seq1, '2': cls.seq2, '3': cls.seq3, '4': cls.seq4, '5': cls.seq5, '6': cls.seq6}
        cls.pop_r1 = {'6': cls.seq6}
        cls.pop_r5 = {'6': cls.seq6, '5': cls.seq5, '4': cls.seq4, '3': cls.seq3, '2': cls.seq2}
        cls.pop_r6 = {'6': cls.seq6, '5': cls.seq5, '4': cls.seq4, '3': cls.seq3, '2': cls.seq2, '1': cls.seq1}

    def test_build_archive_not_existing(self):
        pops = [self.pop_5, self.pop_1]
        for pop in pops:
            self.assertRaises(KeyError, arch.build_archive, pop, len(pop) - 1, self.not_exist_sort_key)
            self.assertRaises(KeyError, arch.build_archive, pop, 0, self.not_exist_sort_key)
            self.assertRaises(KeyError, arch.build_archive, pop, -1, self.not_exist_sort_key)

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
        self.assertRaises(IndexError, arch.build_archive, self.pop_5, -1, self.correct_sort_key)

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
        self.assertRaises(IndexError, arch.build_archive, self.pop_1, -1, self.correct_sort_key)

    def test_build_archive_empty_arch_size_1(self):
        archive = arch.build_archive(self.empty, 1, self.correct_sort_key)
        if len(archive) != 0:
            self.fail('Expected empty archive')

    def test_build_archive_empty_arch_size_0(self):
        archive = arch.build_archive(self.empty, 0, self.correct_sort_key)
        if len(archive) != 0:
            self.fail('Expected empty archive')

    def test_build_archive_empty_arch_size_neg_1(self):
        self.assertRaises(IndexError, arch.build_archive, self.empty, -1, self.correct_sort_key)

    def test_add_dominated_to_archive_pass(self):
        # should add the best dominated members from the population into the archive
        expected = {'5': self.seq5, '6': self.seq6, '4': self.seq4, '3': self.seq3}
        actual = {'5': self.seq5, '6': self.seq6}
        pop = copy.deepcopy(self.pop_6)
        desired_size = 4
        arch.add_dominated_to_archive(actual, pop, desired_size, self.correct_sort_key)
        self.assertEqual(expected, actual, 'Partial-Did not result in the same dict')
        # same thing but with an empty dict
        actual = {}
        arch.add_dominated_to_archive(actual, pop, desired_size, self.correct_sort_key)
        self.assertEqual(expected, actual, 'Empty-Did not result in the same dict')

    def test_add_dominated_to_archive_fail(self):
        # population too small
        archive = population = {}
        desired_size = 10000000
        self.assertRaises(ValueError, arch.add_dominated_to_archive, archive, population, desired_size, self.correct_sort_key)
        # wrong key
        self.assertRaises(KeyError, arch.add_dominated_to_archive, archive, self.pop_6, 1, self.not_exist_sort_key)

    def test_truncate_archive_pass(self):
        # should remove worst (highest scores) until correct size
        actual = copy.deepcopy(self.pop_6)
        desired_size = 0
        expected = {}
        arch.truncate_archive(actual, desired_size, self.correct_sort_key)
        self.assertEqual(expected, actual, 'Can not truncate archive down to 0')
        actual = copy.deepcopy(self.pop_6)
        desired_size = 4
        expected = {'5': self.seq5, '6': self.seq6, '4': self.seq4, '3': self.seq3}
        arch.truncate_archive(actual, desired_size, self.correct_sort_key)
        self.assertEqual(expected, actual, 'Can not truncate archive down to regular case')

    def test_truncate_archive_fail(self):
        # archive smaller than desired size
        self.assertRaises(ValueError, arch.truncate_archive, {}, 10000, self.correct_sort_key)
        # wrong key
        self.assertRaises(KeyError, arch.truncate_archive, self.pop_6, 0, self.not_exist_sort_key)


if __name__ == '__main__':
    unittest.main()

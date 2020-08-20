import unittest

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.string_functions import *
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.test_parameters import algorithm_params


class TestFindingRepeats(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sequences = ['123123123456456456789789789',
                         'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
                         str(algorithm_params['codon_opt_seq'])
                         ]
        cls.low_repeat_sizes = [n for n in range(-1, 1)]
        cls.repeat_sizes = [n for n in range(2, 11)]

    def test_overlapping_repeats(self):
        for sequence in self.sequences:
            for n in self.repeat_sizes:
                result1 = find_number_of_overlapping_repeats(string=sequence, min_repeat_size=n)
                repeats = find_repeats(string=sequence, min_repeat_size=n, overlapping=True)
                result2 = get_number_of_repeats_from_dict(repeats)
                self.assertEqual(result1, result2, "number of overlapping repeats must be the same")

    def test_non_overlapping_repeats(self):
        for sequence in self.sequences:
            for n in self.repeat_sizes:
                result1 = find_number_of_non_overlapping_repeats(string=sequence, min_repeat_size=n)
                repeats = find_repeats(string=sequence, min_repeat_size=n, overlapping=False)
                result2 = get_number_of_repeats_from_dict(repeats)
                self.assertEqual(result1, result2, "number of overlapping repeats must be the same")

    def test_small_repeat_size(self):
        for sequence in self.sequences:
            for n in self.low_repeat_sizes:
                self.assertRaises(ArithmeticError,
                                  find_number_of_non_overlapping_repeats, sequence, n)
                self.assertRaises(ArithmeticError,
                                  find_number_of_overlapping_repeats, sequence, n)
                self.assertRaises(ArithmeticError,
                                  find_repeats, sequence, n)

    def test_large_min_repeat_size(self):
        increase = 5
        for sequence in self.sequences:
            result = find_number_of_non_overlapping_repeats(string=sequence, min_repeat_size=len(sequence) + increase)
            self.assertEqual(result, 0, "a repeat size larger than sequence should yield 0")
            result = find_number_of_overlapping_repeats(string=sequence, min_repeat_size=len(sequence) + increase)
            self.assertEqual(result, 0, "a repeat size larger than sequence should yield 0")
            result = find_repeats(string=sequence, min_repeat_size=len(sequence) + increase)
            self.assertEqual(result, {}, "a repeat size larger than sequence should yield an empty dict")


if __name__ == '__main__':
    unittest.main()

import unittest

from botsynthesis.string_functions import *
from botsynthesis.all_algorithm_parameters import algorithm_params


class TestStringFunctions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sequences = [
            "123123123456456456789789789",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            str(algorithm_params["codon opt seq"]),
        ]
        cls.low_repeat_sizes = [n for n in range(-1, 1)]
        cls.repeat_sizes = [n for n in range(2, 11)]

    def test_overlapping_repeats(self):
        for sequence in self.sequences:
            for n in self.repeat_sizes:
                result1 = find_number_of_overlapping_repeats(
                    string=sequence, min_repeat_size=n
                )
                repeats = find_repeats(
                    string=sequence, min_repeat_size=n, overlapping=True
                )
                result2 = get_number_of_repeats_from_repeats_dict(repeats)
                self.assertEqual(
                    result1,
                    result2,
                    "number of overlapping repeats must be the same",
                )

    def test_non_overlapping_repeats(self):
        for sequence in self.sequences:
            for n in self.repeat_sizes:
                result1 = find_number_of_non_overlapping_repeats(
                    string=sequence, min_repeat_size=n
                )
                repeats = find_repeats(
                    string=sequence, min_repeat_size=n, overlapping=False
                )
                result2 = get_number_of_repeats_from_repeats_dict(repeats)
                self.assertEqual(
                    result1,
                    result2,
                    "number of overlapping repeats must be the same",
                )

    def test_small_repeat_size(self):
        for sequence in self.sequences:
            for n in self.low_repeat_sizes:
                self.assertRaises(
                    ArithmeticError,
                    find_number_of_non_overlapping_repeats,
                    sequence,
                    n,
                )
                self.assertRaises(
                    ArithmeticError,
                    find_number_of_overlapping_repeats,
                    sequence,
                    n,
                )
                self.assertRaises(ArithmeticError, find_repeats, sequence, n)

    def test_large_min_repeat_size(self):
        increase = 5
        for sequence in self.sequences:
            result = find_number_of_non_overlapping_repeats(
                string=sequence, min_repeat_size=len(sequence) + increase
            )
            self.assertEqual(
                result, 0, "a repeat size larger than sequence should yield 0"
            )
            result = find_number_of_overlapping_repeats(
                string=sequence, min_repeat_size=len(sequence) + increase
            )
            self.assertEqual(
                result, 0, "a repeat size larger than sequence should yield 0"
            )
            result = find_repeats(
                string=sequence, min_repeat_size=len(sequence) + increase
            )
            self.assertEqual(
                result,
                {},
                "a repeat size larger than sequence should yield an empty dict",
            )

    def test_find_separated_palindromes_pass(self):
        palindrome = "1QWERTYUIOP788978POIUYTREWQ1ASDFGHJKL:URT:LKJHGFDSA"
        expected = {
            6: [1, "788978", 26],
            3: [28, "URT", 50],
            8: [0, "P788978P", 27],
        }
        actual = find_separated_palindromes(palindrome)
        self.assertEqual(expected, actual)
        palindrome = ""
        expected = {}
        actual = find_separated_palindromes(palindrome)
        self.assertEqual(expected, actual)

    def test_find_num_differences_pass(self):
        string1 = "asdfghjklqwertyuiop"
        string2 = "aadfhhjklewervbuiop"
        expected = 5
        actual = find_num_differences(string1, string2)
        self.assertEqual(expected, actual)

    def test_get_number_of_repeats_from_dict_pass(self):
        self.assertEqual(True, False)


if __name__ == "__main__":
    unittest.main()

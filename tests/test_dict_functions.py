import unittest

from src.dict_functions import *
from Bio.Data import CodonTable


# TODO: Test fail paths
class TestDictFunctions(unittest.TestCase):
    def setUp(self) -> None:
        self.d1 = {1: ['123', '23'], 2: ['mm', 3, 5, 4, 6, 3], 3: [], 5: [1, 2, 3, 4, 5]}
        self.d2 = {1: ['23'], 2: ['mm', 3, 5, 4, 6, 3], 4: [], 5: [1, 2]}

    def test_get_differences_dict_list(self):
        expected = {1: ['123'], 3: [], 4: [], 5: [3, 4, 5]}
        actual = get_differences_dict_list(self.d1, self.d2)
        self.assertEqual(expected, actual)

    def test_get_number_of_differences_dict_list(self):
        expected = 4
        actual = get_number_of_differences_dict_list(self.d1, self.d2)
        self.assertEqual(expected, actual)

    def test_find_keys_of_minimum_value(self):
        d = {1: 78787787787, 3: 3, 4: 7, 'aohd': 92000, 'ajfp': 4, 'we': 5, 'q': 6, 'eoqafi': 3, 'o0aihf': 4,
             'ofajiw': 3}
        expected = [3, 'eoqafi', 'ofajiw']
        actual = find_keys_of_minimum_value(d)
        self.assertEqual(expected, actual)

    def test_sort_dict_by_value_get_list_of_keys(self):
        d = {1: 78787787787, 3: 3, 4: 7, 'aohd': 92000, 'ajfp': 4, 'we': 5, 'q': 6}
        expected = [3, 'ajfp', 'we', 'q', 4, 'aohd', 1]
        actual = sort_dict_by_value_get_list_of_keys(d)
        self.assertEqual(expected, actual)

    def test_sort_dict_by_value(self):
        d = {1: 78787787787, 3: 3, 4: 7, 'aohd': 92000, 'ajfp': 4, 'we': 5, 'q': 6}
        expected = [3, 4, 5, 6, 7, 92000, 78787787787]
        actual = sort_dict_by_value(d)
        self.assertEqual(expected, actual)

    def test_invert_dict(self):
        table = CodonTable.generic_by_id[1].forward_table
        print(table)
        actual = invert_dict(table)

        expected = {}
        for codon, amino in table.items():
            if amino not in expected:
                expected[amino] = set()
            expected[amino].add(codon)

        self.assertEqual(actual, expected, 'could not invert correctly')

if __name__ == '__main__':
    unittest.main()

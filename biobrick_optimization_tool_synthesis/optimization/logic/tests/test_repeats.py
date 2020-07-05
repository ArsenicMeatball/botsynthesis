import unittest

from biobrick_optimization_tool_synthesis.optimization.logic.string_manipulation import *
from biobrick_optimization_tool_synthesis.optimization.logic.test_parameters import algorithm_params


class TestFindingRepeats(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sequence1 = '123123123456456456789789789'
        sequence2 = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        sequence3 = str(algorithm_params['codon_opt_seq'])
    def test_overlapping_repeats(self):
        result1 = find_number_of_overlapping_repeats(string=spec_seq, min_repeat_size=n)
    def test_non_overlapping_repeats(self):

        self.assertEqual(True, False)




    for n in range(3, 11):
        result1 =
        print(result1)
        result2 = find_number_of_non_overlapping_repeats(string=spec_seq, min_repeat_size=n)
        print(result2)
        result3 = find_repeats(string=spec_seq, min_repeat_size=n,
                               overlapping=True)
        print(sorted(result3.keys()))
        result31 = 0
        for v in result3.values():
            result31 += len(v)
        print(result31 - len(result3))
        result4 = find_repeats(string=spec_seq, min_repeat_size=n,
                               overlapping=False)
        print(result4)
        result41 = 0
        for v in result4.values():
            result41 += len(v)
        print(result41 - len(result4))



if __name__ == '__main__':
    unittest.main()

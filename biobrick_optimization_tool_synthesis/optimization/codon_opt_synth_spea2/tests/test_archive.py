import unittest

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.test_parameters import algorithm_params
import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions as fit_func

class TestArchive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass
    # add dominated to archive black box
# empty archive, archive size 1, archive smaller than population, archive == population, archive larger than population
# empty population, population of 1, population smaller than archive, archive == population, population larger than archive
# archive size of -1, 0, 1, 1<size<size of population, 1<size<size of archive, size = archive, size = population, size > archive, size > population
# if nothing goes badly, the size of the archive will INCREASE to the archive size
# else an error should pop up

    # truncate archive

    # build archive black box
# empty population, population of 1, population of many (5)
# archive size of -1, 0, 1, between 1 and size of population, size of population, greater than size of population
# if nothing goes badly, an archive (subset of the population) with size equal to archive size will be returned
# else descriptive error message
    def test_build_archive_(self):
        pass


if __name__ == '__main__':
    unittest.main()

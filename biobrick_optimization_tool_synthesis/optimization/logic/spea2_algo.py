import multiprocessing

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from biobrick_optimization_tool_synthesis.optimization.logic import interpret_sequence, test_parameters
from biobrick_optimization_tool_synthesis.optimization.logic.codon_optimization_host import \
    determine_ideal_codon_optimized_sequence
from biobrick_optimization_tool_synthesis.optimization.logic.fitness_functions import fitness_evals
from biobrick_optimization_tool_synthesis.optimization.logic.mutations import initialize_population
from biobrick_optimization_tool_synthesis.optimization.logic.set_codon_table import fetch_codon_table
from biobrick_optimization_tool_synthesis.optimization.logic.test_parameters import algorithm_params


def spea2_main_loop(params: dict):
    # create population
    population = initialize_population(params)
    # start loop
    generation = 0
    is_converged = False
    while generation < params['generations']:
        # check end conditions
        if is_converged:
            break
        # multiprocessing check fitness
        processes = []
        for func in fitness_evals:
            process = multiprocessing.Process(
                target=func,
                args=tuple([(k, v) for k, v in params.items()])
            )
            processes.append(process)
            process.start()

    # domination
    # fitness


if __name__ == '__main__':
    spea2_main_loop(algorithm_params)

import multiprocessing as mp

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from biobrick_optimization_tool_synthesis.optimization.logic import cleaning, test_parameters
from biobrick_optimization_tool_synthesis.optimization.logic.codon_optimization_host import \
    determine_ideal_codon_optimized_sequence
from biobrick_optimization_tool_synthesis.optimization.logic.fitness_functions import fitness_evals
from biobrick_optimization_tool_synthesis.optimization.logic.mutations import initialize_population
from biobrick_optimization_tool_synthesis.optimization.logic.set_codon_table import fetch_codon_table
from biobrick_optimization_tool_synthesis.optimization.logic.test_parameters import algorithm_params


def spea2_main_loop(params: dict):
    # create population
    params['population'] = initialize_population(params)
    # start loop
    generation = 0
    is_converged = False
    out_q = mp.Queue()
    for generation in range(params['generations']):
        # check end conditions
        if is_converged:
            break
        # multiprocessing check fitness
        processes = []
        for func in fitness_evals:
            process = mp.Process(
                target=func,
                args=(params, out_q)
            )
            processes.append(process)
            process.start()

        # retrieve data
        result = {}
        for i in range(len(processes)):
            result.update(out_q.get())

        # kill processes
        for p in processes:
            p.join()

        # update the main data store
        for eval_type, seq_n_score in result.items():
            for seq, score in seq_n_score.items():
                params['population'][seq][eval_type] = score

    for k, v in params['population'].items():
        print(k, v)
    # domination
    # fitness


if __name__ == '__main__':
    spea2_main_loop(algorithm_params)

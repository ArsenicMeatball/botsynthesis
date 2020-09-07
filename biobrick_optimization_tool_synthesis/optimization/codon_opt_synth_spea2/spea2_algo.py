import logging
import multiprocessing as mp

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.archive import build_archive
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions import fitness_evals, \
    calculate_fitness
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.mutations import initialize_population, \
    generate_population_from_archive
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.test_parameters import algorithm_params
# /home/arsenic/.cache/JetBrains/PyCharm2020.2/snapshots/BOTS_development1.pstat

def spea2_main_loop(params: dict) -> dict:
    """
    Main loop of spea2 algorithm
    https://pdfs.semanticscholar.org/6672/8d01f9ebd0446ab346a855a44d2b138fd82d.pdf
    Considered one of, if not the best Multi-objective optimization algorithm for highly dimensional problems.
    :param params: dict containing everything necessary to run the program
    :return: dict containing the final archive
    """
    # create population
    params['population'] = initialize_population(params)
    archive_sequences_per_gen = {}
    # main loop
    is_converged = False
    out_q = mp.Queue()
    for generation in range(params['generations']):
        # check end conditions
        if is_converged:
            break
        # multiprocessing check fitness
        processes = []
        for fitness_eval in fitness_evals:
            process = mp.Process(
                target=fitness_eval,
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
        for eval_type, id_score in result.items():
            for seq_id, score in id_score.items():
                params['population'][seq_id][eval_type] = score

        # calculate fitness
        calculate_fitness(params['population'])
        # build archive
        params['archive'] = build_archive(params['population'], params['archive size'])
        # generate population of next generation
        params['population'] = generate_population_from_archive(params)
        print(generation)
    return params['archive']


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)
    result = spea2_main_loop(algorithm_params)

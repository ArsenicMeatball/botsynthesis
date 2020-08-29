import multiprocessing as mp

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions import fitness_evals, \
    calculate_fitness
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.mutations import initialize_population
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.test_parameters import algorithm_params


def spea2_main_loop(params: dict):
    # create population
    params['population'], sequences_set = initialize_population(params)
    # main loop
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
        for eval_type, id_score in result.items():
            for seq_id, score in id_score.items():
                params['population'][seq_id][eval_type] = score

        # calculate fitness
        calculate_fitness(params['population'])

        for k1, kv in params['population'].items():
            print('{0} :'.format(k1))
            for k2, v in kv.items():
                print(' {0} : {1}'.format(k2, v))


if __name__ == '__main__':
    spea2_main_loop(algorithm_params)

import logging
import multiprocessing as mp
import pickle

# /home/arsenic/.cache/JetBrains/PyCharm2020.2/snapshots/BOTS_development1.pstat
from botsynthesis.utils.constants import FITNESS_KEY
from tests.all_algorithm_parameters import algorithm_params
from botsynthesis.spea2.archive import build_archive
from botsynthesis.spea2.mutations import (
    initialize_population,
    generate_population_from_archive,
)
import botsynthesis.spea2.fitness_functions as fit_func


def spea2_main_loop(params: dict) -> dict:
    """
    Main loop of spea2 algorithm
    https://pdfs.semanticscholar.org/6672/8d01f9ebd0446ab346a855a44d2b138fd82d.pdf
    Considered one of, if not the best Multi-objective optimization algorithm
    for highly dimensional problems.
    :param params: dict containing everything necessary to run the program
    :return: dict containing the final archive
    """
    # create population
    params["population"] = initialize_population(
        params["population size"],
        params["codon opt seq"],
        params["mutation %"],
        params["codon_table"],
    )
    # archive_sequences_per_gen = {}
    # main loop
    is_converged = False
    out_q = mp.Queue()
    for generation in range(params["generations"]):
        # check end conditions
        if is_converged:
            break
        # multiprocessing check fitness
        processes = []
        for fitness_eval in fit_func.fitness_evals:
            process = mp.Process(target=fitness_eval, args=(params, out_q))
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
                params["population"][seq_id][eval_type] = score

        # calculate fitness
        fit_func.calculate_fitness(params["population"])
        # build archive
        params["archive"] = build_archive(
            params["population"],
            params["archive size"],
            FITNESS_KEY,
        )
        # generate population of next generation
        params["population"] = generate_population_from_archive(
            params["archive"],
            params["mating pool size"],
            params["num crossover sites"],
            params["mutation %"],
            params["population size"],
            "fitness",
            params["codon_table"],
        )
        print(generation)
    return params["archive"]


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)
    result = spea2_main_loop(algorithm_params)
    print(result)
    with open("../../playground/output.pkl", "wb") as p:
        pickle.dump(result, p)

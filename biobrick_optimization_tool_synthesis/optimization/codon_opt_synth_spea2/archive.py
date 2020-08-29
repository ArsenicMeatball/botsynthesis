import logging

import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions as fit_func


def build_archive(population: dict, archive_size: int) -> dict:
    """
    best members of the population go into the archive
    :param archive_size: the size of the archive
    :param population: the union between the population and the "archive"
    :return: archive dict of best members for next generation
    """
    if len(population) <= archive_size:
        logging.warning(
            'Archive larger ({0} items) than supplied population ({1} items)!'.format(archive_size, len(population))
        )
        return population
    archive = {}
    logging.info('Adding all non_dominated sequences into archive')
    for seq_id in population.keys():
        if population[seq_id][fit_func.__FITNESS_KEY__] < 1:
            archive[seq_id] = population[seq_id]
    # get archive to correct size
    if len(archive) > archive_size:
        truncate_archive(archive, archive_size)
    elif len(archive) < archive_size:
        add_dominated_to_archive(archive, population, archive_size)
    return archive


def add_dominated_to_archive(archive, population, archive_size):
    """
    get best members into the archive until it is full
    - expects population size > archive size
    :param archive: the archive we are building
    :param population: the population we are pulling members from
    :param archive_size:
    :return: None, just updates archive
    """
    sorted_population = sorted(population.keys(), key=lambda x: population[x][fit_func.__FITNESS_KEY__])
    print(sorted_population)
    while len(archive) < archive_size:
        break


def truncate_archive(archive, archive_size):
    pass


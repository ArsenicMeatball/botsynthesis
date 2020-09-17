import logging


def build_archive(population: dict, archive_size: int, sort_key: str, non_dominated_upper_limit: int = 1) -> dict:
    """
    best members of the population go into the archive
    :param archive_size: the size of the archive
    :param population: the union between the population and the "archive"
    :param sort_key: key by which to sort the population
    :param non_dominated_upper_limit: default 1, all scores below this go into archive automatically,
        but if size too large, some get truncated
    :return: archive dict of best members for next generation
    """
    logging.info('building archive')
    if len(population) <= archive_size:
        logging.warning(
            'Archive larger ({0} items) than supplied population ({1} items)!'.format(archive_size, len(population))
        )
        return population
    archive = {}
    logging.info('Adding all non_dominated sequences into archive')
    for seq_id in population.keys():
        if population[seq_id][sort_key] < non_dominated_upper_limit:
            archive[seq_id] = population[seq_id]
    # get archive to correct size
    if len(archive) > archive_size:
        truncate_archive(archive, archive_size, sort_key)
    elif len(archive) < archive_size:
        add_dominated_to_archive(archive, population, archive_size, sort_key)
    return archive


def add_dominated_to_archive(archive: dict, population: dict, archive_size: int, sort_key: str):
    """
    get best members into the archive until the right size is achieved
    - expects population size > archive size
    :param archive: the archive we are building
    :param population: the population we are pulling members from
    :param archive_size: the final size we are trying to get to
    :param sort_key: key by which to sort the population
    :return: None, just updates archive
    """
    logging.info('adding best dominated members from population to archive')
    sorted_keys = sorted(population.keys(), key=lambda x: population[x][sort_key])
    while len(archive) < archive_size:
        next_item = sorted_keys[len(archive)]
        archive[next_item] = population[next_item]


def truncate_archive(archive: dict, archive_size: int, sort_key: str):
    """
    removes worst members from archive until the right size is achieved
    - expects population size > archive size
    :param archive: the archive we are truncating
    :param archive_size: the final size we are trying to get to
    :param sort_key: key by which to sort the population
    :return: None, just updates archive
    """
    logging.info('truncating worst non-dominated members from archive')
    sorted_keys = sorted(archive.keys(), key=lambda x: archive[x][sort_key])
    while len(archive) > archive_size:
        archive.pop(sorted_keys[-1])
        sorted_keys.pop(-1)

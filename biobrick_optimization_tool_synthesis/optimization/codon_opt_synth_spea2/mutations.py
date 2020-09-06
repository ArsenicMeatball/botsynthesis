import copy
import logging
import random
import uuid

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions as fit_func

__SEQUENCE_KEY__ = 'sequence'


def mutate_codon(sequence: str, codon_usage: dict, idx=0, translation_dict: dict = None):
    """

    :param sequence:
    :param idx:
    :param codon_usage:
    :param translation_dict:
    :return:
    """

    codon = str(sequence[idx:idx + 3])
    amino = Seq(codon, IUPAC.unambiguous_dna).translate()
    if translation_dict is not None:
        amino = Seq(codon, IUPAC.unambiguous_dna).translate(translation_dict)
    if len(codon_usage[amino]) < 1:
        raise Exception(
            'amino acid recovered is not compatible with codon usage table given. \n'
            'please provide a translation table compatible with the codon_usage table'
        )

    new_codon = random.choice(list(codon_usage[amino].keys()))
    if len(codon_usage[amino]) > 1:
        while new_codon == codon:  # if the same try again for higher chance of actual mutation
            new_codon = random.choice(list(codon_usage[amino].keys()))
            # print('New codon: ', new_codon)
    return new_codon


def mutate_seq(sequence: str, param: dict) -> str:
    """
        mutate a DNA Seq based on the mutation probability, returns a different DNA Seq
    :param sequence:
    :param param:
    :return:
    """
    new_sequence = ''
    for idx in range(0, len(sequence), 3):
        # either add the og codon or mutated boy
        new_sequence += (
            str(sequence[idx:idx + 3])
            if random.randint(1, 100) > param['mutation %']
            else str(mutate_codon(sequence, param['codon usage'], idx))
        )
        # print('Current sequence: ', new_sequence)

    return new_sequence


def initialize_population(algorithm_parameters: dict) -> tuple:
    """
        create population of mutants based on the codon optimized sequence
    :param algorithm_parameters: dict containing the parameters for the SPEA2 algorithm
    :return: tuple(dict representing the population, set representing the sequences in the population)
    """
    population = {}
    sequences = set()
    attempts = 0
    while len(population) < algorithm_parameters['population size'] and \
            attempts < algorithm_parameters['max init population attempts']:
        seq = mutate_seq(algorithm_parameters['codon opt seq'], algorithm_parameters)
        if seq not in sequences:
            sequences.add(seq)
            seq_id = uuid.uuid4()
            population[seq_id] = {__SEQUENCE_KEY__: seq}
        else:
            attempts += 1
            logging.info('failed to create a new sequence')
    return population, sequences


def tournament_selection_without_replacement(population: dict, n_ary: int = 2):
    """
    Get winner of a tournament (key of winner) based on the fitness value
    :param n_ary: number of individuals to conduct tournament on, default = binary
    :param population: all the individuals we are selecting from
    :return: the key of the winner of tournament
    """
    if n_ary < 2:
        raise ValueError('tournament selection must be binary or larger, currently {0}'.format(n_ary))
    population_keys = list(population.keys())
    key_of_minimum = population_keys[random.randint(0, len(population_keys))]
    for _ in range(1, n_ary):
        idx = random.randint(0, len(population))
        # TODO: make if statement more general, not just fitness value
        if population[population_keys[idx]][fit_func.__FITNESS_KEY__] < \
                population[key_of_minimum][fit_func.__FITNESS_KEY__]:
            key_of_minimum = population_keys[idx]
    return key_of_minimum


def generate_mating_pool_from_archive(archive: dict, mating_pool_size: int) -> dict:
    if mating_pool_size >= len(archive):
        return archive
    mating_pool = {}
    while len(archive) < mating_pool_size:
        key_to_add = tournament_selection_without_replacement(archive, 2)
        mating_pool[key_to_add] = archive[key_to_add]
    return mating_pool


def recombine(sequence1: str, sequence2: str, number_of_sites: int) -> str:
    # random positions with no duplicates
    crossover_sites = {random.randint(0, len(sequence1)) for _ in range(number_of_sites)}
    crossover_sites = list(crossover_sites)
    crossover_sites.sort()
    # assuming sequences are randomly chosen, it is ok to pick 1 to go first
    current_sequence = sequence1
    last_pos = 0
    new_seq = ''
    for site in crossover_sites:
        new_seq.join(current_sequence[last_pos:site])
        current_sequence = sequence1 if current_sequence == sequence2 else sequence2
        last_pos = site
    return new_seq


def generate_population_from_archive(params: dict) -> dict:
    """

    :param params: dict that contains:
        archive: dict which contains the sequences
        population size: int which tells the size of final population
        num crossover sites: int
        mutation %: int for percent mutation chance
        codon usage: dict to know the mapping of different codons
    :return:
    """
    # don't touch archive
    archive_copy = copy.deepcopy(params['archive'])
    # create mating pool from subset of archive
    mating_pool = generate_mating_pool_from_archive(archive_copy, params['mating pool size'])
    # until population is full
    #   choose 2 randomly selected individuals in mating pool
    #       apply recombination to create a new child
    #       run child under mutations
    new_population = {}
    mating_pool_keys = list(mating_pool.keys())

    while len(new_population) < params['population_size']:
        idx1 = random.randint(0, len(mating_pool_keys))
        idx2 = random.randint(0, len(mating_pool_keys))
        if idx1 != idx2:
            parent_a = mating_pool[mating_pool_keys[idx1]]
            parent_b = mating_pool[mating_pool_keys[idx2]]
            child_seq = recombine(parent_a[__SEQUENCE_KEY__], parent_b[__SEQUENCE_KEY__], params['num crossover sites'])
            child = {__SEQUENCE_KEY__: child_seq}
        else:
            child = mating_pool[mating_pool_keys[idx1]]
        child = mutate_seq(child, params)
        new_population[uuid.uuid4()] = child
    return new_population

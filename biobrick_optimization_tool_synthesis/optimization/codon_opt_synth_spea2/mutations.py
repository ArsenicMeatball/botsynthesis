import copy
import logging
import random
import uuid

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Data import CodonTable

import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions as fit_func
import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.dict_functions as dictf
__SEQUENCE_KEY__ = 'sequence'


def mutate_codon(codon: str, codon_table: CodonTable = CodonTable.generic_by_id[1]):
    """ Takes a codon and attempts to mutate it
    :param codon_table:
    :param codon:
    :return:
    """
    if len(codon) != 3:
        raise ValueError('Codon {0} is not of size 3'.format(codon))  # not handled by biopython
    # get our amino acid to figure out what codon possibilities there are
    amino = codon_table.forward_table[codon]
    full_back_table = dictf.invert_dict(codon_table.forward_table)
    possible_codons = full_back_table[amino]
    if len(possible_codons) == 1:
        logging.info('wasted a call on mutate codon - codon is unique and cannot be changed')
        return codon
    # remove original codon
    possible_codons.remove(codon)
    # pick another
    new_codon = random.choice(list(possible_codons))
    return new_codon


def mutate_seq(sequence: str, param: dict) -> str:
    """
        mutate a DNA Seq based on the mutation probability, returns a different DNA Seq
    :param sequence:
    :param param:
    :return:
    """
    if sequence == '':
        raise ValueError('Sequence is empty')
    new_sequence = ''
    logging.debug('Mutation % {0}, Codon Usage {1}, Unmutated Sequence {2}'.format(
        param['mutation %'], param['codon usage'], sequence))
    codons = []
    for idx in range(0, len(sequence), 3):
        # either add the og codon or mutated boy
        codon = str(sequence[idx:idx + 3])
        if random.randint(1, 100) > param['mutation %']:
            codon = mutate_codon(codon, param['codon usage'])
        codons.append(codon)
    new_sequence.join(codons)
    logging.debug('Mutated sequence {0}'.format(new_sequence))
    return new_sequence


def initialize_population(algorithm_parameters: dict) -> dict:
    """
        create population of mutants based on the codon optimized sequence
    :param algorithm_parameters: dict containing the parameters for the SPEA2 algorithm
    :return: dict representing the population
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
            logging.debug('failed to create a "new" sequence')
    return population


def tournament_selection_without_replacement(population: dict, n_ary: int = 2, minima: bool = True,
                                             fitness_key_name: str = fit_func.__FITNESS_KEY__) -> str:
    """
    Get winner of a tournament (key of winner) based on the fitness value
    :param fitness_key_name:
    :param minima: look for lowest value if True(default), False looks for largest value
    :param n_ary: number of individuals to conduct tournament on, default = binary
    :param population: all the individuals we are selecting from
    :return: the key of the winner of tournament
    """
    if n_ary < 2:
        raise ValueError('tournament selection must be binary or larger, currently {0}'.format(n_ary))
    population_keys = list(population.keys())
    # choose a random one to start
    current_key = random.choice(population_keys)
    for _ in range(1, n_ary):
        contender_key = random.choice(population_keys)
        val_current = population[current_key][fitness_key_name]
        val_contender = population[contender_key][fitness_key_name]
        replace = val_contender < val_current if minima else val_contender > val_current
        if replace:
            current_key = contender_key
    # key that won tournament
    return current_key


def generate_mating_pool_from_archive(archive: dict, mating_pool_size: int) -> dict:
    if len(archive) == 0:
        raise ValueError('archive cannot be empty')
    if mating_pool_size >= len(archive):
        return archive
    mating_pool = {}
    while len(mating_pool) < mating_pool_size:
        key_to_add = tournament_selection_without_replacement(archive, 2)
        mating_pool[key_to_add] = archive[key_to_add]
    return mating_pool


def recombine_dna_sequence(sequence1: str, sequence2: str, number_of_sites: int) -> str:
    # recombination we must ensure at codon lengths only
    if number_of_sites < 1:
        raise ValueError('number of sites needs to be greater than 0')
    if len(sequence1) != len(sequence2):
        raise ValueError('Cannot recombine sequences of different sizes')
    # randomly select all the crossover sites
    all_possible_crossover_sites = [num for num in range(0, len(sequence1), 3)]
    # pare down the crossover sites (may result in less than number of desired sites but that is ok. RNGesus has spoken
    crossover_sites = {random.choice(all_possible_crossover_sites) for _ in range(number_of_sites)}
    crossover_sites = list(crossover_sites)
    crossover_sites.sort()
    # assuming sequences are randomly chosen, it is ok to pick 1 to go first
    current_sequence = sequence1
    last_pos = 0
    new_seq = ''
    chunks = []
    for site in crossover_sites:
        chunks.append(current_sequence[last_pos:site])
        current_sequence = sequence1 if current_sequence == sequence2 else sequence2
        last_pos = site
    # add the last bit
    chunks.append(current_sequence[last_pos:])
    # add them all together
    new_seq = new_seq.join(chunks)
    if len(new_seq) != len(sequence1):
        raise ValueError(
            'Recombined sequence {0} is not the same size as its ancestor {1}\nBuilt from chunks {2}\nParents:\n{3}\n{4}'.format(
                len(new_seq),
                len(sequence1),
                chunks,
                sequence1,
                sequence2,
            )
        )
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
    sequences_set = set()
    attempts = 0
    max_attempts = 2 * params['population size']
    while len(new_population) < params['population size'] and attempts < max_attempts:
        idx1 = random.randint(0, len(mating_pool_keys) - 1)
        idx2 = random.randint(0, len(mating_pool_keys) - 1)
        if idx1 != idx2:
            parent_a = mating_pool[mating_pool_keys[idx1]]
            parent_b = mating_pool[mating_pool_keys[idx2]]
            child_seq = recombine_dna_sequence(parent_a[__SEQUENCE_KEY__], parent_b[__SEQUENCE_KEY__],
                                               params['num crossover sites'])
        else:
            child_seq = mating_pool[mating_pool_keys[idx1]][__SEQUENCE_KEY__]
        child_seq = mutate_seq(child_seq, params)
        if child_seq not in sequences_set:
            sequences_set.add(child_seq)
            child = {__SEQUENCE_KEY__: child_seq}
            new_population[uuid.uuid4()] = child
        else:
            attempts += 1
    return new_population

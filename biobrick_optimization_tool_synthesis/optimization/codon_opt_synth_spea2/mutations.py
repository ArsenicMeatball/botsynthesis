import logging
import random
import uuid

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions as fit_func

__SEQUENCE_KEY__ = 'sequence'


def mutate_codon(sequence: Seq, codon_usage: dict, idx=0, translation_dict: dict = None):
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


def mutate_seq(sequence: Seq, param: dict) -> Seq:
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
            if random.randint(1, 100) > param['mutation_%']
            else str(mutate_codon(sequence, param['codon_usage'], idx))
        )
        # print('Current sequence: ', new_sequence)

    return Seq(new_sequence, IUPAC.unambiguous_dna)


def initialize_population(algorithm_parameters: dict) -> tuple:
    """
        create population of mutants based on the codon optimized sequence
    :param algorithm_parameters: dict containing the parameters for the SPEA2 algorithm
    :return: tuple(dict representing the population, set representing the sequences in the population)
    """
    population = {}
    sequences = set()
    attempts = 0
    while len(population) < algorithm_parameters['population_size'] and \
            attempts < algorithm_parameters['max init population attempts']:
        seq = str(mutate_seq(algorithm_parameters['codon_opt_seq'], algorithm_parameters))
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

import random
import copy
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from biobrick_optimization_tool_synthesis.optimization.logic import cleaning, test_parameters
from biobrick_optimization_tool_synthesis.optimization.logic.codon_optimization_host import \
    determine_ideal_codon_optimized_sequence
from biobrick_optimization_tool_synthesis.optimization.logic.set_codon_table import fetch_codon_table
from biobrick_optimization_tool_synthesis.optimization.logic.test_parameters import algorithm_params


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


def initialize_population(algorithm_parameters: dict) -> dict:
    """
        create population of mutants based on the codon optimized sequence
    :param algorithm_parameters: dict containing the parameters for the SPEA2 algorithm
    :return: dict representing the population
    """
    population = {}
    while len(population) < algorithm_parameters['population_size']:
        population[str(mutate_seq(algorithm_parameters['codon_opt_seq'], algorithm_parameters))] = {}
        # print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    return population


if __name__ == '__main__':
    pop = initialize_population(algorithm_params)
    for k, v in pop.items():
        print(k, v)

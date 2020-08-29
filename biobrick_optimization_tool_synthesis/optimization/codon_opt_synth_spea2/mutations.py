import random
import uuid

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.test_parameters import algorithm_params

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
    while len(population) < algorithm_parameters['population_size']:
        seq = str(mutate_seq(algorithm_parameters['codon_opt_seq'], algorithm_parameters))
        if seq not in sequences:
            sequences.add(seq)
            seq_id = uuid.uuid4()
            population[seq_id] = {__SEQUENCE_KEY__: seq}
    return population, sequences


if __name__ == '__main__':
    pop = initialize_population(algorithm_params)
    for k, v in pop.items():
        print(k, v)

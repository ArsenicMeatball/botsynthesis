from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from biobrick_optimization_tool_synthesis.optimization.logic import interpret_sequence, test_parameters
from biobrick_optimization_tool_synthesis.optimization.logic.codon_optimization_host import \
    determine_ideal_codon_optimized_sequence
from biobrick_optimization_tool_synthesis.optimization.logic.set_codon_table import fetch_codon_table


def initialize_population(population_size, problem_size, ancestor_sequence, probability_mutation, codon_use_table):
    # create a population based on the expression optimized sequence

    population = []
    mutable_seq = ancestor_sequence.tomutable()
    if problem_size < population_size:
        size = population_size
    else:
        size = problem_size
    for idx in range(size):
        population.append(SequenceContainer(mutable_seq.toseq()))
        mutate_sequence(population[-1], codon_use_table, probability_mutation)

    return population


if __name__ == '__main__':
    sequence = Seq(interpret_sequence.clean_sequence(test_parameters.valid_protein), IUPAC.protein)
    codon_usage = fetch_codon_table()
    codon_optimized = determine_ideal_codon_optimized_sequence(sequence, codon_usage)
    algorithm_parameters = {
        'population_size': 20,
        'num_param_considered': 8,
        'ancestor_sequence': codon_optimized,
        'codon_usage': codon_usage,
        'mutation_%': 5,
        'archive_size': 10,
        'crossover_%': 15,
        'gc_parameters': {
            'min': 0.1,
            'max': 0.9,
            'window_size': 20
        },
        'restriction_sites': 'EcoRI'
    }

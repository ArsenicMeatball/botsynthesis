from typing import Union

from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions import score_names

neither_equal_case = 'neither'
dominates_case = 'dominates'
dominated_case = 'dominated'
strength_key = 'strength'
raw_fitness_key = 'raw fitness'


def compare_two_solutions_for_dominance(point_a: dict, point_b: dict, coordinates: list = score_names) -> tuple:
    """ Determines if one of the points is dominated:
    Domination (for minima):
        x dominates y, if for all coordinates, x's coordinate is smaller or equal to y's coordinate
    eg:
        (2, 10) vs (5, 6)
            2 < 5, left wins
            10 > 6, right wins
        therefore neither dominates
        and would return neither, neither
        (1, 2) vs (100, 1000)
            1 < 100, left wins
            2 < 1000, left wins
        for ALL coordinates, left is smaller, therefore left dominates right as a minima
        and would return dominated, dominates
    :param coordinates: What parameters we are checking for
    :param point_a: the first point to compare
    :param point_b: the second point to compare
    :return: tuple of bools, True is dominated, False is non-dominated
    """
    tally_a = 0
    tally_b = 0
    for coordinate in coordinates:
        if point_a[coordinate][0] < point_b[coordinate][0]:
            tally_a += 1
        elif point_a[coordinate][0] > point_b[coordinate][0]:
            tally_b += 1
    # neither dominates
    if tally_a > 0 and tally_b > 0:
        return neither_equal_case, neither_equal_case
    # always equal case
    elif tally_a == tally_b == 0:
        return neither_equal_case, neither_equal_case
    # b dominates
    elif tally_a == 0:
        return dominated_case, dominates_case
    # a dominates
    elif tally_b == 0:
        return dominates_case, dominated_case
    raise RuntimeError('should not be able to get here')


def parse_dominated_solutions_for_strength(left_x_right: dict, population: dict) -> dict:
    """
    Parse dict for
    strength value : how many solutions it dominates
    :param left_x_right:
    :param population:
    :return: population but with proper domination values appended to dict
    """
    for seq1, seq2_dict in left_x_right.items():
        strength = 0
        for seq2 in seq2_dict.keys():
            if left_x_right[seq1][seq2] == dominates_case:
                strength += 1
        population[seq1][strength_key] = strength
    return population


def parse_dominated_solutions_for_raw_fitness(left_x_right: dict, population: dict) -> dict:
    """
    Parse dict for:
    raw fitness value: the sum of the strengths of solutions which dominate it
    eg:
        a dominates b and c (S=2)
        b dominates c (S=1)
        c never dominates (S=0)
        therefore the raw fitness (rf) will be:
        rf_a = 0
        rf_b = S_a = 2
        rf_c = S_a + S_b = 2 + 1 = 3

    :param left_x_right:
    :param population:
    :return:
    """
    for seq1, seq2_dict in left_x_right.items():
        raw_fitness = 0
        for seq2 in seq2_dict.keys():
            if left_x_right[seq1][seq2] == dominated_case:
                raw_fitness += population[seq2][strength_key]
        population[seq1][raw_fitness_key] = raw_fitness
    return population


def find_dominated_solutions(population: dict) -> dict:
    """

    :param population: (dict)
    :return: (dict)
    """
    # test all pairs to see if they get dominated O(dn^2)
    left_x_right = {}
    visited = set()
    for seq1 in population.keys():
        for seq2 in population.keys():
            # prevent key errors
            if seq1 not in left_x_right:
                left_x_right[seq1] = {seq2: neither_equal_case}
            if seq2 not in left_x_right:
                left_x_right[seq2] = {seq1: neither_equal_case}

            if (seq1, seq2) not in visited:
                left = population[seq1]
                right = population[seq2]
                left_x_right[seq1][seq2], left_x_right[seq2][seq1] = compare_two_solutions_for_dominance(
                    left, right
                )
                visited.add((seq1, seq2))
                visited.add((seq2, seq1))
    # update population
    population = parse_dominated_solutions_for_strength(left_x_right, population)
    population = parse_dominated_solutions_for_raw_fitness(left_x_right, population)

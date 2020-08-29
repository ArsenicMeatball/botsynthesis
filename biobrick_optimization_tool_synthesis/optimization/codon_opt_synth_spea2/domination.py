import biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions as fit_funcs

__NEITHER_EQUAL_CASE__ = 'neither'
__DOMINATES_CASE__ = 'dominates'
__DOMINATED_CASE__ = 'dominated'
__STRENGTH_KEY__ = 'strength'
__RAW_FITNESS_KEY__ = 'raw fitness'


def calculate_raw_fitness(population):
    """
    Calculates strength and raw fitness values for every member of the population
    :param population: dict of sequences and values
    :return: None, updates population
    """
    left_x_right = find_dominated_solutions(population)
    calculate_strength(left_x_right, population)
    parse_for_raw_fitness(left_x_right, population)


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
                left_x_right[seq1] = {seq2: __NEITHER_EQUAL_CASE__}
            if seq2 not in left_x_right:
                left_x_right[seq2] = {seq1: __NEITHER_EQUAL_CASE__}

            if (seq1, seq2) not in visited:
                left = population[seq1]
                right = population[seq2]
                left_x_right[seq1][seq2], left_x_right[seq2][seq1] = compare_two_solutions_for_dominance(
                    left, right, fit_funcs.__SCORE_NAMES__
                )
                visited.add((seq1, seq2))
                visited.add((seq2, seq1))
    return left_x_right


def compare_two_solutions_for_dominance(
        point_a: dict,
        point_b: dict,
        coordinates: list,
) -> tuple:
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
        return __NEITHER_EQUAL_CASE__, __NEITHER_EQUAL_CASE__
    # always equal case
    elif tally_a == tally_b == 0:
        return __NEITHER_EQUAL_CASE__, __NEITHER_EQUAL_CASE__
    # b dominates
    elif tally_a == 0:
        return __DOMINATED_CASE__, __DOMINATES_CASE__
    # a dominates
    elif tally_b == 0:
        return __DOMINATES_CASE__, __DOMINATED_CASE__
    raise RuntimeError('should not be able to get here')


def calculate_strength(left_x_right: dict, population: dict) -> dict:
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
            if left_x_right[seq1][seq2] == __DOMINATES_CASE__:
                strength += 1
        population[seq1][__STRENGTH_KEY__] = strength
    return population


def parse_for_raw_fitness(left_x_right: dict, population: dict):
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
    :return: None, updates population
    """
    for seq1, seq2_dict in left_x_right.items():
        raw_fitness = 0
        for seq2 in seq2_dict.keys():
            if left_x_right[seq1][seq2] == __DOMINATED_CASE__:
                raw_fitness += population[seq2][__STRENGTH_KEY__]
        population[seq1][__RAW_FITNESS_KEY__] = raw_fitness

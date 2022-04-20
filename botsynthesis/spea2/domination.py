from botsynthesis.utils.constants import NEITHER_EQUAL_CASE, SCORE_NAMES, \
    DOMINATED_CASE, DOMINATES_CASE, STRENGTH_KEY, RAW_FITNESS_KEY


def calculate_raw_fitness(population):
    """
    Calculates strength and raw fitness values for every member of the
    population
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
    for seq_id1 in population.keys():
        for seq_id2 in population.keys():
            # prevent key errors
            if seq_id1 not in left_x_right:
                left_x_right[seq_id1] = {seq_id2: NEITHER_EQUAL_CASE}
            if seq_id2 not in left_x_right:
                left_x_right[seq_id2] = {seq_id1: NEITHER_EQUAL_CASE}

            if (seq_id1, seq_id2) not in visited:
                left = population[seq_id1]
                right = population[seq_id2]
                (
                    left_x_right[seq_id1][seq_id2],
                    left_x_right[seq_id2][seq_id1],
                ) = compare_two_solutions_for_dominance(
                    left, right, SCORE_NAMES
                )
                visited.add((seq_id1, seq_id2))
                visited.add((seq_id2, seq_id1))
    return left_x_right


def compare_two_solutions_for_dominance(
        point_a: dict,
        point_b: dict,
        coordinates: list,
) -> tuple:
    """Determines if one of the points is dominated:
    Domination (for minima):
        x dominates y, if for all coordinates, x's coordinate is smaller or
        equal to y's coordinate
    eg:
        (2, 10) vs (5, 6)
            2 < 5, left wins
            10 > 6, right wins
        therefore neither dominates
        and would return neither, neither
        (1, 2) vs (100, 1000)
            1 < 100, left wins
            2 < 1000, left wins
        for ALL coordinates, left is smaller, therefore left dominates right
        as a minima
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
        return NEITHER_EQUAL_CASE, NEITHER_EQUAL_CASE
    # always equal case
    elif tally_a == tally_b == 0:
        return NEITHER_EQUAL_CASE, NEITHER_EQUAL_CASE
    # b dominates
    elif tally_a == 0:
        return DOMINATED_CASE, DOMINATES_CASE
    # a dominates
    elif tally_b == 0:
        return DOMINATES_CASE, DOMINATED_CASE
    raise RuntimeError("should not be able to get here")


def calculate_strength(left_x_right: dict, population: dict) -> dict:
    """
    Parse dict for
    strength value : how many solutions it dominates
    :param left_x_right:
    :param population:
    :return: population but with proper domination values appended to dict
    """
    for seq_id1, seq2_dict in left_x_right.items():
        strength = 0
        for seq_id2 in seq2_dict.keys():
            if left_x_right[seq_id1][seq_id2] == DOMINATES_CASE:
                strength += 1
        population[seq_id1][STRENGTH_KEY] = strength
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
    for seq_id1, seq2_dict in left_x_right.items():
        raw_fitness = 0
        for seq_id2 in seq2_dict.keys():
            if left_x_right[seq_id1][seq_id2] == DOMINATED_CASE:
                raw_fitness += population[seq_id2][STRENGTH_KEY]
        population[seq_id1][RAW_FITNESS_KEY] = raw_fitness

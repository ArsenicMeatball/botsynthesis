from typing import Union
from biobrick_optimization_tool_synthesis.optimization.codon_opt_synth_spea2.fitness_functions import score_names


def compare_two_solutions_for_dominance(point_a: dict, point_b: dict, coordinates: list = score_names) -> tuple:
    """ Determines if one of the points is dominated:
    Domination (for minima):
        x dominates y, if for all coordinates, x's coordinate is smaller or equal to y's coordinate
    eg:
        (2, 10) vs (5, 6)
            2 < 5, left wins
            10 > 6, right wins
        therefore neither dominates
        and would return False, False
        (1, 2) vs (100, 1000)
            1 < 100, left wins
            2 < 1000, left wins
        for ALL coordinates, left is smaller, therefore left dominates right as a minima
        and would return False, True
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
    if tally_a > 0 and tally_b > 0:
        return False, False
    if tally_a == tally_b == 0:
        return False, False
    if tally_a == 0:
        return False, True
    if tally_b == 0:
        return True, False
    raise RuntimeError('should not be able to get here')


def find_dominated_solutions(points: dict) -> tuple:
    # test all pairs to see if they get dominated O(dn^2)
    dominated = []
    for point in points:
        point['dominated'] = False
    for point_a in points:
        if not point_a['dominated']:
            for point_b in points:
                if not point_b['dominated']:
                    point_a['dominated'], point_b['dominated'] = compare_two_solutions_for_dominance(point_a, point_b)
                    if point_a['dominated']:
                        dominated.append(point_a)
                    elif point_b['dominated']:
                        dominated.append(point_b)
    non_dominated = [point not in dominated for point in points]
    return dominated, non_dominated

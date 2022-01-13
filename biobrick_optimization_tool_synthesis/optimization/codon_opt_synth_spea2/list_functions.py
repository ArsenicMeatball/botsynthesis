def differences_between_two_lists(list1: list, list2: list) -> list:
    """gets differences between two lists
    :param list1: (list) of elements
    :param list2: (list) of elements
    :return (set) differences
    """
    return list(set(list1).symmetric_difference(set(list2)))


def number_of_differences_between_two_lists(list1: list, list2: list) -> int:
    """gets number of differences between two lists
    :param list1: (list) of elements
    :param list2: (list) of elements
    :return (int) num differences
    """
    return len(differences_between_two_lists(list1, list2))

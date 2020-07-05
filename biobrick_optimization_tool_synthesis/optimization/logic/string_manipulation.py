import re


def find_number_of_non_overlapping_repeats(string: str, min_repeat_size=10):
    if min_repeat_size < 1:
        raise ArithmeticError("Can't find repeats smaller than 1")
    if len(string) < min_repeat_size:
        return 0
    tracker = set()
    result = 0
    for idx in range(len(string) - min_repeat_size):
        substring = string[idx:idx + min_repeat_size]
        if substring not in tracker:
            tracker.add(substring)
            result += string.count(substring) - 1
    return result


def find_number_of_overlapping_repeats(string: str, min_repeat_size=10):
    if min_repeat_size < 1:
        raise ArithmeticError("Can't find repeats smaller than 1")
    if len(string) < min_repeat_size:
        return 0
    maximum = len(string) - min_repeat_size + 1
    strings_of_size = {string[idx:idx + min_repeat_size] for idx in range(maximum)}
    return maximum - len(strings_of_size)


def find_repeats(string: str, min_repeat_size=10, overlapping=True):
    if min_repeat_size < 1:
        raise ArithmeticError("Can't find repeats smaller than 1")
    if len(string) < min_repeat_size:
        return {}
    result = {}
    visited = set()
    for idx in range(len(string) - min_repeat_size):
        substring = string[idx:idx + min_repeat_size]
        if substring not in visited:
            visited.add(substring)
            pattern = substring if not overlapping else '(?=%s)' % substring
            matches = [m.start() for m in re.finditer(pattern, string)]
            if len(matches) > 1:
                result[substring] = matches
    return result

import re


def get_number_of_repeats_from_dict(d: dict) -> int:
    r = 0
    for v in d.values():
        r += len(v)
    return r - len(d)


def find_number_of_non_overlapping_repeats(string: str, min_repeat_size=10) -> int:
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


def find_number_of_overlapping_repeats(string: str, min_repeat_size=10) -> int:
    if min_repeat_size < 1:
        raise ArithmeticError("Can't find repeats smaller than 1")
    if len(string) < min_repeat_size:
        return 0
    maximum = len(string) - min_repeat_size + 1
    strings_of_size = {string[idx:idx + min_repeat_size] for idx in range(maximum)}
    return maximum - len(strings_of_size)


def find_repeats(string: str, min_repeat_size=10, overlapping=True) -> dict:
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

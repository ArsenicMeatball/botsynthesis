import re
import logging


def get_number_of_repeats_from_repeats_dict(d: dict) -> int:
    r = 0
    for v in d.values():
        r += len(v)
    return r - len(d)


def find_number_of_non_overlapping_repeats(
    string: str, min_repeat_size=10
) -> int:
    if min_repeat_size < 1:
        raise ArithmeticError("Can't find repeats smaller than 1")
    if len(string) < min_repeat_size:
        return 0
    tracker = set()
    result = 0
    for idx in range(len(string) - min_repeat_size):
        substring = string[idx : idx + min_repeat_size]
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
    strings_of_size = {
        string[idx : idx + min_repeat_size] for idx in range(maximum)
    }
    return maximum - len(strings_of_size)


def find_repeats(string: str, min_repeat_size=10, overlapping=True) -> dict:
    if min_repeat_size < 1:
        raise ArithmeticError("Can't find repeats smaller than 1")
    if len(string) < min_repeat_size:
        return {}
    result = {}
    visited = set()
    for idx in range(len(string) - min_repeat_size):
        substring = string[idx : idx + min_repeat_size]
        if substring not in visited:
            visited.add(substring)
            pattern = substring if not overlapping else "(?=%s)" % substring
            matches = [m.start() for m in re.finditer(pattern, string)]
            if len(matches) > 1:
                result[substring] = matches
    return result


def find_separated_palindromes(
    string: str,
    min_separation_size: int = 3,
    max_separation_size: int = 30,
    palindrome_size: int = 10,
) -> dict:
    result = {}
    for separation in range(min_separation_size, max_separation_size + 1):
        result[separation] = list()
        if palindrome_size * 2 + separation <= len(
            string
        ):  # check if sequence large enough
            for idx in range(
                palindrome_size - 1, len(string) - palindrome_size - separation
            ):  # need enough room
                counter = 0
                while counter < palindrome_size:
                    if (
                        string[idx - counter]
                        == string[idx + separation + 1 + counter]
                    ):
                        counter += 1
                    else:
                        break
                if counter == palindrome_size:
                    # give the start, separation, and end
                    result[separation] += [
                        idx - palindrome_size + 1,
                        string[idx + 1 : idx + separation + 1],
                        idx + separation + palindrome_size,
                    ]
                    logging.debug(
                        "Palindrome (L | Sep | R): {} | {} | {}".format(
                            string[idx - palindrome_size + 1 : idx + 1],
                            string[idx + 1 : idx + separation + 1],
                            string[
                                idx
                                + separation
                                + 1 : idx
                                + separation
                                + palindrome_size
                                + 1
                            ],
                        )
                    )
        if result[separation] == list():
            # remove if empty
            result.pop(separation)
    return result


def find_num_differences(str1: str, str2: str) -> int:
    return sum(1 for a, b in zip(str1, str2) if a != b)

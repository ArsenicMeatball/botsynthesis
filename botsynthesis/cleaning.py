from Bio import Restriction


def clean_sequence(sequence: str):
    return ("".join(sequence.split())).upper()


def get_rest_enzymes_from_list(restriction_enzymes: list):
    return Restriction.RestrictionBatch(
        [Restriction.AllEnzymes.get(enz) for enz in restriction_enzymes]
    )


def get_rest_enzymes_from_string(restriction_enzymes: str):
    return get_rest_enzymes_from_list(restriction_enzymes.split())

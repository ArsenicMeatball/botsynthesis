# uses python codon tables which uses the
# http://www.kazusa.or.jp/codon/readme_codon.html
# Codon usage tabulated from the international DNA sequence databases:
# status for the year 2000.
# Nakamura, Y., Gojobori, T. and Ikemura, T. (2000) Nucl. Acids Res. 28, 292.

import python_codon_tables as pct


def fetch_codon_table(organism_id=316407) -> dict:
    """
        grab a dict representing the codon_table w/ codon usage of the organism
        ('amino acid', {'codon' : ratio, ... })\
        ('A', {'GCA': 0.21, 'GCC': 0.27, 'GCG': 0.36, 'GCT': 0.16})
    :param organism_id: the id of the organism in kasuza, though e_coli

    will work
    :return: (dict)
    """
    return pct.get_codons_table(organism_id)


if __name__ == "__main__":
    print(fetch_codon_table(316407))
    print(fetch_codon_table("316407"))
    print(fetch_codon_table("e_coli"))
    print("done")

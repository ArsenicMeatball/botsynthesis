from Bio.Data import CodonTable

table = CodonTable.generic_by_id[1]
print(table.forward_table)
print(table.back_table)
back_table_non_naive = {}
for codon, amino in table.forward_table.items():
    if amino in back_table_non_naive:
        back_table_non_naive[amino].append(codon)
    else:
        back_table_non_naive[amino] = [codon]
print(back_table_non_naive)

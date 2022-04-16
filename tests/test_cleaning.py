# import unittest
#
# from Bio import Restriction
# from Bio.Alphabet import IUPAC
# from Bio.Seq import Seq
#
# import botsynthesis.cleaning as cleaning
#
#
# class Cleaning(unittest.TestCase):
#     def test_determine_sequence_type_pass(self):
#         expected = cleaning.INDETERMINATE, [True, True, True]
#         sequence = "A"
#         actual = cleaning.determine_sequence_type(sequence)
#         self.assertEqual(
#             expected,
#             actual,
#             "failed to determine correct sequence type - indeterminate",
#         )
#
#         expected = cleaning.INVALID
#         sequence = "$!^&E_)*(^8678969386"
#         actual = cleaning.determine_sequence_type(sequence)
#         self.assertEqual(
#             expected,
#             actual,
#             "failed to determine correct sequence type - invalid",
#         )
#
#         expected = cleaning.RNA
#         sequence = "ACUG"
#         actual = cleaning.determine_sequence_type(sequence)
#         self.assertEqual(
#           expected, actual, "failed to determine correct sequence type - RNA"
#         )
#
#         expected = cleaning.AA
#         sequence = "PRTK"
#         actual = cleaning.determine_sequence_type(sequence)
#         self.assertEqual(
#           expected, actual, "failed to determine correct sequence type - AA"
#         )
#
#     def test_determine_sequence_type_fail(self):
#         expected = cleaning.INDETERMINATE, [True, True, False]
#         sequence = "ACTG"
#         actual = cleaning.determine_sequence_type(sequence)
#         self.assertEqual(
#             expected,
#             actual,
#             "failed to determine sequence that could be AA or DNA",
#         )
#
#         sequence = {"aaaaa": "aaaaaa", "bbbbbbb": "bbbbb"}
#         self.assertRaises(
#             TypeError, cleaning.determine_sequence_type(sequence)
#         )
#
#     def test_clean_sequence_pass(self):
#         expected = "ACTGGCGTACGTCATGCATCTACACGATCG"
#         sequence = "act ggc gt acg tca tg ca tc tac ac g at cg"
#         actual = cleaning.clean_sequence(sequence)
#       self.assertEqual(expected, actual, "failed to clean sequence properly")
#
#     def test_clean_sequence_fail(self):
#         sequence = {"aaaaa": "aaaaaa", "bbbbbbb": "bbbbb"}
#         self.assertRaises(AttributeError, cleaning.clean_sequence, sequence)
#
#     def test_turn_string_sequence_into_amino_pass(self):
#         expected = Seq("TGVRHASTRS", IUPAC.protein)
#         sequence = "ACTGGCGTACGTCATGCATCTACACGATCG"
#         sequence_type = cleaning.DNA
#         actual = cleaning.turn_string_sequence_into_amino(
#             sequence, sequence_type
#         )
#         self.assertEqual(
#             expected, actual, "failed to convert DNA string into AMINO Seq"
#         )
#
#         expected = Seq("TGVRHASTRS", IUPAC.protein)
#         sequence = "ACTGGCGTACGTCATGCATCTACACGATCG".replace("T", "U")
#         sequence_type = cleaning.RNA
#         actual = cleaning.turn_string_sequence_into_amino(
#             sequence, sequence_type
#         )
#         self.assertEqual(
#             expected, actual, "failed to convert RNA string into AMINO SEQ"
#         )
#
#         expected = Seq("TGVRHASTRS", IUPAC.protein)
#         sequence = "TGVRHASTRS"
#         sequence_type = cleaning.AA
#         actual = cleaning.turn_string_sequence_into_amino(
#             sequence, sequence_type
#         )
#         self.assertEqual(
#             expected, actual, "failed to convert AA string into AMINO SEQ"
#         )
#
#         expected = Seq("TGVRHASTRS", IUPAC.protein)
#         sequence = "ACTGGCGTACGTCATGCATCTACACGATCG".replace("T", "U")
#         actual = cleaning.turn_string_sequence_into_amino(sequence)
#         self.assertEqual(
#             expected,
#             actual,
#             "failed to convert an undocumented RNA string into AMINO SEQ",
#         )
#
#         expected = Seq("TGVRHASTRS", IUPAC.protein)
#         sequence = "TGVRHASTRS"
#         actual = cleaning.turn_string_sequence_into_amino(sequence)
#         self.assertEqual(
#             expected,
#             actual,
#             "failed to convert an undocumented AA string into AMINO SEQ",
#         )
#
#     def test_turn_string_sequence_into_amino_fail(self):
#         sequence = "ACTGGCGTACGTCATGCATCTACACGATCG"
#         self.assertRaises(
#             TypeError, cleaning.turn_string_sequence_into_amino, sequence
#         )
#         sequence = "$^!%E()&*%$E"
#         self.assertRaises(
#             TypeError, cleaning.turn_string_sequence_into_amino, sequence
#         )
#
#     def test_get_rest_enzymes_from_list_pass(self):
#       input_list = ["NdeI", "XhoI", "HpaI", "PstI", "EcoRV", "NcoI", "BamHI"]
#         expected = Restriction.RestrictionBatch(
#             [Restriction.AllEnzymes.get(enz) for enz in input_list]
#         )
#         actual = cleaning.get_rest_enzymes_from_list(input_list)
#         self.assertEqual(
#           expected, actual, "failed to pull the correct restriction enzymes"
#         )
#
#     def test_get_rest_enzymes_from_list_fail(self):
#         input_list = [
#             "REEEEEEEEEEEEEEEEE",
#             "XhoI",
#             "HpaI",
#             "PstI",
#             "EcoRV",
#             "NcoI",
#             "BamHI",
#         ]
#         self.assertRaises(
#             ValueError, cleaning.get_rest_enzymes_from_list, input_list
#         )
#
#     def test_get_rest_enzymes_from_string_pass(self):
#         input_list = "XhoI HpaI PstI EcoRV NcoI BamHI"
#         expected = Restriction.RestrictionBatch(
#             [Restriction.AllEnzymes.get(enz) for enz in input_list.split()]
#         )
#         actual = cleaning.get_rest_enzymes_from_string(input_list)
#         self.assertEqual(
#           expected, actual, "failed to pull the correct restriction enzymes"
#         )
#
#     def test_get_rest_enzymes_from_string_fail(self):
#         input_list = "REEEEEEEEEEEEEEEEE XhoI HpaI PstI EcoRV NcoI BamHI"
#         self.assertRaises(
#             ValueError, cleaning.get_rest_enzymes_from_string, input_list
#         # )
#
#
# if __name__ == "__main__":
#     unittest.main()

import unittest
import dendropy
import os

from pangolin.scripts.lineage_assigner.lineage_finder import LineageFinder, all_equal, trim_to_common_ancestor, \
    get_basal_lineage
from pangolin.scripts.lineage_assigner.utils import collapse_nodes

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, "test", 'data', 'lineage_finder')


class LineageTests(unittest.TestCase):
    def test_all_equal(self):
        self.assertTrue(all_equal(["B.1.1", "B.1.1", "B.1.1"]))
        self.assertFalse(all_equal(["B.2.1", "B.1.1"]))

    def test_trim(self):
        self.assertEqual(trim_to_common_ancestor(["B.1.2", "B.1.3", "B.1.4"]), "B.1")
        self.assertEqual(trim_to_common_ancestor(["B.2.6", "B.1.5", "B.3.4"]), "B")

    def test_basal_lineage(self):
        simple = ["B.1.2", "B.1.3", "B.1.4"]
        a_little_tricky = ["A", "B.2"]
        extra_ticky = ["A.1", "A.1.2", "B.3", "C.3", "C.3"]

        self.assertEqual(get_basal_lineage(simple), "B.1")
        self.assertEqual(get_basal_lineage(a_little_tricky), "A")
        self.assertEqual(get_basal_lineage(extra_ticky), "A.1")

    def test_sibling(self):
        tree = dendropy.Tree.get(path="%s/representative_sequences.aln.fasta.nexus.tree" % data_dir, schema="nexus",
                                 preserve_underscores=True)
        finder = LineageFinder(tree, "Iceland/222/2020|EPI_ISL_417837|B.1.8|Iceland|||2020-03-16")

        self.assertEqual(finder.get_lineage(2, "|"), ["B.1.8", "77"])

    def test_nested_no_sibling(self):
        tree = dendropy.Tree.get(path="%s/representative_sequences.aln.fasta.nexus.tree" % data_dir, schema="nexus",
                                 preserve_underscores=True)
        finder = LineageFinder(tree, "France/IDF2256/2020|EPI_ISL_416498|B.1.4|France|||2020-03-11")

        self.assertEqual(finder.get_lineage(2, "|"), ["B.1.4", "99"])

    def test_between_clades_fall_back(self):
        tree = dendropy.Tree.get(path="%s/representative_sequences.aln.fasta.nexus.tree" % data_dir, schema="nexus",
                                 preserve_underscores=True)
        finder = LineageFinder(tree, "USA/CZB-RR057-013/2020|EPI_ISL_417937|B.1.3|USA|||2020-03-18")

        self.assertEqual(finder.get_lineage(2, "|"), ["B.1", "37"])


    def test_empty_subtype_in_name(self):
        tree = dendropy.Tree.get(path="%s/EPI_ISL_413581X.aln.fasta.nexus.tree" % data_dir, schema="nexus",
                                 preserve_underscores=True)

        collapse_nodes(tree, lambda x: x.edge.length == 0)
        finder = LineageFinder(tree, "Netherlands___Oss_1363500___2020|EPI_ISL_413581X|2020-02-29")

        self.assertEqual(finder.get_lineage(2, "|"), ["B", "96.00"])

    def test_polytomy_no_tip_Siblings(self):
        tree = dendropy.Tree.get(path="%s/EPI_ISL_416626.aln.fasta.nexus.tree" % data_dir, schema="nexus",
                                 preserve_underscores=True)
        collapse_nodes(tree, lambda x: x.edge.length == 0)
        finder = LineageFinder(tree, "EPI_ISL_416626||")
        self.assertEqual(finder.get_lineage(2, "|"), ["B","89.00"])



if __name__ == '__main__':
    unittest.main()

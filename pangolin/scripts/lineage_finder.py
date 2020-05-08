#!/usr/bin/env python3

import warnings


class LineageFinder:
    def __init__(self, tree, query_taxon, index, separator):
        self.tree = tree
        self.root = tree.seed_node
        self.index = index
        self.separator = separator
        query_node = self.tree.find_node_with_taxon_label(query_taxon)
        if query_node is None:
            raise KeyError("Taxon %s not found in tree" % query_taxon)

        self.query_node_parent = query_node.parent_node
        self.removed_query_node = self.query_node_parent.remove_child(query_node)
        self.at_root = False
        self.lineage = None
        self.label = None
        self.run()

    def run(self):
        self.get_lineage()
        self.get_label()

    def annotate_tips_from_label(self, index, separator):
        annotations = {}
        for tip in self.tree.leaf_node_iter():
            trait = {}
            value = tip.taxon.label.split(separator)[index]
            trait["lineage"] = value if len(value) > 0 else None
            annotations[tip.taxon.label] = trait

        self.annotate_tips(annotations)

    def annotate_tips(self, annotations):
        for tip in annotations:
            self.annotate_node(tip, annotations[tip])

    def annotate_node(self, tip_label, annotations):
        node = self.tree.find_node_with_taxon(lambda taxon: True if taxon.label == tip_label else False)
        if node is None:
            warnings.warn("Taxon: %s not found in tree" % tip_label)
        else:
            for a in annotations:
                if a != "taxon":
                    setattr(node, a, annotations[a])
                    node.annotations.add_bound_attribute(a)

    def lineage_parsimony(self, node):
        if node.is_leaf():
            return node.annotations.get_value("lineage")

        child_lineages = []
        for child in node.child_node_iter():
            child_lineages.append(self.lineage_parsimony(child))

        if all_equal(child_lineages):
            imputed_lineage = child_lineages[0]
        else:
            imputed_lineage = get_basal_lineage(child_lineages)

        if node.num_child_nodes() > 1 or (node == self.query_node_parent and node.child_nodes()[0].is_leaf()):
            setattr(node, "lineage", imputed_lineage)
            node.annotations.add_bound_attribute("lineage")

        return imputed_lineage

    def get_lineage(self):
        self.annotate_tips_from_label(self.index, self.separator)
        self.lineage_parsimony(self.root)
        lineage = self.query_node_parent.annotations.get_value("lineage")

        if lineage is None:
            if self.query_node_parent is not self.tree.seed_node:
                lineage = self.query_node_parent.parent_node.annotations.get_value("lineage")
            else:
                # the parent node is the root and and the ancestral lineage is A
                lineage = "A"
        if lineage == "A":
            self.at_root = True
        self.lineage = lineage

    def get_label(self):

        if self.at_root:
            # get all non A lineages
            nonALineage = []
            current_arlt = 100
            current_boot = 100

            for node in self.tree.preorder_node_iter(lambda n: n.annotations.get_value("lineage") is not None and (
                    n.annotations.get_value("lineage") == "B" or (
                    len(n.annotations.get_value("lineage").split(".")) == 2 and n.annotations.get_value("lineage")[
                0] == "A"))):
                lineage = node.annotations.get_value("lineage")
                if lineage not in nonALineage:
                    nonALineage.append(lineage)
                    alrt, boot = node.label.split("/")[-2:]
                    current_arlt -= (100 - float(alrt))
                    current_boot -= (100 - float(boot))

            self.alrt = int(current_arlt)
            self.boot = int(current_boot)
        else:
            lineage_mrca = \
                [node for node in
                 self.tree.preorder_node_iter(lambda n: n.annotations.get_value("lineage") == self.lineage)][0]
            label = lineage_mrca.label  # if lineage_mrca.label is not None else lineage_mrca.annotations.get_value("label")

            alrt, boot = label.split("/")[-2:]

            self.alrt = int(float(alrt))
            self.boot = int(float(boot))


def all_equal(lineage_list):
    for i in range(0, len(lineage_list) - 1):
        if lineage_list[i] != lineage_list[i + 1]:
            return False
    return True


def get_basal_lineage(lineage_list):
    #     Check we're on the same level A,B,C if not select the most basal
    base_clades = [base.split('.')[0] for base in lineage_list]
    if not all_equal(base_clades):
        current_lineages = [lineage for lineage in lineage_list if lineage[0] == sorted(base_clades)[0]]
    else:
        current_lineages = lineage_list

    most_basal_lineage_length = min([len(lineage) for lineage in current_lineages])
    standardized_depth_lineages = [lineage[0:most_basal_lineage_length] for lineage in current_lineages]
    return trim_to_common_ancestor(standardized_depth_lineages)


def trim_to_common_ancestor(lineage_list):
    if len(lineage_list) == 1 or all_equal(lineage_list):
        return lineage_list[0]
    else:
        return trim_to_common_ancestor([lineage[0:-2] for lineage in lineage_list])


def get_annotations(taxon_key, annotation_list):
    annotation_dict = {}
    for row in annotation_list:
        annotation_dict[row[taxon_key]] = row
    return annotation_dict

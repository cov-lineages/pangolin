#!/usr/bin/env python

import warnings

class LineageFinder:
    def __init__(self, tree, query_taxon):
        self.tree = tree
        self.root = tree.seed_node
        query_node = self.tree.find_node_with_taxon_label(query_taxon)
        if query_node is None:
            raise KeyError("Taxon %s not found in tree" % query_taxon)

        self.query_node_parent = query_node.parent_node
        self.removed_query_node = self.query_node_parent.remove_child(query_node)
        pass

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

    def get_lineage(self, index, separator):
        self.annotate_tips_from_label(index, separator)
        self.lineage_parsimony(self.root)

        label = self.query_node_parent.label if self.query_node_parent.label is not None else self.query_node_parent.annotations.get_value(
                    "label")

        if self.query_node_parent.annotations.get_value("lineage") is not None:
            return [self.query_node_parent.annotations.get_value("lineage"),
                    label]
        else:
            return [self.query_node_parent.parent_node.annotations.get_value("lineage"),
                    label]


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

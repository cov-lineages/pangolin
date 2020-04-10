# return parameter if not a string
import dendropy


def parse_tree(file, file_format):
    return dendropy.Tree.get(path=file, schema=file_format.lower(), preserve_underscores=True)


def write_tree(tree, file, format):
    tree.write(path=file, schema=format)


def collapse_nodes(tree, predicate):
    for node in tree.preorder_node_iter(predicate):
        if not node.is_leaf():
            node.edge.collapse(adjust_collapsed_head_children_edge_lengths=True)


def prepare_tree(input_file, file_format, collapse):
    tree = parse_tree(input_file, file_format)
    if collapse:
        collapse_nodes(tree, lambda node: node.edge.length == 0)
    return tree

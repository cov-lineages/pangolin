# return parameter if not a string
import dendropy


def parse_tree(file, format):
    return dendropy.Tree.get(path=file, schema=format.lower(), preserve_underscores=True)


def write_tree(tree, file, format):
    tree.write(path=file, schema=format)


def collapse_nodes(tree, threshold):
    for node in tree.postorder_node_iter(lambda n: not n.is_leaf()):
        if not node == tree.seed_node:
            if node.edge_length <= threshold:
                node.edge.collapse(adjust_collapsed_head_children_edge_lengths=True)


def prepare_tree(options, input_file_name=None):
    input_path = input_file_name if input_file_name is not None else options.input
    tree = parse_tree(input_path, options.format)
    if options.collapse:
        collapse_nodes(tree, options.collapse)
    return tree
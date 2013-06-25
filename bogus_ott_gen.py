import sys
import dendropy
from dendropy.treesplit import encode_splits,split_to_list
from dendropy.utility.error import DataParseError
def life_ott(fo):
    output = sys.stdout
    dataset = dendropy.DataSet()
    try:
        dataset.read(stream=fo, schema="Newick")
    except DataParseError as dfe:
        raise ValueError(str(dfe))
    tree_list = dataset.tree_lists[0]
    parent_id = '805080'
    branch_counter = 0
    tree_labels = set()
    for tree in tree_list:
        encode_splits(tree)
        tree_mask = tree.seed_node.edge.split_bitmask
        assert tree_mask is not None
        tree_tax = set(split_to_list(tree_mask))
        split_list = []
        for node in tree.leaf_iter():
            tree_labels.add(node.taxon.label)
    output.write('805080\t|\t\t|\tlife\t|\tno rank\t|\tncbi:1,gbif:0\t|\t\t|\t\t|\t\t|\t\n')
    for i in tree_labels:
        name, ottoid = i.split('@')
        output.write(ottoid+'\t|\t'+parent_id+'\t|\t'+name+'\t|\tspecies\t|\tncbi:1\t|\t\t|\t\n')

if __name__ == '__main__':
    try:
        filename = sys.argv[1]
    except:
        fo = sys.stdin
    else:
        fo = open(filename, 'rU')
    life_ott(fo)

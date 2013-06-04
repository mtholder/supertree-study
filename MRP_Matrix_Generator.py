#!/usr/bin/env python
import sys
from dendropy.treesplit import encode_splits,split_to_list
from cStringIO import StringIO
from dendropy import DataSet
from dendropy.utility.error import DataParseError
from dendropy.utility.textutils import escape_nexus_token

if __name__ == '__main__':
    output = sys.stdout
    fo = open(sys.argv[1], "rU")
    dataset = DataSet()
    try:
        dataset.read(stream=fo, schema="Newick")
    except DataParseError as dfe:
        raise ValueError(str(dfe))
    if len(dataset.taxon_sets) != 1:
        raise ValueError("Expecting one set of taxa in %s" % f)
    if len(dataset.tree_lists) != 1:
        raise ValueError("Expecting one tree in %s" % f)
    taxon_set = dataset.taxon_sets[0]
    tree_list = dataset.tree_lists[0]
    number_of_taxon = len(taxon_set)
    branch_counter = 0
    code_list = [StringIO() for i in taxon_set]
    for tree in tree_list:
        encode_splits(tree)
        tree_mask = tree.seed_node.edge.split_bitmask
        assert tree_mask is not None
        tree_tax = set(split_to_list(tree_mask))
        print tree_tax
        split_list = []
        for node in tree.postorder_internal_node_iter():
            if node.parent_node is not None:
                branch_counter +=1
                split_set = set(split_to_list(node.edge.split_bitmask))
                split_list.append(split_set)
        for i,stream in enumerate(code_list):
            if i in tree_tax:
                for split in split_list:
                    if i in split:
                        stream.write('1')
                    else:
                        stream.write('0')
            else:
                stream.write('?'*len(split_list))
                        
    output.write("""#NEXUS
    Begin Data;
    Dimensions ntax = %d nchar = %d;
    Format datatype=standard symbols="01" Missing = ?;
    Matrix \n""" % (number_of_taxon, branch_counter))
    
    escaped_names = []
    max_name_length = 0
    for t in taxon_set:
        escaped_names.append(escape_nexus_token(t.label))
        max_name_length = max(max_name_length, len(escaped_names[-1]))
    fmt = '%%-%ds %%s\n' %(max_name_length)
    for i,stream in enumerate(code_list):
        output.write(fmt %(escaped_names[i],stream.getvalue()))

    output.write(""";
    END;
    """)
        


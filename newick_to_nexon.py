import sys
import json
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
    otu_counter = 0
    node_counter= 0
    edge_counter=0
    for n,tree in enumerate(tree_list):
        otu_list = []
        node_list = [None]
        edge_list = []
        otus = {'otu': otu_list}
        trees = {'tree':[{'node': node_list,
                          'edge': edge_list,
                          '@id': str(n)}
                        ]}

        nexon = {'nexml':{'otus': otus,
                          'trees': trees
                }        }
        for node in tree.preorder_node_iter():
            if len(node.child_nodes()) == 0:
                taxon = node.taxon
                if 'nexon_otu' in taxon.__dict__:
                    o = taxon.nexon_otu
                else:
                    o = str(otu_counter)
                    otu_counter += 1
                    taxon.nexon_otu = o
                d = {'@id':o,
                     '@label':taxon.label,
                     'meta': {'$': int(taxon.label[1:]), 
                              '@property': 'ot:ottolid'
                             }
                    }
                otu_list.append(d)
            edge = node.edge
            parent = edge.tail_node
            if parent is not None:
                if 'nexon_node_id' in node.__dict__:
                    n = node.nexon_node_id
                else:
                    n = 'node'+ str(node_counter )
                    node_counter += 1
                    node.nexon_node_id = n

                if 'nexon_node_id' in parent.__dict__:
                    p = parent.nexon_node_id
                else:
                    p = 'node'+ str(node_counter)
                    node_counter += 1
                    parent.nexon_node_id = p

                if 'nexon_edge_id' in edge.__dict__:
                    e = edge.nexon_edge_id
                else:
                    e = str(edge_counter)
                    edge_counter += 1
                    edge.nexon_edge_id = e
                d = {'@id': n}
                if len(node.child_nodes()) == 0:
                    d['@otu'] = o
                node_list.append(d)
                d = {'@id:': e, '@source': p, '@target': n}
                if node.edge_length is not None:
                    d['@length'] = str(node.edge_length)
                edge_list.append(d)
        d = {'@id': tree.seed_node.nexon_node_id,'@root':'true'}
        node_list[0]= d
        json.dump(nexon, output, sort_keys=True, indent=5)
        output.write('\n')

            


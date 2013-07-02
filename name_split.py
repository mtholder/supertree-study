#!/usr/bin/env python
import dendropy
import sys


if __name__ == '__main__':
	taxa = dendropy.TaxonSet()
	true_tree = dendropy.Tree.get_from_path(sys.argv[1],"Newick",taxon_set = taxa,suppress_internal_node_taxa=False) # true tree (bigtree)
	for node in true_tree.postorder_node_iter():
		try:
			node.taxon.label = node.taxon.label.split('@')[0]
		except:
			pass
	true_tree.write(sys.stdout,schema = "Newick",suppress_rooting = True)
#!/usr/bin/env python
import dendropy
import sys
from dendropy.treesplit import encode_splits

def prune_tree_to_included(tree, included):
	to_exclude = [node for node in tree.leaf_nodes() if node.taxon not in included]
	for nd in to_exclude:
		tree.prune_subtree(nd)

if __name__ == '__main__':
	taxa = dendropy.TaxonSet()
	true_tree = dendropy.Tree.get_from_path(sys.argv[1],"Newick",taxon_set = taxa) # true tree (bigtree)
	mrp_tree = dendropy.Tree.get_from_path(sys.argv[2], "Nexus",taxon_set = taxa) # MRP tree
	sas_tree = dendropy.Tree.get_from_path(sys.argv[3],"Newick",taxon_set = taxa) # SAS tree
	included = set([node.taxon for node in mrp_tree.leaf_nodes()])
	prune_tree_to_included(sas_tree, included)
	prune_tree_to_included(true_tree, included)
	encode_splits(true_tree)
	encode_splits(mrp_tree)
	encode_splits(sas_tree)
	mrp_distance = true_tree.symmetric_difference(mrp_tree)
	sas_distance = true_tree.symmetric_difference(sas_tree)
	print "SAS:"
	sas_tree.write(sys.stdout,schema = "Newick")
	print 
	print "MRP:"
	mrp_tree.write(sys.stdout,schema = "Newick")
	print
	print "true_tree:"
	true_tree.write(sys.stdout,schema = "Newick")
	print "the symmetric distances are: %d for the MRP algorithm and %d for the SAS algorithm"%(mrp_distance,sas_distance)

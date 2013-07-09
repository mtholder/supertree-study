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
	mrp_tree = dendropy.Tree.get_from_path(sys.argv[2],"Nexus",taxon_set = taxa) # MRP tree
	mrp_con = dendropy.Tree.get_from_path(sys.argv[3],"Nexus",taxon_set = taxa) # MRP tree
	sas_tree = dendropy.Tree.get_from_path(sys.argv[4],"Newick",taxon_set = taxa) # SAS tree
	included = set([node.taxon for node in mrp_tree.leaf_nodes()])
	print "number of leaves",len(included)
	prune_tree_to_included(sas_tree, included)
	prune_tree_to_included(true_tree, included)
	encode_splits(true_tree)
	encode_splits(mrp_tree)
	encode_splits(sas_tree)
	mrp_distance = true_tree.false_positives_and_negatives(mrp_tree)
	sas_distance = true_tree.false_positives_and_negatives(sas_tree)
	mrp_to_mrpcon = mrp_tree.false_positives_and_negatives(mrp_con)
	#print sas_tree.as_ascii_plot() 
	#print mrp_tree.as_ascii_plot()
	#print true_tree.as_ascii_plot()
	print "mrp to mrp con",mrp_to_mrpcon[0],mrp_to_mrpcon[1]
	print "true MRP",mrp_distance[0],mrp_distance[1]
	print "true SAS",sas_distance[0],sas_distance[1]
	mrp_to_sas = mrp_tree.false_positives_and_negatives(sas_tree)
	print "MRP SAS ",mrp_to_sas[0],mrp_to_sas[1]
	tree_list = dendropy.TreeList([mrp_tree,sas_tree],taxon_set = taxa)
	con_tree = tree_list.consensus(min_freq=0.99)
	prune_tree_to_included(con_tree,included)
	con_dist = true_tree.false_positives_and_negatives(con_tree)
	#print con_tree.as_ascii_plot()
	print "true con",con_dist[0],con_dist[1]
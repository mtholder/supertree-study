#!/usr/bin/env python
import dendropy
import sys

if __name__ == '__main__':
	taxa = dendropy.TaxonSet()
	true_tree = dendropy.Tree.get_from_path(sys.argv[1],"Newick",taxon_set = taxa) # true tree (bigtree)
	tree1 = dendropy.Tree.get_from_path(sys.argv[2], "Nexus",taxon_set = taxa) # MRP tree
	tree2 = dendropy.Tree.get_from_path(sys.argv[3],"Newick",taxon_set = taxa) # SAS tree
	mrp_distance = true_tree.symmetric_difference(tree1)
	sas_distance = true_tree.symmetric_difference(tree2)
	print "the symmetric distances are: %d for the MRP algorithm and %d for the SAS algorithm"%(mrp_distance,sas_distance)

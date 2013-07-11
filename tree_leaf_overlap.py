#!/usr/bin/env python
import dendropy
import sys
#from dendropy.treesplit import encode_splits

if __name__ == '__main__':
	taxa = dendropy.TaxonSet()
	tree_list = dendropy.TreeList.get_from_path(sys.argv[1],"Newick",taxon_set = taxa) # true tree (bigtree)
	n = len(tree_list)
	print "Tree Number"+'\t'+"Number of leaves"
	for i in range(0,n): 
		#Prints number of leaves on each tree individually
		included = set([node.taxon for node in tree_list[i].leaf_nodes()])
		print "Tree"+str(i+1),'\t\t',len(included)
	print "Number of overlapping leaves:"
	for i in range(0,n-1):
		#Prints number of overlapping branches e.g. t1 to t2, t1 to t2..
		tree_i= tree_list[i]
		if hasattr(tree_i, 'included'):
			included_i = tree_i.included
		else:
			included_i = set([node.taxon for node in tree_list[i].leaf_nodes()])
			tree_i.included = included_i		
		for j in range(i+1,n):
			tree_j= tree_list[j]
			if hasattr(tree_j, 'included'):
				included_j = tree_j.included
			else:
				included_j = set([node.taxon for node in tree_list[j].leaf_nodes()])
				tree_j.included = included_j	 
			#print included_i.intersection(included_j)
			print "Tree"+str(i+1)+" to Tree"+str(j+1),'\t',len(included_i.intersection(included_j))
	
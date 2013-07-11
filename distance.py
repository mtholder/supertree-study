#!/usr/bin/env python
import dendropy
import sys
from dendropy.treesplit import encode_splits
import numpy as np
import matplotlib.pyplot as plt

def prune_tree_to_included(tree, included):
	to_exclude = [node for node in tree.leaf_nodes() if node.taxon not in included]
	for nd in to_exclude:
		tree.prune_subtree(nd)
def bar_graph(column1,column2):
	N = 5
	ind = np.arange(N)  # the x locations for the groups
	width = 0.25       # the width of the bars
	fig = plt.figure()
	ax = fig.add_subplot(111)
	rects1 = ax.bar(ind, column1, width, color='r')#, yerr=menStd)

	#womenStd =   (3, 5, 2, 3, 3)
	rects2 = ax.bar(ind+width, column2, width, color='y')#, yerr=womenStd)

	# add some
	ax.set_ylabel('Value')
	ax.set_title('Method Comparison')
	ax.set_xticks(ind+width)
	plt.yticks(np.arange(0,200,10))
	ax.set_xticklabels( ('MRPtoMRPCON', 'MRP', 'SAS', 'MRPtoSAS', 'CON') )

	ax.legend( (rects1[0], rects2[0]), ('False Positive', 'False Negative') )

	def autolabel(rects):
	    # attach some text labels
	    for rect in rects:
	        height = rect.get_height()
	        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
	                ha='center', va='bottom')

	autolabel(rects1)
	autolabel(rects2)

	#plt.show()
	#fig.savefig( 'distance.png' )
	fig.savefig( 'run'+str(sys.argv[5])+'/'+'distance.png' )

if __name__ == '__main__':
	taxa = dendropy.TaxonSet()
	true_tree = dendropy.Tree.get_from_path(sys.argv[1],"Newick",taxon_set = taxa) # true tree (bigtree)
	mrp_tree = dendropy.Tree.get_from_path(sys.argv[2],"Nexus",taxon_set = taxa) # MRP tree
	mrp_con = dendropy.Tree.get_from_path(sys.argv[3],"Nexus",taxon_set = taxa) # MRP tree
	sas_tree = dendropy.Tree.get_from_path(sys.argv[4],"Newick",taxon_set = taxa) # SAS tree
	included = set([node.taxon for node in mrp_tree.leaf_nodes()])
	#print "number of leaves",len(included)
	#print mrp_tree.leaf_nodes()
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
	col1 = mrp_to_mrpcon[0],mrp_distance[0],sas_distance[0],mrp_to_sas[0],con_dist[0]
	col2 = mrp_to_mrpcon[1],mrp_distance[1],sas_distance[1],mrp_to_sas[1],con_dist[1]
	with open(sys.argv[6], 'a') as the_file:
		the_file.write(str(col1[1])+'\t'+str(col1[2])+'\t'+str(col2[1])+'\t'+str(col2[2])+'\n')
	bar_graph(col1,col2)

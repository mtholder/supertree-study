#!/usr/bin/env python
import dendropy
import sys
from dendropy.treesplit import encode_splits
import numpy as np
import matplotlib.pyplot as plt

def prune_tree_to_included(tree, included):
	'''
	This function prunes down the true tree to to what the MRP tree/TAG tree
	holds.
	'''
	to_exclude = [node for node in tree.leaf_nodes() if node.taxon not in included]
	for nd in to_exclude:
		tree.prune_subtree(nd)

def bar_graph(column1,column2):
	'''
	This function graphs false positives and false negatives
	in the form of a bar graph. 
	'''
	N = 6
	ind = np.arange(N)  # the x locations for the groups
	width = 0.35       # the width of the bars
	fig = plt.figure()
	ax = fig.add_subplot(111)
	rects1 = ax.bar(ind, column1, width, color='r')
	rects2 = ax.bar(ind+width, column2, width, color='y')

	ax.set_ylabel('Number of False Positives and Negatives',fontsize = 18)
	ax.set_title('Method Comparison',fontsize = 18)
	ax.set_xticks(ind+width)
	plt.yticks(np.arange(0,250,25),fontsize = 16)
	ax.patch.set_linewidth(4.0)
	plt.figure(figsize=(170, 170))
	ax.set_xticklabels( ('MRPto\nMRPCON','TAGto\nMRPCON', 'MRP', 'TAG', 'MRPtoTAG', 'CON'),fontsize = 16)
	ax.legend( (rects1[0], rects2[0]), ('False Positive', 'False Negative'),2,fontsize = 15 )

	def autolabel(rects):
	    # attach some text labels
	    for rect in rects:
	        height = rect.get_height()
	        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
	                ha='center', va='bottom',fontsize = 16)

	autolabel(rects1)
	autolabel(rects2)

	#plt.show()
	#fig.savefig( 'distance.png' )
	fig.savefig( 'run'+str(sys.argv[5])+'/'+'distance.png' ) # save as png file
	fig.savefig( 'run'+str(sys.argv[5])+'/'+'distance.eps' ) # save as eps file
if __name__ == '__main__':
	''' 
	The main method is what calculates the distances 
	(e.g. the number of false positives and negatives)
	between MRP and TAG to that of the true tree and the
	consensus tree. Also returns the results in the form of a bar
	graph called from the function 'bar_graph()' written above.  
	'''
	taxa = dendropy.TaxonSet()
	true_tree = dendropy.Tree.get_from_path(sys.argv[1],"Newick",taxon_set = taxa) # true tree (bigtree)
	mrp_tree = dendropy.Tree.get_from_path(sys.argv[2],"Nexus",taxon_set = taxa) # MRP tree
	mrp_con = dendropy.Tree.get_from_path(sys.argv[3],"Nexus",taxon_set = taxa) # MRP tree
	sas_tree = dendropy.Tree.get_from_path(sys.argv[4],"Newick",taxon_set = taxa) # SAS tree
	to_prune = []
	for node in mrp_tree.leaf_nodes():
		if node.taxon.label == 'roottaxon':
			to_prune.append(node)
	assert(len(to_prune) == 1)
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
	sas_to_mrpcon = sas_tree.false_positives_and_negatives(mrp_con)
	#print sas_tree.as_ascii_plot() 
	#print mrp_tree.as_ascii_plot()
	#print true_tree.as_ascii_plot()
	print "mrp to mrp con",mrp_to_mrpcon[0],mrp_to_mrpcon[1]
	print "sas to mrp con",sas_to_mrpcon[0],sas_to_mrpcon[1]
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
	col1 = mrp_to_mrpcon[0],sas_to_mrpcon[0],mrp_distance[0],sas_distance[0],mrp_to_sas[0],con_dist[0]
	col2 = mrp_to_mrpcon[1],sas_to_mrpcon[1],mrp_distance[1],sas_distance[1],mrp_to_sas[1],con_dist[1]
	with open(sys.argv[6], 'a') as the_file:
		'''
		writes MRP and TAG false positves and false negatives to a text file to graph as averages at
		a later time
		'''
		the_file.write(str(col1[2])+'\t'+str(col1[3])+'\t'+str(col2[2])+'\t'+str(col2[3])+'\n')
	bar_graph(col1,col2) # graph results in bar graph function.d

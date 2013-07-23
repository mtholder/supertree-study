#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
def collect_false_pos_and_negs():
	'''
	This function collects all of the false positives and negatives and takes an average of them.
	Will be pipelined at a later date..
	'''
	data = np.loadtxt(sys.argv[1])
	a = data.mean(0)
	mrp_false_pos = a[0]
	sas_false_pos = a[1]
	mrp_false_neg = a[2]
	sas_false_neg = a[3]
	print mrp_false_pos,'\t',sas_false_pos,'\t',mrp_false_neg,'\t',sas_false_neg
	with open(sys.argv[2], 'a') as f:
		f.write(str(mrp_false_pos)+'\t'+str(sas_false_pos)+'\t'+str(mrp_false_neg)+'\t'+str(sas_false_neg)+'\t#for '+str(sys.argv[3])+' trees\n')

def plot_false_pos_and_negs():
	'''
	This function plots an average of false positives and negatives calculated from the function above.
	Will be pipelined at a later time..
	'''
	data = np.loadtxt(sys.argv[1])
	#print data.T[0].split('\t')
	columns = [col for col in data.T]
	x_noise = [5,10,20,40,80,200]
	MRP_false_positive = columns[0] # y1
	SAS_false_positive = columns[1] # y2
	fig1 = plt.figure()
	plt.plot(x_noise, MRP_false_positive,marker='o',color='b', label='MRP')
	plt.plot(x_noise, SAS_false_positive, marker='o', color='r', label='TAG')
	plt.xlabel('Number of Trees',fontsize = 17)
	plt.tick_params(axis='both', which='major', labelsize=17)
	plt.ylabel('Number of False Positives',fontsize = 17)
	plt.title('Number of Trees Vs. Number of False Positives',fontsize = 17)
	plt.legend(loc=1,fontsize = 16,)
	fig1.savefig('false_pos.png')

	fig2 = plt.figure()
	MRP_false_negative = columns[2]
	SAS_false_negative = columns[3]
	plt.plot(x_noise, MRP_false_negative,marker='o',color='b', label='MRP')
	plt.plot(x_noise, SAS_false_negative, marker='o', color='r', label='TAG')
	plt.xlabel('Number of Trees',fontsize = 17)
	plt.ylabel('Number of False Negatives',fontsize = 17)
	plt.tick_params(axis='both', which='major', labelsize=17)
	plt.title('Number of Trees Vs. Number of False Negatives',fontsize = 17)
	plt.legend(loc=1,fontsize = 16,)
	fig2.savefig('false_neg.png')

if __name__ == '__main__':
	'''
	calls functions named above.
	'''
	#collect_false_pos_and_negs()
	plot_false_pos_and_negs()
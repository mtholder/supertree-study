#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
def spr_avg():
	'''
	takes the averages of 10 SPR values (e.g. avg of all 10 runs)
	for each different spr value one at a time. Works manually for now, 
	will be pipelined into the shell script at a later time..
	'''
	data = np.loadtxt(sys.argv[1])
	a = data.mean(0)
	mrp_false_pos = a[0]
	sas_false_pos = a[1]
	mrp_false_neg = a[2]
	sas_false_neg = a[3]
	print mrp_false_pos,'\t',sas_false_pos,'\t',mrp_false_neg,'\t',sas_false_neg
	

	with open(sys.argv[2], 'a') as f:
		f.write(str(mrp_false_pos)+'\t'+str(sas_false_pos)+'\t'+str(mrp_false_neg)+'\t'+str(sas_false_neg)+'\t#for SPR:'+str(sys.argv[3])+'\n')


def spr_plots():
	'''
	Graphs all the spr averages done via the spr_avg() function. Works manually for now.
	Will be pipelined at a later time. 
	'''
	data = np.loadtxt(sys.argv[1])
	#print data.T[0].split('\t')
	columns = [col for col in data.T]
	x_noise = [0,1,2,3,4,5,6]
	MRP_false_positive = columns[0] # y1
	SAS_false_positive = columns[1] # y2
	fig1 = plt.figure()
	plt.plot(x_noise, MRP_false_positive,marker='o',color='b', label='MRP')
	plt.plot(x_noise, SAS_false_positive, marker='o', color='r', label='TAG')
	plt.xlabel('Increase in SPR',fontsize = 17)
	plt.ylabel('Number of False Positives',fontsize = 17)
	plt.title('Increase in SPR Vs. Number of False Positives',fontsize = 17)
	plt.tick_params(axis='both', which='major', labelsize=17)
	plt.legend(loc=2,)
	fig1.savefig('SPR_false_pos.png')
	fig1.savefig('SPR_false_pos.eps')
	fig2 = plt.figure()
	MRP_false_negative = columns[2]
	SAS_false_negative = columns[3]
	plt.plot(x_noise, MRP_false_negative,marker='o',color='b', label='MRP')
	plt.plot(x_noise, SAS_false_negative, marker='o', color='r', label='TAG')
	plt.xlabel('Increase in SPR',fontsize = 17)
	plt.ylabel('Number of False Negatives',fontsize = 17)
	plt.tick_params(axis='both', which='major', labelsize=17)
	plt.title('Increase in SPR Vs. Number of False Negatives',fontsize = 17)
	plt.legend(loc=2,)
	fig2.savefig('SPR_false_neg.png')
	fig2.savefig('SPR_false_neg.eps')


if __name__ == '__main__':
	'''
	calls the functions above individually..
	Pipelined at a later time.
	'''
	#spr_avg()
	spr_plots()
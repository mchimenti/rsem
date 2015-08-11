#! /usr/bin/env python
#Author: Michael S Chimenti
#Date: Aug 11 2015
#Data: ~/iihg/RNA_seq/xue_2015_gene_names

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress

####
#import Bowtie results from normalized counts csv
####


norm_xing_all = pd.read_csv("normcount_xing_all.csv", usecols = range(9))

##average the two replicates

norm_xing_all["cntrl_ave"] = (norm_xing_all.Xing_Ctrl1 + norm_xing_all.Xing_Ctrl2)/2
norm_xing_all["tko_ave"] = (norm_xing_all.Xing_Tko1 + norm_xing_all.Xing_Tko3)/2
norm_xing_all["tlko_ave"] = (norm_xing_all.Xing_TLko1 + norm_xing_all.Xing_TLko2)/2
norm_xing_all["p45ko_ave"] = (norm_xing_all.Xing_p45ko1 + norm_xing_all.Xing_p45ko2)/2

##subset the dataframe

norm_xing_ave = norm_xing_all[['Unnamed: 0','cntrl_ave','tko_ave', 'tlko_ave','p45ko_ave']]

new_colnames = norm_xing_ave.columns.values
new_colnames[0] = 'gene_id'
norm_xing_ave.columns = new_colnames

####
#import rsem expected counts and FPKM/TPM results for each condition 
####

p45ko1_rsem = pd.read_table("mm9_p45ko1.genes.results", sep='\t', usecols=range(7))
p45ko1_final = p45ko1_rsem.merge(norm_xing_ave, how="inner", on="gene_id")


cntrl_rsem = pd.read_table("mm9_cntrl1.genes.results", sep='\t', usecols=range(7))
cntrl_final = cntrl_rsem.merge(norm_xing_ave, how="inner", on="gene_id")


tko1_rsem = pd.read_table("mm9_tko1.genes.results", sep='\t', usecols=range(7))
tko1_final = tko1_rsem.merge(norm_xing_ave, how="inner", on="gene_id")


tlko1_rsem = pd.read_table("mm9_tlko1.genes.results", sep='\t', usecols=range(7))
tlko1_final = tlko1_rsem.merge(norm_xing_ave, how="inner", on="gene_id")

###
#Make the counts vs. FPKM and counts(bowtie) vs. expected counts (RSEM) plots
###

####  p45ko
#
#ax = plt.plot()
#
#p45ko1_final.plot(x="expected_count", y="p45ko_ave", kind="scatter", xlim=(0,100000), ylim=(0,100000), 
#    title="p45ko RSEM Expected Counts vs. Bowtie Norm Counts")
#
#plt.show()   
# 
#slope, intercept, rval, pval, stderr = linregress(p45ko1_final.expected_count, p45ko1_final.p45ko_ave)
#print "p45ko slope, rval:", slope, rval
#
#ax = plt.plot()
#    
#p45ko1_final.plot(x="FPKM", y="p45ko_ave", kind="scatter", xlim=(0,3000), ylim=(0,100000), 
#    title="p45ko RSEM FPKM vs Bowtie Norm Counts")
#
#plt.show()
#
####  cntrl
#
#ax = plt.plot()
#
#cntrl_final.plot(x="expected_count", y="cntrl_ave", kind="scatter", xlim=(0,100000), ylim=(0,100000), 
#    title="Control RSEM Expected Counts vs. Bowtie Norm Counts")
#
#plt.show()   
# 
#slope, intercept, rval, pval, stderr = linregress(p45ko1_final.expected_count, p45ko1_final.p45ko_ave)
#print "cntrl slope, rval:", slope, rval
#
#ax = plt.plot()
#    
#cntrl_final.plot(x="FPKM", y="cntrl_ave", kind="scatter", xlim=(0,3000), ylim=(0,100000), 
#    title="Control RSEM FPKM vs Bowtie Norm Counts")
#
#plt.show()
#
####  tko
#
#ax = plt.plot()
#
#tko1_final.plot(x="expected_count", y="tko_ave", kind="scatter", xlim=(0,100000), ylim=(0,100000), 
#    title="tko RSEM Expected Counts vs. Bowtie Norm Counts")
#
#plt.show()   
# 
#slope, intercept, rval, pval, stderr = linregress(p45ko1_final.expected_count, p45ko1_final.p45ko_ave)
#print "tko slope, rval:", slope, rval
#
#ax = plt.plot()
#    
#tko1_final.plot(x="FPKM", y="tko_ave", kind="scatter", xlim=(0,3000), ylim=(0,100000), 
#    title="tko RSEM FPKM vs Bowtie Norm Counts")
#
#plt.show()

###  tlko
#
#ax = plt.plot()
#
#tlko1_final.plot(x="expected_count", y="tlko_ave", kind="scatter", xlim=(0,100000), ylim=(0,100000), 
#    title="tlko RSEM Expected Counts vs. Bowtie Norm Counts")
#
#plt.show()   
# 
#slope, intercept, rval, pval, stderr = linregress(p45ko1_final.expected_count, p45ko1_final.p45ko_ave)
#print "tlko slope, rval:", slope, rval

ax = plt.plot()
    
tlko1_final.plot(x="FPKM", y="tlko_ave", kind="scatter", title="tlko RSEM FPKM vs Bowtie Norm Counts", xlim=(0,3000), ylim=(0,100000)) 

labels = tlko1_final["gene_id"]


for label, x, y in zip(labels, tlko1_final.FPKM, tlko1_final.tlko_ave):
    if x > 1500 or y > 45000:
        plt.annotate(label, xy=(x,y), xytext=(-5,5), textcoords = 'offset points', fontsize = 10, ha = 'left', va = 'bottom',
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'angle3'))
         

plt.show()
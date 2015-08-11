#Aug 5, 2015  Analysis of DESeq and EBSeq RNASeq results on Xing Xue dataset (Mouse mm9 refseq genes)

import pandas as pd

cpm_xing_all = pd.read_csv("cpm_Xing_all.csv")   #Import DESeq results master table from Tom's CSV file

######
#### Analysis of cntrl vs p45ko for EBSeq/DESeq results
######

cpm_xing_p45ko = cpm_xing_all[cpm_xing_all.ctlVp45 == True]
new_colnames = cpm_xing_p45ko.columns.values
new_colnames[0] = 'Gene'
cpm_xing_p45ko.columns = new_colnames

cntrl_p45ko = pd.read_table("cntrl_p45ko_DEfound", skiprows=0, sep=' ', header = False, names=["Num",'DEfound'])
cntrl_p45ko.drop('Num', axis=1, inplace=True)

sum(cpm_xing_p45ko.Gene.isin(cntrl_p45ko.DEfound))  #number of genes found in common = 114 (44% of DESeq)

cpm_xing_p45ko_overlap = cpm_xing_p45ko[cpm_xing_p45ko.Gene.isin(cntrl_p45ko.DEfound)]
cpm_xing_p45ko_overlap.to_csv("cpm_xing_cntrlVp45ko_DESeq_EBSeq_overlap_genes.csv", sep = ',')


######
#### Analysis of cntrl vs tlko data for EBSeq/DESeq
######

cpm_xing_tlko = cpm_xing_all[cpm_xing_all.ctlVtlko == True]   #DESeq analysis results
new_colnames = cpm_xing_tlko.columns.values
new_colnames[0] = 'Gene'
cpm_xing_tlko.columns = new_colnames

cntrl_tlko = pd.read_table("cntrl_tlko_DEfound", skiprows=0, sep=' ', header = False, names=["Num",'DEfound']) #EBSeq
cntrl_tlko.drop('Num', axis=1, inplace=True)

sum(cpm_xing_tlko.Gene.isin(cntrl_tlko.DEfound))   #number of genes in common = 2779 (90% of DESeq)

cpm_xing_tlko_overlap = cpm_xing_tlko[cpm_xing_tlko.Gene.isin(cntrl_tlko.DEfound)]
cpm_xing_tlko_overlap.to_csv("cpm_xing_cntrlVtlko_DESeq_EBSeq_overlap_genes.csv", sep = ',')

######
#### Analysis of cntrl vs tko data for EBSeq/DEseq
######

cpm_xing_tko = cpm_xing_all[cpm_xing_all.ctlVtko == True]   #DESeq analysis results
new_colnames = cpm_xing_tko.columns.values
new_colnames[0] = 'Gene'
cpm_xing_tko.columns = new_colnames

cntrl_tko = pd.read_table("cntrl_tko_DEfound", skiprows=0, sep=' ', header = False, names=["Num",'DEfound']) #EBSeq
cntrl_tko.drop('Num', axis=1, inplace=True)

sum(cpm_xing_tko.Gene.isin(cntrl_tko.DEfound))   #number of genes in common = 1853 (86% of DESeq)

cpm_xing_tko_overlap = cpm_xing_tko[cpm_xing_tko.Gene.isin(cntrl_tko.DEfound)]
cpm_xing_tko_overlap.to_csv("cpm_xing_cntrlVtko_DESeq_EBSeq_overlap_genes.csv", sep = ',')

#####
###Analysis of p45 vs tko data for EBSeq/DESeq
#####

cpm_xing_tko_v_p45 = cpm_xing_all[cpm_xing_all.tkoVp45 == True]   #DESeq analysis results
new_colnames = cpm_xing_tko_v_p45.columns.values
new_colnames[0] = 'Gene'
cpm_xing_tko_v_p45.columns = new_colnames

tko_v_p45 = pd.read_table("p45ko_tko_DEfound", skiprows=0, sep=' ', header = False, names=["Num",'DEfound']) #EBSeq
tko_v_p45.drop('Num', axis=1, inplace=True)

print "The number of genes in common between methods on p45 v tko is ", sum(cpm_xing_tko_v_p45.Gene.isin(tko_v_p45.DEfound))   #number of genes in common = 1405

cpm_xing_tko_v_p45_overlap = cpm_xing_tko[cpm_xing_tko.Gene.isin(cntrl_tko.DEfound)]
cpm_xing_tko_v_p45_overlap.to_csv("cpm_xing_tkoVp45_DESeq_EBSeq_overlap_genes.csv", sep = ',')

#####
###Analysis of tko vs tlko data for EBSeq/DESeq
#####

cpm_xing_tko_v_tlko = cpm_xing_all[cpm_xing_all.tlkoVtko == True]   #DESeq analysis results
new_colnames = cpm_xing_tko_v_tlko.columns.values
new_colnames[0] = 'Gene'
cpm_xing_tko_v_tlko.columns = new_colnames

tko_v_tlko = pd.read_table("tko_tlko_DEfound", skiprows=0, sep=' ', header = False, names=["Num",'DEfound']) #EBSeq
tko_v_tlko.drop('Num', axis=1, inplace=True)

print "The number of genes in common between methods on tko v tlko is ", sum(cpm_xing_tko_v_tlko.Gene.isin(tko_v_tlko.DEfound))   #number of genes in common = 3258

cpm_xing_tko_v_tlko_overlap = cpm_xing_tko[cpm_xing_tko.Gene.isin(cntrl_tko.DEfound)]
cpm_xing_tko_v_tlko_overlap.to_csv("cpm_xing_tkoVtlko_DESeq_EBSeq_overlap_genes.csv", sep = ',')



import pandas as pd

def strip_comma(oldfile, newfile):
    """strip commas out of UCSC genes to IDs file"""
    with open(newfile, 'w') as outfile, open(oldfile, 'r') as infile:
        for line in infile:
            genename = line.split()[0]
            geneid = line.split()[1].strip(',')
            newline = genename + " " + geneid + "\n"
            outfile.write(newline)
        
        
def gene_df(filename):
    """Read in the genes to UCSC Ids file as a dataframe"""
    gene_df = pd.read_table('mm9_RefSeq_to_UCSCgeneid_nocomma.txt', sep=' ', names=['symbol','id'], header=0)
    gene_df = gene_df.set_index('id')
    return gene_df
    

def gene_names(filename, gene_df):
    """Read in genes.results file and join on transcript_ids"""
    df = pd.read_table(filename, sep='\t', header=0)
    df['transcript_id(s)'] = df['transcript_id(s)'].map(lambda r: r.split(',')[0])  #only keep one ID to join on
    df = df.join(genes, on='transcript_id(s)')   #join on genes index and "transcript_id" column
    df['gene_id'] = df['symbol']
    df = df.drop('symbol', axis=1)
    df = df[df.gene_id.notnull()]     #drop NAs (~13)
    df = df[df.gene_id.duplicated == 0]  #drop duplicate gene symbols  (~900)
    df.to_csv(filename + '.mod', sep='\t', index=False, header=True)


#####

if __name__ == '__main__':
    
    strip_comma('mm9_RefSeq_to_UCSCgeneid.txt', 'mm9_RefSeq_to_UCSCgeneid_nocomma.txt')
    
    genes = gene_df('mm9_RefSeq_to_UCSCgeneid_nocomma.txt')
    
    files = ['mouse_cntrl_1.genes.results','mouse_cntrl_2.genes.results', 'mouse_tko_1.genes.results','mouse_tko3.genes.results',
        'mouse_tlko1.genes.results', 'mouse_tlko2.genes.results','mouse_p45ko1.genes.results','mouse_p45ko2.genes.results']
    
    for filename in files:
        gene_names(filename, genes)     
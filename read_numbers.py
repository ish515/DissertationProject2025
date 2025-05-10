# Import necessary packages
import itertools as it
import numpy as np
import pandas as pd
from Bio import SeqIO

# Import the FASTA sequence
for seq_record in SeqIO.parse("NC_000913_new.fasta", "fasta"):
	sequence = seq_record.seq

# Import the ribosome profiling data
colnames = ['nucleotide', 'read_num']

read_numbers = pd.read_csv('Mohammad_2016.dens', sep=r"\s{2,}", names=colnames, engine='python', header=None)

# Import the start and stop codon positions
colnames = ['Start', 'Stop', 'FR', 'name']

genes = pd.read_csv('ecoli_CDS_name.dat', sep='\s+', names=colnames, header=None)

# Import the gene copy numbers and lengths
colnames = ['name', 'copy_num', 'length']

gene_copy_num = pd.read_csv('gene_counts.dat', sep='\s+', names=colnames, header=None)

# Merge dataframes by gene name
df_merged = genes.merge(gene_copy_num[['name', 'copy_num', 'length']])

# For ease, rename df_merged and use throughout
genes = df_merged

# Divide code into forward and reverse genes
forward = genes[genes["FR"] == 1]

reverse = genes[genes["FR"] == -1]

# Set codon of interest
coi = 'CTG'

# Ouput read numbers by codon
codon_read_numbers = []
reads_over_length =[]

for i in range(len(genes)):
    if genes.iloc[i]['FR'] == 1:
        
        # For forward genes
        ## Pull out the start and stop positions for the gene
        start = int(forward.Start[i])
        stop = int(forward.Stop[i])
        
        ## Subset the RP data by gene position
        gene_read_numbers = read_numbers.read_num[start-1:stop]
        
        ## Take the max of the group of three (one codon)
        codons = [gene_read_numbers[j:j+3] for j in range(0, len(gene_read_numbers), 3)]
        
        codons_max = list(map(max, codons))
        
        ## Subset FASTA sequence by start and stop positions
        gene_seq = sequence[start-1:stop]
        
        ## Group together each codon
        gene_seq_codons = [gene_seq[j:j+3] for j in range(0, len(gene_seq), 3)]
        
        ## Subset codon list by codon of interest
        fil = [gene_seq_codons[n] == coi for n in range(0, len(gene_seq_codons))]
        codon_indices = list(it.compress(range(len(fil)), fil))
        
        ## Subset read numbers by codon indices
        codon_reads = [codons_max[j] for j in codon_indices]
        
        ## Append to full data
        codon_read_numbers.append(codon_reads)
        
        ## Divide read numbers by gene length
        reads_by_len = [codons_max[j] / int(genes.length[i]) for j in range(len(codons_max))]
        
        ## Append to full data
        reads_over_length.append(reads_by_len)
        
    else:               
        # For reverse genes
        ## Pull out the start and stop positions for the gene
        start = int(reverse.Start[i])
        stop = int(reverse.Stop[i])
        
        ## Reverse the order of the ribosome density values
        reverse_read_numbers = list(reversed(read_numbers.read_num))
        
        ## Sum ribosome density of each codon of the gene for plus1 strand
        ### Subset the RP data by gene position
        read_number_list = reverse_read_numbers[start-1:stop]
        
        ## Take the max read number of the group of three (one codon)
        codons = [read_number_list[j:j+3] for j in range(0, len(read_number_list), 3)]
        
        codons_max = list(map(max, codons))
        
        ## Subset FASTA sequence by start and stop positions
        gene_seq = sequence[start-1:stop]

        ## Find the reverse complement
        reverse_comp_gene_seq = gene_seq.reverse_complement()

        ## Group together each codon
        gene_seq_codons = [reverse_comp_gene_seq[j:j+3] for j in range(0, len(reverse_comp_gene_seq), 3)]
        
        ## Concatenate strings of nucleotides into sets of 3
        for r in range(len(gene_seq_codons)):
                gene_seq_codons[r] = gene_seq_codons[r][0] + gene_seq_codons[r][1] + gene_seq_codons[r][2]
        
        ## Subset codon list by codon of interest
        fil = [gene_seq_codons[n] == coi for n in range(0, len(gene_seq_codons))]
        codon_indices = list(it.compress(range(len(fil)), fil))
        
        ## Subset read numbers by codon indices
        codon_reads = [codons_max[j] for j in codon_indices]
        
        ## Append to full data
        codon_read_numbers.append(codon_reads)
        
        ## Divide read numbers by gene length
        reads_by_len = [codons_max[j] / int(genes.length[i]) for j in range(len(codons_max))]
        
        ## Append to full data
        reads_over_length.append(reads_by_len)
        
# Flatten the list of lists into one list
codon_reads = list(it.chain(*codon_read_numbers))
codon_reads_bylen = list(it.chain(*reads_over_length))

# Save the data to file
file_name = 'codon_'+coi+'_reads.csv'
np.savetxt(file_name, codon_reads, delimiter=', ', fmt='% s')

file_name = 'codon_'+coi+'_reads_bylen.csv'
np.savetxt(file_name, codon_reads_bylen, delimiter=', ', fmt='% s')

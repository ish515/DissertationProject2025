# Import necessary packages
import matplotlib.pyplot as plt
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

# Calculate T for every gene
## Calulate number of codons in every gene
num_codons = genes.length/3

## T = num. codons / 18
T = num_codons / 18

########## LOOKING AT A GENE OF INTEREST ##########

## State the name of the gene of interest
name = 'yaaA'

# Print the gene info
print(genes[genes["name"] == name])

# Change index
i = 5

############ For forward genes ##########
#### Pull out the start and stop positions for the gene
##start = int(forward.Start[i])
##stop = int(forward.Stop[i])
##
#### Subset FASTA sequence by start and stop positions
##gene_seq = sequence[start-1:stop]
##
#### Group together each codon
##gene_seq_codons = [gene_seq[j:j+3] for j in range(0, len(gene_seq), 3)]
##
#### Subset read data by gene position
##gene_read_numbers = read_numbers[start-1:stop]
##
##read_number_list = []
##
##for k in range(len(gene_read_numbers)):
##        read_number_list.append(float(gene_read_numbers.iloc[k][1]))
#########################################
                
########## For reverse genes ##########
## Pull out the start and stop positions for the gene
start = int(reverse.Start[i])
stop = int(reverse.Stop[i])

## Subset FASTA sequence by start and stop positions
gene_seq = sequence[start-1:stop]

## Find the reverse complement
reverse_comp_gene_seq = gene_seq.reverse_complement()

## Group together each codon
gene_seq_codons = [reverse_comp_gene_seq[j:j+3] for j in range(0, len(reverse_comp_gene_seq), 3)]

## Reverse the order of the numbers of reads
reverse_read_numbers = list(reversed(read_numbers.read_num))

## Sum ribosome density of each codon of the gene
### Subset the RP data by gene position
gene_read_numbers = reverse_read_numbers[start-1:stop]

## Pull out just the read values
read_number_list = []

for k in range(len(gene_read_numbers)):
        read_number_list.append(float(gene_read_numbers[k]))
#######################################

########## For all genes ##########
## Take the max of the group of three (one codon)
codons = [read_number_list[i:i+3] for i in range(0, len(read_number_list), 3)]

codons_max = list(map(max, codons))

# Calculate average translation time of each codon in the gene
## Sum read numbers per codon across the entire gene
summed_max = sum(codons_max)

## Divide each codon_max value by the sum
codons_max_over_sum = [x / summed_max for x in codons_max]

## Multiply by T for that gene
tau = [x * int(T[i]) for x in codons_max_over_sum]

# Save the data to file
file_name = 'tau_gene_'+str(name)+'_data.csv'
np.savetxt(file_name,
           tau,
           delimiter=", ",
           fmt='% s')

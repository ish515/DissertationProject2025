# Import necessary packages
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats

# Import the map file containing the sequence of fragments
colnames = ['len', 'FR', 'ID', 'pos', 'sequence', 'block', 'zero', 'base']

ecoli_map = pd.read_csv('e_coli_1.map', sep='\s+', names=colnames, header=None)

# Import the full E.coli nucleotide sequence
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

## Merge dataframes by gene name
df_merged = genes.merge(gene_copy_num[['name', 'copy_num', 'length']])

## For ease, rename df_merged and use throughout
genes = df_merged

# Split the map dataset into forward and reverse fragments
forward = ecoli_map[ecoli_map["FR"] == '+']
reverse = ecoli_map[ecoli_map["FR"] == '-']

# Define a function to sum one column of a list of lists
def sum_column(nums, C):

    # Calculate the sum of the specified column 'C' using a generator expression
    result = sum(row[C] for row in nums)
    return result

###### FORWARD FRAGMENTS ######
# Subset the final 2 nucleotides from each fragment and first 2 from the next fragment
## Output the position of the final nucleotide of each fragment
Y_index_f = []

for i in range(len(forward)):
               Y_index_f.append(forward.len.iloc[i] + forward.pos.iloc[i])

## Remove Y_index values > len(sequence)
Y_index_f_inseq = [i for i in Y_index_f if i < len(sequence)-2]

###### REVERSE FRAGMENTS ######
# Subset the final 2 nucleotides from each fragment and first 2 from the next fragment
## Output the position of the final nucleotide of each fragment
Y_index_r = []

for i in range(len(reverse)):
               Y_index_r.append(reverse.pos.iloc[i])

## Remove Y_index values > len(sequence)
Y_index_r_inseq = [i for i in Y_index_r if i < len(sequence)-2]

## Subset the two nucleotides before cleavage
beforef = []
x1f = []
Yf = []
for x in Y_index_f_inseq:
        x1 = sequence[x-2]
        Y = sequence[x-1]
        beforef.append(x1)
        beforef.append(Y)
        x1f.append(x1)
        Yf.append(Y)

beforer = []
x1r = []
Yr = []
for x in Y_index_r_inseq:
        x1 = sequence[x+1]
        Y = sequence[x]
        x1str = ''.join(x1)
        Ystr = ''.join(Y)
        complementx1 = Seq(str(x1str)).complement()
        complementY = Seq(str(Ystr)).complement()
        beforer.append(str(complementx1))
        beforer.append(str(complementY))
        x1r.append(str(complementx1))
        Yr.append(str(complementY))

## Append the two lists together
before = beforef + beforer
x1 = x1f + x1r
Y = Yf + Yr

## Divide the lists up into pairs of nucleotides
before_pairs = [before[m:m+2] for m in range(0, len(before), 2)]

## Join nucleotides into single string
before_pairs_tg = [''.join(before_pairs[x]) for x in range(len(before_pairs))]

## Subset the two nucleotides after cleavage
afterf = []
Zf = []
x2f = []
for x in Y_index_f_inseq:
        Z = sequence[x]
        x2 = sequence[x+1]
        afterf.append(Z)
        afterf.append(x2)
        Zf.append(Z)
        x2f.append(x2)

afterr = []
Zr = []
x2r = []
for x in Y_index_r_inseq:
        Z = sequence[x-1]
        x2 = sequence[x-2]
        Zstr = ''.join(Z)
        x2str = ''.join(x2)
        complementZ = Seq(str(Zstr)).complement()
        complementx2 = Seq(str(x2str)).complement()
        afterr.append(str(complementZ))
        afterr.append(str(complementx2))
        Zr.append(str(complementZ))
        x2r.append(str(complementx2))

## Append the two lists together
after = afterf + afterr
Z = Zf + Zr
x2 = x2f + x2r

## Divide the lists up into pairs of nucleotides
after_pairs = [after[m:m+2] for m in range(0, len(after), 2)]

## Join nucleotides into single string
after_pairs_tg = [''.join(after_pairs[x]) for x in range(len(after_pairs))]

# Count the number of occurrences of each possible sequence of nucleotides
num_occ_before = [(j, before_pairs_tg.count(j)) for j in set(before_pairs_tg)]
num_occ_after = [(j, after_pairs_tg.count(j)) for j in set(after_pairs_tg)]
num_occ_x1 = [(j, x1.count(j)) for j in set(x1)]
num_occ_Y = [(j, Y.count(j)) for j in set(Y)]
num_occ_Z = [(j, Z.count(j)) for j in set(Z)]
num_occ_x2 = [(j, x2.count(j)) for j in set(x2)]

# Calculate percentage of each type of nucleotide pair before and after cleavage
## Before
before_total = sum_column(num_occ_before, 1)
before_percent = [num_occ_before[i][1]/before_total for i in range(len(num_occ_before))]

## After
after_total =  sum_column(num_occ_after, 1)
after_percent = [num_occ_after[i][1]/after_total for i in range(len(num_occ_after))]

## x1
x1_total = sum_column(num_occ_x1, 1)
x1_percent = [num_occ_x1[i][1]/x1_total for i in range(len(num_occ_x1))]

## Y
Y_total = sum_column(num_occ_Y, 1)
Y_percent = [num_occ_Y[i][1]/Y_total for i in range(len(num_occ_Y))]

## Z
Z_total = sum_column(num_occ_Z, 1)
Z_percent = [num_occ_Z[i][1]/Z_total for i in range(len(num_occ_Z))]

## x2
x2_total = sum_column(num_occ_x2, 1)
x2_percent = [num_occ_x2[i][1]/x2_total for i in range(len(num_occ_x2))]

# Append percentages to number of occurrences lists
before_results = [[*i, j*100] for i,j in zip(num_occ_before, before_percent)]
after_results = [[*i, j*100] for i,j in zip(num_occ_after, after_percent)]
x1_results = [[*i, j*100] for i,j in zip(num_occ_x1, x1_percent)]
Y_results = [[*i, j*100] for i,j in zip(num_occ_Y, Y_percent)]
Z_results = [[*i, j*100] for i,j in zip(num_occ_Z, Z_percent)]
x2_results = [[*i, j*100] for i,j in zip(num_occ_x2, x2_percent)]

# Save the results to file
## Before
file_name = 'before_nucleotides.csv'
np.savetxt(file_name, before_results, delimiter=', ', fmt='% s')

## After
file_name = 'after_nucleotides.csv'
np.savetxt(file_name, after_results, delimiter=', ', fmt='% s')

## x1
file_name = 'x1_nucleotides.csv'
np.savetxt(file_name, x1_results, delimiter=', ', fmt='% s')

## Y
file_name = 'Y_nucleotides.csv'
np.savetxt(file_name, Y_results, delimiter=', ', fmt='% s')

## Z
file_name = 'Z_nucleotides.csv'
np.savetxt(file_name, Z_results, delimiter=', ', fmt='% s')

## x2
file_name = 'x2_nucleotides.csv'
np.savetxt(file_name, x2_results, delimiter=', ', fmt='% s')

# Divide data into forward and reverse genes
forward_genes = genes[genes["FR"] == 1]

reverse_genes = genes[genes["FR"] == -1]

# Calculate occurrence of each dinucleotide sequence within all genes
pairs = []
for i in range(len(genes)):
        if genes.iloc[i]["FR"] == 1:
                ## Pull out the start and stop positions for the gene
                start = int(forward_genes.Start[i])
                stop = int(start-1+forward_genes.length[i])
                
                ## Subset the sequence by start and stop positions
                gene_seq = sequence[start-1:stop]
                
                ## Create an empty vector to hold the pairs of nucleotides of this gene
                gene_pairs = []
                
                ## Create a list of each possible pair of nucleotides
                for j in range(len(gene_seq)-1):
                        first = gene_seq[j]
                        second = gene_seq[j+1]
                        gene_pairs.append(first)
                        gene_pairs.append(second)
                
                ## Append list to full list
                pairs.append(gene_pairs)
                
        else:
                ## Pull out the start and stop positions for the gene
                start = int(reverse_genes.Start[i])
                stop = int(start-1+reverse_genes.length[i])
                
                ## Subset the sequence by start and stop positions
                gene_seq = sequence[start-1:stop]

                ## Reverse the gene sequence
                rev_gene_seq = list(reversed(gene_seq))
                
                ## Create an empty vector to hold the pairs of nucleotides of this gene
                gene_pairs = []
                
                ## Create a list of each possible pair of nucleotides
                for j in range(len(rev_gene_seq)-1):
                        first = rev_gene_seq[j]
                        second = rev_gene_seq[j+1]
                        gene_pairs.append(first)
                        gene_pairs.append(second)
                
                ## Append list to full list
                pairs.append(gene_pairs)

## Create a flat list from the list of lists
flat_pairs = [
        x
        for xs in pairs
        for x in xs
]

## Divide the lists up into pairs of nucleotides
pairs_paired = [flat_pairs[m:m+2] for m in range(0, len(flat_pairs), 2)]

## Join nucleotides into single string
pairs_tg = [''.join(pairs_paired[x]) for x in range(len(pairs_paired))]

# Count the number of occurrences of each possible sequence of nucleotides
num_occ_pairs = [(j, pairs_tg.count(j)) for j in set(pairs_tg)]

# Calculate percentage of each type of nucleotide pair before and after cleavage
pairs_total = sum_column(num_occ_pairs, 1)
pairs_percent = [num_occ_pairs[i][1]/pairs_total for i in range(len(num_occ_pairs))]

# Append percentages to number of occurrences lists
pairs_results = [[*i, j*100] for i,j in zip(num_occ_pairs, pairs_percent)]

# Save the results to file
file_name = 'pairs_percentages.csv'
np.savetxt(file_name, pairs_results, delimiter=', ', fmt='% s')

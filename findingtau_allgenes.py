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

# Divide data into forward and reverse genes
forward = genes[genes["FR"] == 1]

reverse = genes[genes["FR"] == -1]

# Calculate T for every gene
## Calulate number of codons in every gene
num_codons = genes.length/3

## T = num. codons / 18
T = num_codons / 18

########## LOOKING AT A CODON OF INTEREST ##########

# State codon of interest
coi = 'AGA'

max_tau = []
for i in range(len(genes)):
       if genes.iloc[i]["FR"] == 1:
               # Pull out the start and stop positions for the gene
               start = int(forward.Start[i])
               stop = int(start-1+forward.length[i])
               
               # Subset the RP data by gene position
               gene_read_numbers = read_numbers.read_num[start-1:stop]

               # Take the max of the group of three (one codon)
               codons = [gene_read_numbers[j:j+3] for j in range(0, len(gene_read_numbers), 3)]
               
               codons_max = list(map(max, codons))
               
               # Output the gene sequence
               gene_seq = sequence[start-1:stop]
               
               # Group together each codon
               gene_seq_codons = [gene_seq[m:m+3] for m in range(0, len(gene_seq), 3)]
               
               # Subset codon list by codon of interest
               fil = [gene_seq_codons[n] == coi for n in range(0, len(gene_seq_codons))]
               codon_indices = list(it.compress(range(len(fil)), fil))
               
               # Calculate average translation time of each codon in the gene
               ## Sum read numbers per codon across the entire gene
               summed_max = sum(codons_max)
               
               # Set tau to zero if summed_max is zero
               if summed_max == 0:
                       tau = 0
               else:
                       # Divide each codon_max value by the sum
                       codons_max_over_sum = [x / summed_max for x in codons_max]
                       
                       # Multiply by T for that gene
                       tau = [y * int(T[i]) for y in codons_max_over_sum]
                       
                       # Subset tau values by codon indices
                       codon_tau = [tau[p] for p in codon_indices]
                       
                       # Append all tau values together
                       max_tau.append(codon_tau)
       else:
               # Pull out the start and stop positions for the gene
               start = int(reverse.Start[i])
               stop = int(start-1+reverse.length[i])
               
               # Reverse the order of the ribosome density values
               reverse_read_numbers = list(reversed(read_numbers.read_num))
               
               # Sum ribosome density of each codon of the gene for plus1 strand
               ## Subset the RP data by gene position
               read_number_list = reverse_read_numbers[start-1:stop]
               
               # Take the max read number of the group of three (one codon)
               codons = [read_number_list[j:j+3] for j in range(0, len(read_number_list), 3)]
               
               codons_max = list(map(max, codons))
               
               # Output the gene sequence
               gene_seq = sequence[start-1:stop]
               
               # Find the reverse complement
               rev_complement = gene_seq.reverse_complement()
               
               # Group together each codon
               gene_seq_codons = [rev_complement[m:m+3] for m in range(0, len(rev_complement), 3)]
               
               # Concatenate strings of nucleotides into sets of 3
               for r in range(len(gene_seq_codons)):
                       gene_seq_codons[r] = gene_seq_codons[r][0] + gene_seq_codons[r][1] + gene_seq_codons[r][2]
               
               # Subset codon list by codon of interest
               fil = [gene_seq_codons[n] == coi for n in range(0, len(gene_seq_codons))]
               codon_indices = list(it.compress(range(len(fil)), fil))
               
               # Calculate average translation time of each codon in the gene
               ## Sum read numbers per codon across the entire gene
               summed_max = sum(codons_max)
               
               # Set tau to zero if summed_max is zero
               if summed_max == 0:
                       tau = 0
               else:
                       # Divide each codon_max value by the sum
                       codons_max_over_sum = [x / summed_max for x in codons_max]
                       
                       # Multiply by T for that gene
                       tau = [y * int(T[i]) for y in codons_max_over_sum]
                       
                       # Subset tau values by codon indices
                       codon_tau = [tau[p] for p in codon_indices]
                       
                       # Append all tau values together
                       max_tau.append(codon_tau)

# Flatten the list of lists into one list
codon_tau = list(it.chain(*max_tau))

# Save the data to file
file_name = 'codon_'+coi+'_tau_all.csv'
np.savetxt(file_name, codon_tau, delimiter=', ', fmt='% s')

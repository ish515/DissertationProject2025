# Statistical Analysis of Dwell Times from Ribosome Profiling Data in _E.coli_
The following workflow is used to analyse _E.coli_ ribosome profiling data to estimate codon dwell times. This repository contains alternative pipelines that estimate dwell times by fitting a Poisson distribution to either codon translation times or ribosome read numbers.

## Calculating Percentage of Each Nucleotide at Cleavage Site
1. Python script `nuclease_bias.py` calculates the percentages of nucleotide pairs present at the 5' and 3' sides of the nuclease cleavage site. It also calculates the total percentage of each nucleotide pair across all _E.coli_ genes.
2. R script `cleavagestats.R` visualises the results as a set of tables.

## Calculating Average Translation Time of Codons in a Gene
1. Python script `findingtau.py` takes as input the name of a gene of interest and calculates the average translation time at each codon across that gene.
2. R script `plottingtau.R` visualises the result as a plot of codon position by average translation time.

## Calculating Average Codon Dwell Time by Fitting a Poisson Rate Parameter to Read Number Data
1. Python script `read_numbers.py` takes as input the nucleotide sequence of a codon of interest and outputs the read number values for each instance of that codon across all _E.coli_ genes.
2. R script `findinglambda_reads.R` fits the data to a Poisson distribution to estimate a rate parameter representing the average number of reads for each codon. The script also visualises the results as a series of plots and applies the Kolmogorov-Smirnov test to compare the cumulative distribution functions of the observed data and the fitted Poisson distribution.

## Calculating Average Codon Dwell Time by Fitting a Poisson Rate Parameter to Translation Time Data
1. Python script `findingtau_allgenes.py` takes as input the nucleotide sequence of a codon of interest and calculates the translation time for every instance of that codon across all _E.coli_ genes.
2. R script `findinglambda_time.R` fits the data to a Poisson distribution to estimate a rate parameter representing the average translation time of each codon. The script also visualises the results as a series of plots and applies the Kolmogorov-Smirnov test to compare the cumulative distribution functions of the observed data and the fitted Poisson distribution.

# Author
Isabella Hodgson

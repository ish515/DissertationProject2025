# Load necessary packages
library(MASS)
library(actuar)
library(gt)

# Colourblind-friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Create a vector containing all 64 codons
codon = c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
          "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
          "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
          "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
          "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG",
          "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
          "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
          "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")

# Create vectors to store dwell times, lambda estimate, mean, and SD
dwell_times = c()
lambda_values = c()
mean_values = c()
var_values = c()
length_values = c()
ks_stat = c()
p_values = c()

# Iterate through all 64 codons
for (i in 1:length(codon)){
  # Set codon of interest
  coi = codon[i]
  
  # Import the reads data 
  file_name = paste('codon_', coi, '_reads.csv', sep='')
  read_data = read.csv(file_name)
  
  # Output length of data
  length_values[i] <- length(read_data[,1])
  
  # Save mean and SD of data
  mean_values[i] <- mean(read_data[,1])
  var_values[i] <- var(read_data[,1])
  
  # Calculate the CDF of the data
  read_cdf = ecdf(read_data[,1])
  
  # Fit a Poisson distribution to the data
  fit <- fitdistr(read_data[,1], "poisson")
  
  # Calculate the CDF of this Poisson distribution
  read_seq <- seq(0, max(read_data[,1])-1)
  cdf_poisson <- ppois(read_seq, lambda=fit$estimate)
  
  # Save the following plot
  file_name = paste('codon_', coi, '_cdf_reads.png', sep='')
  png(file_name)
  par(mar=c(5, 5, 4, 2))
  
  # Plot the CDF of the data and the fitted Poisson distribution
  plot(read_seq, ppois(read_seq, lambda=fit$estimate), type="l", col=cbPalette[6], lwd=2,
       xlab="Average Number of Reads", ylab="Cumulative Probability",
       cex.axis=1.5, cex.lab=1.5)
  lines(read_cdf(0:max(read_data[,1])), col="black", lwd=2)
  
  # Add a legend
  legend("bottomright", title='CDF', legend=c("Poisson Fit", "Data"), 
         col=c(cbPalette[6], "black"), lwd=2, bty="n", inset=c(0.05, 0.05),
         cex=1.5)
  
  dev.off()
  
  # Plot their densities
  plot(density(read_data[,1]), col=cbPalette[2], lwd=2,
       xlab="Average Number of Reads", ylab="Probability Density", main="")
  lines(dpois(read_seq, lambda=fit$estimate), col=cbPalette[6], lwd=2)
  
  # Add a legend
  legend("topright", title='Density', legend=c("Data", "Poisson Fit"), col=c(cbPalette[2], cbPalette[6]), lwd=2, bty="n")
  
  dev.off()
  
  # Compare the two CDFs using the Kolmogorov-Smirnov test
  ks_test <- ks.test(read_cdf(0:max(read_data[,1])), cdf_poisson)
  ks_stat[i] <- ks_test$statistic
  p_values[i] <- ks_test$p.value
  
  # Save the following plot
  file_name = paste('codon_', coi, '_histogramwithdensity.png', sep='')
  png(file_name)
  par(mar=c(5, 5, 4, 2))
  
  # Plot a histogram of the data
  hist(read_data[,1], prob=TRUE, breaks=10000, xlim=c(0,80), col="grey",
       xlab="Average Number of Reads", ylab="Probability Density", main="",
       cex.axis=1.5, cex.lab=1.5)

  # Add density estimates to the histogram
  lines(density(read_data[,1]), col=cbPalette[2], lwd=2)
  lines(dpois(read_seq, lambda=fit$estimate), col=cbPalette[6], lwd=2)

  # Add a legend
  legend("topright", title='Density', legend=c("Data", "Poisson Fit"),
         col=c(cbPalette[2], cbPalette[6]), lwd=2, bty="n",
         cex=1.5)

  dev.off()
  
  # Save lambda
  lambda_values[i] <- fit$estimate
}

# Save the following plot
png("mean_vs_sd_reads.png")
par(mar=c(5, 5, 4, 2))

# Plot mean values against SD values
plot(mean_values, var_values, xlab="Mean (no. of reads)", ylab=expression("Variance (no. of reads"^2*")"), 
     xlim=c(0,20), ylim=c(0,15000), cex.axis=1.5, cex.lab=1.5)

dev.off()

# Create a vector containing all 64 codons but swap T for U
codon = c("UUU", "UUC", "UUA", "UUG", "UCU", "UCC", "UCA", "UCG",
          "UAU", "UAC", "UAA", "UAG", "UGU", "UGC", "UGA", "UGG",
          "CUU", "CUC", "CUA", "CUG", "CCU", "CCC", "CCA", "CCG",
          "CAU", "CAC", "CAA", "CAG", "CGU", "CGC", "CGA", "CGG",
          "AUU", "AUC", "AUA", "AUG", "ACU", "ACC", "ACA", "ACG",
          "AAU", "AAC", "AAA", "AAG", "AGU", "AGC", "AGA", "AGG",
          "GUU", "GUC", "GUA", "GUG", "GCU", "GCC", "GCA", "GCG",
          "GAU", "GAC", "GAA", "GAG", "GGU", "GGC", "GGA", "GGG")

# Create a table of lambda values
lambda_table = data.frame(codon, lambda_values)

# Order the table
lambda_table = lambda_table[order(lambda_table$lambda_values),]

# Create a table of length values
length_table = data.frame(codon, length_values)

# Order the table
length_table = length_table[order(length_table$length_values),]

# Create a table of ks results
ks_table = data.frame(codon, ks_stat, p_values)

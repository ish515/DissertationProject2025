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
  
  # Import the tau data 
  file_name = paste('results/codon_', coi, '_tau_all.csv', sep='')
  tau_data = read.csv(file_name)
  
  # Output length of data
  length_values[i] <- length(tau_data[,1])
  
  # Save mean and SD of data
  mean_values[i] <- mean(tau_data[,1])
  var_values[i] <- var(tau_data[,1])
  
  # Calculate the CDF of the data
  tau_cdf = ecdf(tau_data[,1])
  
  # Pull out just cdf values in tau_cdf
  tau_cdf_only = tau_cdf(tau_data[,1])
  
  # Discretise the CDF to find the probability mass function
  pmf <- discretise(tau_cdf, method='upper', from=0, to=max(tau_data[,1]), step=0.025)
  
  # Normalise the PMF to ensure it sums to 1
  pmf <- pmf / sum(pmf)
  
  # Calculate the CDF from the PMF
  cdf_data <- cumsum(pmf)
  
  # Fit a Poisson distribution to the data
  fit <- fitdistr(tau_data[,1], "poisson")
  
  # Calculate the CDF of this Poisson distribution
  tau_seq <- seq(0, max(tau_data[,1])-0.025, by=0.025)
  cdf_poisson <- ppois(tau_seq, lambda=fit$estimate)
  
  # Save the following plot
  file_name = paste('codon_', coi, '_cdf.png', sep='')
  png(file_name)
  par(mar=c(5, 5, 4, 2))
  
  # Plot the CDF of the data and the fitted Poisson distribution
  plot(tau_seq, ppois(tau_seq, lambda=fit$estimate), type="l", col=cbPalette[6], lwd=2,
       xlab="Average Translation Time (s)", ylab="Cumulative Probability",
       cex.axis=1.5, cex.lab=1.5)
  lines(tau_seq, cdf_data, type="l", col="black", lwd=2)
  
  # Add a legend
  legend("bottomright", title='CDF', legend=c("Poisson Fit", "Data"), 
         col=c(cbPalette[6], "black"), lwd=2, bty="n", cex=1.5)
  
  dev.off()
  
  # Plot their densities
  plot(density(tau_data[,1]), col=cbPalette[2], lwd=2, xlim=c(0, 0.2), ylim=c(0, 20),
       xlab="Average Translation Time", ylab="Probability Density", main="")
  lines(density(dpois(tau_seq, lambda=fit$estimate)), col=cbPalette[6], lwd=2)
  
  # Add a legend
  legend("topright", title='Density', legend=c("Data", "Poisson Fit"), 
         col=c(cbPalette[2], cbPalette[6]), lwd=2, bty="n", cex=1.5)
  
  dev.off()
  
  # Compare the two CDFs using the Kolmogorov-Smirnov test
  ks_test <- ks.test(cdf_data, cdf_poisson)
  ks_stat[i] <- ks_test$statistic
  p_values[i] <- ks_test$p.value
  
  # Compare the observed and theoretical data using the Kolmogorov-Smirnov test
  ks_test2 <- ks.test(tau_data[,1], dpois(length(tau_data[,1]), lambda=fit$estimate))
  print(ks_test2$statistic)
  
  # Save the following plot
  file_name = paste('codon_', coi, '_histogramwithdensity.png', sep='')
  png(file_name)
  par(mar=c(5, 5, 4, 2))
  
  # Plot a histogram of the data
  hist(tau_data[,1], prob=TRUE, breaks=1000, xlim=c(0, 0.2), col="grey",
       xlab="Average Translation Time (s)", ylab="Probability Density", main="",
       cex.axis=1.5, cex.lab=1.5)
  
  # Add density estimates to the histogram
  lines(density(tau_data[,1]), col=cbPalette[2], lwd=2)
  lines(density(dpois(tau_seq, lambda=fit$estimate)), col=cbPalette[6], lwd=2)
  
  # Add a legend
  legend("topright", title='Density', legend=c("Data", "Poisson Fit"), 
         col=c(cbPalette[2], cbPalette[6]), lwd=2, bty="n",
         cex=1.5)
  
  dev.off()
  
  # Save lambda
  lambda_values[i] <- fit$estimate
  
  # Use the estimated lambda to calculate average codon dwell time
  dwell_time = fit$estimate
  
  # Append the dwell time to the dwell_times vector
  dwell_times[i] <- dwell_time
}

# Save the following plot
png("mean_vs_sd.png")
par(mar=c(5, 5, 4, 2))

# Plot mean values against SD values
plot(mean_values, var_values, xlab="Mean (s)", ylab=expression("Variance (s"^2*")"),
     xlim=c(0,0.2), ylim=c(0,0.2), cex.axis=1.5, cex.lab=1.5)

# Add the line x=y
abline(0, 1, col="red", lwd=2)

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

# Create a table of the average dwell times
dwell_times_table = data.frame(codon, dwell_times*1000)

# Create a vector containing amino acid codes
amino_acids = c("Phe", "Phe", "Leu", "Leu", "Ser", "Ser", "Ser", "Ser",
                "Tyr", "Tyr", "STOP", "STOP", "Cys", "Cys", "STOP", "Trp",
                "Leu", "Leu", "Leu", "Leu", "Pro", "Pro", "Pro", "Pro",
                "His", "His", "Gln", "Gln", "Arg", "Arg", "Arg", "Arg",
                "Ile", "Ile", "Ile", "Met", "Thr", "Thr", "Thr", "Thr",
                "Asn", "Asn", "Lys", "Lys", "Ser", "Ser", "Arg", "Arg",
                "Val", "Val", "Val", "Val", "Ala", "Ala", "Ala", "Ala",
                "Asp", "Asp", "Glu", "Glu", "Gly", "Gly", "Gly", "Gly")

# Append to the table of the average dwell times
dwell_times_table$amino_acid = amino_acids

# Order the data frame
dwell_times_table = dwell_times_table[order(dwell_times_table$dwell_times),]

# Output the table
table <- gt(dwell_times_table) %>%
  fmt_number(decimals = 3) %>%
  cols_label(codon = md("**Codon**"), 
             starts_with('dwell_times') ~ md("**Average Dwell Time (ms)**"),
             amino_acid = md("**Amino Acid**")) %>%
  gt_split(row_every_n = 32)

# Save the tables
grp_pull(table, 1) %>% gtsave(file="dwell_times_table1new.png")
grp_pull(table, 2) %>% gtsave(file="dwell_times_table2new.png")

# Save the data frame as a CSV file
write.csv(dwell_times_table, file="dwell_times_tablenew.csv", row.names=FALSE)

# Create a table of lambda values
lambda_table = data.frame(codon, lambda_values)

# Order the table
lambda_table = lambda_table[order(lambda_table$lambda_values),]

# Create a table of length values
length_table = data.frame(codon, length_values)

# Order the table
length_table = length_table[order(length_table$length_values),]

# Save the length table
length_table_print <- gt(length_table) %>%
  cols_label(codon = md("**Codon**"), 
             starts_with('length_values') ~ md("**No. of Reads**")) %>%
  gt_split(row_every_n = 16)

# Save the tables
grp_pull(length_table_print, 1) %>% gtsave(file="length_table1new.png")
grp_pull(length_table_print, 2) %>% gtsave(file="length_table2new.png")
grp_pull(length_table_print, 3) %>% gtsave(file="length_table3new.png")
grp_pull(length_table_print, 4) %>% gtsave(file="length_table4new.png")

# Create a table of ks results
ks_table = data.frame(codon, ks_stat, p_values)

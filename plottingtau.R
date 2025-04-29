# Load necessary packages
library(ggplot2)

# Import tau data for 4 genes
tau_ftsI = read.csv('tau_gene_ftsI_data.csv')
tau_secM = read.csv('tau_gene_secM_data.csv')
tau_cof = read.csv('tau_gene_cof_data.csv')
tau_yaaA = read.csv('tau_gene_yaaA_data.csv')

# Plot the data for each gene
plot_ftsI <- ggplot(data=tau_ftsI, aes(x=c(1:length(tau_ftsI[,1])), 
                                       y=tau_ftsI[,1])) +
  geom_line(color='blue') +
  xlab('Codon Number') +
  ylab('Average Translation Time') +
  ylim(0, 3) +
  scale_x_continuous(breaks=seq(0, 600, by=50), 
                     labels=seq(0, 600, by=50)) +
  theme_classic(base_size = 14)

plot_secM <- ggplot(data=tau_secM, aes(x=c(1:length(tau_secM[,1])), y=tau_secM[,1])) +
  geom_line(color='blue') +
  xlab('Codon Number') +
  ylab('Average Translation Time') +
  ylim(0, 5) +
  scale_x_continuous(breaks=seq(0, 170, by=50), 
                     labels=seq(0, 170, by=50)) +
  theme_classic(base_size = 14)

plot_cof <- ggplot(data=tau_cof, aes(x=c(1:length(tau_cof[,1])), y=tau_cof[,1])) +
  geom_line(color='blue') +
  xlab('Codon Number') +
  ylab('Average Translation Time') +
  ylim(0, 8) +
  scale_x_continuous(breaks=seq(0, 300, by=50), 
                     labels=seq(0, 300, by=50)) +
  theme_classic(base_size = 14)

plot_yaaA <- ggplot(data=tau_yaaA, aes(x=c(1:length(tau_yaaA[,1])), y=tau_yaaA[,1])) +
  geom_line(color='blue') +
  xlab('Codon Number') +
  ylab('Average Translation Time') +
  ylim(0, 6) +
  scale_x_continuous(breaks=seq(0, 250, by=50), 
                     labels=seq(0, 250, by=50)) +
  theme_classic(base_size = 14)

# Save the plots
ggsave(plot_ftsI, filename='tau_ftsI.png', device='png', width=6, height=4)
ggsave(plot_secM, filename='tau_secM.png', device='png', width=6, height=4)
ggsave(plot_cof, filename='tau_cof.png', device='png', width=6, height=4)
ggsave(plot_yaaA, filename='tau_yaaA.png', device='png', width=6, height=4)


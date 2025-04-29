# Load necessary packages
library(gt)

# Import the cleavage site data
x1 = read.csv("x1_nucleotides.csv", header=FALSE)
Y = read.csv("Y_nucleotides.csv", header=FALSE)
Z = read.csv("Z_nucleotides.csv", header=FALSE)
x2 = read.csv("x2_nucleotides.csv", header=FALSE)
before = read.csv("before_nucleotides.csv", header=FALSE)
after = read.csv("after_nucleotides.csv", header=FALSE)

# Order the data
x1 = x1[order(x1$V3, decreasing = TRUE),]
Y = Y[order(Y$V3, decreasing = TRUE),]
Z = Z[order(Z$V3, decreasing = TRUE),]
x2 = x2[order(x2$V3, decreasing = TRUE),]
before = before[order(before$V3, decreasing = TRUE),]
after = after[order(after$V3, decreasing = TRUE),]

# Edit before and after so all Ts are Us
before$V1 = gsub("T", "U", before$V1)
after$V1 = gsub("T", "U", after$V1)

# Output the table
table1 <- gt(before) %>%
  fmt_number(columns = c(V3), decimals = 3) %>%
  cols_label(V1 = md("**Upstream of Cleavage**"),
             V2 = md("**Total RNAseq Count**"),
             V3 = md("**Percentage(%)**"))
table2 <- gt(after) %>%
  fmt_number(columns = c(V3), decimals = 3) %>%
  cols_label(V1 = md("**Downstream of Cleavage**"),
             V2 = md("**Total RNAseq Count**"),
             V3 = md("**Percentage(%)**"))

# Save the table
gtsave(table1, file="cleavage_stats_table_before.png")
gtsave(table2, file="cleavage_stats_table_after.png")

# Import the expected cleavage site data
expected = read.csv("expected_genes.csv", header=FALSE)

# Import the genome-wide data
genome_pairs = read.csv("pairs_percentages.csv", header=FALSE)

# Edit so that all Ts are Us
genome_pairs$V1 = gsub("T", "U", genome_pairs$V1)

# Merge genome-wide data with cleavage site data
merged_data_before = merge(before, genome_pairs, by.x = "V1", by.y = "V1", all.x = TRUE)
merged_data_after = merge(after, genome_pairs, by.x = "V1", by.y = "V1", all.x = TRUE)

# Order the merged data
merged_data_before = merged_data_before[order(merged_data_before$V3.x, decreasing = TRUE),]
merged_data_after = merged_data_after[order(merged_data_after$V3.x, decreasing = TRUE),]

# Remove total counts columns from both
merged_data_before = merged_data_before[, -c(4)]
merged_data_after = merged_data_after[, -c(4)]

# Add a column for ratio of RNAseq percentage to genome-wide percentage
merged_data_before$V4 <- merged_data_before$V3.x / merged_data_before$V3.y
merged_data_after$V4 <- merged_data_after$V3.x / merged_data_after$V3.y

# Reorder tables by ratio
merged_data_before = merged_data_before[order(merged_data_before$V4, decreasing = TRUE),]
merged_data_after = merged_data_after[order(merged_data_after$V4, decreasing = TRUE),]

# Output the table
table3 <- gt(merged_data_before) %>%
  fmt_number(columns = c(V3.x,V3.y,V4), decimals = 3) %>%
  cols_label(V1 = md("**Upstream of Cleavage**"),
             V2.x = md("**Total RNAseq Count**"),
             V3.x = md("**Percentage(%)**"),
             V3.y = md("**Genome-wide Percentage(%)**"),
             V4 = md("**Ratio**"))
table4 <- gt(merged_data_after) %>%
  fmt_number(columns = c(V3.x,V3.y,V4), decimals = 3) %>%
  cols_label(V1 = md("**Downstream of Cleavage**"),
             V2.x = md("**Total RNAseq Count**"),
             V3.x = md("**Percentage(%)**"),
             V3.y = md("**Genome-wide Percentage(%)**"),
             V4 = md("**Ratio**"))

# Save the table
gtsave(table3, file="cleavage_stats_table_before_plustotal.png")
gtsave(table4, file="cleavage_stats_table_after_plustotal.png")


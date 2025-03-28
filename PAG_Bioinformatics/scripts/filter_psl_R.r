library(GenomicRanges)

#install.packages("dplyr")

library(dplyr)

psl_file <- "C:/Sylvain_Billiard/Masters_Thesis_Work/pollen_allergen_genes/data/genomes/psldata/blat_Pan_040.hap1_20_09.psl"
psl_data <- read.table(psl_file, header = FALSE, skip = 5)

# Define column names for PSL file
colnames(psl_data) <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert",
                        "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", "qStart", "qEnd",
                        "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts")

# additional columns ( hit_length, aln_length, avg_exon_length, percent_identity, percent_coverage )
# brief description of new columns
# hit_length = tEnd - tStart (length of the hit)
# aln_length = qEnd - qStart (length of the alignment i.e. length of the query is aligned)
# avg_exon_length = mean of blockSizes (average length of exons)
# percent_identity = matches / (matches + misMatches + repMatches + nCount) * 100 
# percent_coverage = aln_length / qSize * 100  (percentage of query covered by the alignment)

psl_data$hit_length <- psl_data$tEnd - psl_data$tStart 
psl_data$aln_length <- psl_data$qEnd - psl_data$qStart 
psl_data$sum_exon_length <- sapply(strsplit(as.character(psl_data$blockSizes), ","), function(x) sum(as.numeric(x)))
psl_data$percent_identity <- psl_data$matches / (psl_data$matches + psl_data$misMatches + psl_data$repMatches + psl_data$nCount) * 100
psl_data$percent_coverage <- psl_data$aln_length / psl_data$qSize * 100

# make a new table with only the columns we need and in specific order

psl_data_stats <- psl_data %>% select(qName, qSize, aln_length, tName, hit_length, blockCount, sum_exon_length, percent_identity, percent_coverage)

# use this table to plot the stats

install.packages("ggplot2")
library(ggplot2)

p <- ggplot(psl_data, aes(x = hit_length, y = sum_exon_length, color = as.factor(blockCount))) +
  geom_point(alpha = 0.9, size = 4) +  # Increase point size and reduce transparency
  scale_color_viridis_d(name = "Number of exons") +
  scale_x_log10() +  # Log scale for x-axis
  labs(title = "Hit Length vs.Sum Exon Length",
       x = "Hit Length (bp)",
       y = "Sum Exon Length (bp)") +
  theme_minimal(base_size = 25) +  # Increase base font size
  theme(
    axis.title = element_text(size = 30),  # Bigger axis labels
    axis.text = element_text(size = 20),   # Bigger axis ticks
    axis.ticks = element_line(size = 1.5), # Thicker tick lines
    axis.ticks.length = unit(0.4, "cm"),   # Increase tick length
    plot.title = element_text(size = 35, face = "bold"),  # Bigger title
    plot.margin = margin(30, 30, 30, 30)   # Adjust margins
  )

# Display the updated plot
p



# Plot 2: Percent Identity vs. Percent Coverage

p2 <- ggplot(psl_data, aes(x = percent_coverage, y = percent_identity, color = as.factor(blockCount))) +
  geom_point(alpha = 0.9, size = 3) +  # Increase point size and reduce transparency
  scale_color_viridis_d(name = "Number of Exons") + 
  labs(#title = "Percent Coverage vs Percent Identity",
       x = "Percent Coverage (%)",
       y = "Percent Identity (%)") +
  theme_minimal(base_size = 15) +  # Increase base font size
  theme(
    axis.title = element_text(size = 20),  # Bigger axis labels
    axis.text = element_text(size = 15),   # Bigger axis ticks
    axis.ticks = element_line(size = 1.5), # Thicker tick lines
    axis.ticks.length = unit(0.4, "cm"),   # Increase tick length
    plot.title = element_text(size = 25, face = "bold"),  # Bigger title
    plot.margin = margin(30, 30, 30, 30)   # Adjust margins
  )

# Display the updated plot
p2


# filter the psl_data_stats table to only include hits with percent_coverage > 90 and percent_identity > 75 and  260 < avg_exon_length < 340
psl_data_filtered <- psl_data %>% filter(percent_coverage > 75, percent_identity > 92.5, sum_exon_length > 400, hit_length < 2000)
psl_data_filtered <- psl_data_filtered %>% filter(grepl("^chr", tName))

psl_data_filtered
summary(psl_data_filtered)

psl_data_filtered %>% select(percent_identity,percent_coverage)

# write the filtered data to a new file
write.table(psl_data_filtered, file = "C:/Sylvain_Billiard/Masters_Thesis_Work/pollen_allergen_genes/data/genomes/psldata/blat_Pan_040_filtered_20_09.psl", sep = "\t", row.names = FALSE, quote = FALSE)


# get the entries where tName starts with "chr"




# This script takes a per_protein_conversion_ratio.csv file (ouputed from 
# ptm_ratio.py) and generates a regular barplot representing the conversion 
# ratio percentage per sample (raw file), indicating the number of PSMs employed 
# to calculate the ratio.

# Programmed by Guillermo Carrillo Martín - December 2023

library(ggplot2)
suppressPackageStartupMessages(library(dplyr))

# Load the csv files and perform some rearrangements
per_protein_csv_path <- commandArgs(trailingOnly = TRUE)[1] #Input value
per_protein_df <- read.csv(per_protein_csv_path)
per_protein_df$Conversion_ratio <- round(per_protein_df$Conversion_ratio * 100, 2) #Change ratio data to percentages

# Plot the results
per_protein_plot <- ggplot(data = per_protein_df, aes(x = Leading_razor_protein, y = Conversion_ratio)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#9CC8ED", colour = "black") +
  labs(x = "Protein", y = "Conversion ratio") +
  ylim(-5,105) +
  geom_text(aes(label = Abundance_count, x = Leading_razor_protein, fontface = "bold", y = -5)) +
  geom_text(aes(label = Conversion_ratio), vjust = -0.5) +
  facet_wrap(.~Raw_file, ncol = 3, scales = "free_x") +
  theme_bw()

path_output_folder <- dirname(per_protein_csv_path)
ggsave(file.path(path_output_folder, "per_protein_conversion_ratio.png"), plot = per_protein_plot, width = 8, height = 6)
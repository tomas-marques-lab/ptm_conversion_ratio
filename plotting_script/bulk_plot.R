# This script takes a bulk_conversion_ratio.csv file (ouputed from ptm_ratio.py)
# and generates a regular barplot representing the conversion ratio percentage 
# per sample (raw file), indicating the number of PSMs employed to calculate the
# ratio.

# Programmed by Guillermo Carrillo Martín - December 2023

library(ggplot2)

# Load the csv files and perform some rearrangements
bulk_csv_path <- commandArgs(trailingOnly = TRUE)[1]
bulk_df <- read.csv(bulk_csv_path)
bulk_df$Conversion_ratio <- round(bulk_df$Conversion_ratio * 100, 2) #Change ratio data to percentages

# Plot the results
bulk_plot <- ggplot(data = bulk_df, aes(x = Raw_file, y = Conversion_ratio)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#9CC8ED", colour = "black") +
  labs(x = "Raw file", y = "Conversion ratio") +
  ylim(-5,105) +
  geom_text(aes(label = Abundance_count, x = Raw_file, fontface = "bold", y = -5)) +
  geom_text(aes(label = Conversion_ratio), vjust = -0.5) +
  theme_bw()

path_output_folder <- dirname(bulk_csv_path)
ggsave(file.path(path_output_folder, "bulk_conversion_ratio.png"), plot = bulk_plot, width = 8, height = 6)
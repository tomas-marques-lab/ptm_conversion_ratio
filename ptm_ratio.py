"""
This script takes an evidence.txt MaxQuant file and calculates a PTM conversion ratio
between the modified / unmodified version of the amino acid. The normalization of the
process can be done using PSM or Intensity information and the ratio can be computed 
per raw file and per protein. See the "ptm_rate_calculator()" function to know more 
about the output and the way the ratio is calculated.

The script also computes a 1000 bootstrap replicates within each group (raw file or
protein). See the "ptm_bootstrap_calculator" to know more about the way this bootstrap
method works.

The script has been coded using an evidence.txt file from MaxQuant v2.0.2 and tested with
the following PTMs: [M(Oxidation (M)), P(Oxidation (P)), (Gln->pyro-Glu)Q, (Glu->pyro-Glu)E, 
N(Deamidation (NQ)), Q(Deamidation (NQ)), R(Arg->Orn), T(Phospho (ST)), S(Phospho (ST))]
"""

# Script information - Wrote in Python 3.9.12 - January 2024
__author__ = "Guillermo Carrillo Martín"
__maintainer__ = "Guillermo Carrillo Martín"
__email__ = "guillermo.carrillo@upf.edu"

import pandas as pd
import numpy as np
import argparse
import re
import os
import sys

# MAIN FUNCTION
def main():
  # Parse the input values
  evidence_path, output_folder_path, ptm_pattern, amino_acid, modification_mark, abundance_column, do_ratio_per_protein, do_remove_contaminants, do_pyro, do_bootstrap = parser()
  
  print("Computing the {} conversion ratio".format(ptm_pattern))
  if do_pyro:
      print("  PYRO SPECIAL CASE: Only the N-terminal amino acids are employed to calculate the {} ratio".format(ptm_pattern))  
  
  # Upload the evidence table, check for format errors and filter the data frame 
  evidence_df = pd.read_csv(evidence_path, sep=None, engine="python")
  format_errors(evidence_df, abundance_column, modification_mark, evidence_path)
  evidence_df = filtering(evidence_df, amino_acid, abundance_column, do_remove_contaminants)

  # Create a folder to store the results
  if not os.path.exists(output_folder_path):
    os.makedirs(output_folder_path)

  # Calculate bulk (and per protein) conversion ratio
  if do_bootstrap:
    bulk_conversion_ratio_df, bulk_bootstrap_df= ptm_rate_calculator(evidence_df, amino_acid, modification_mark, ptm_pattern, abundance_column, do_pyro, do_bootstrap)
    bulk_conversion_ratio_df.to_csv("{}/bulk_conversion_ratio.csv".format(output_folder_path))
    bulk_bootstrap_df.to_csv("{}/bulk_bootstrap_replicates.csv".format(output_folder_path), index=False)
    if do_ratio_per_protein:
      per_proten_conversion_ratio_df, per_protein_bootstrap_df = ptm_rate_calculator(evidence_df, amino_acid, modification_mark, ptm_pattern, abundance_column, do_pyro, do_bootstrap, do_per_protein=True)
      per_proten_conversion_ratio_df.to_csv("{}/per_protein_conversion_ratio.csv".format(output_folder_path))
      per_protein_bootstrap_df.to_csv("{}/per_protein_bootstrap_replicates.csv".format(output_folder_path), index=False)
  
  elif not do_bootstrap:
    bulk_conversion_ratio_df = ptm_rate_calculator(evidence_df, amino_acid, modification_mark, ptm_pattern, abundance_column, do_pyro, do_bootstrap)
    bulk_conversion_ratio_df.to_csv("{}/bulk_conversion_ratio.csv".format(output_folder_path))
    if do_ratio_per_protein:
      per_proten_conversion_ratio_df = ptm_rate_calculator(evidence_df, amino_acid, modification_mark, ptm_pattern, abundance_column, do_pyro, do_bootstrap, do_per_protein=True)
      per_proten_conversion_ratio_df.to_csv("{}/per_protein_bootstrap_replicates.csv".format(output_folder_path))

# FUNCTIONS
def parser():
  """
  This function parses and inputs the required arguments from the terminal to the python script
  *OUTPUT
  - evidence_path; A string indicating the path to the evidence.txt file
  - output_folder_path; A string indicating the path to locate the folder with the result 
    (default: evidence.txt containing folder)
  - ptm_pattern; A string indicating the amino acid + PTM pattern, indicated in the way it's wrote 
    in the 'Modified sequence' column e.g. "Q(Deamidation (NQ))" or "(Glu->pyro-Glu)E"
  - amino_acid; A string indicating the unmodified aminoacid e.g. "Q" or "E"
  - modification_mark; A string indicating the PTM modification mark e.g. "(Deamidation (NQ))" 
    or "(Glu->pyro-Glu)"
  - abundance_column; A string indicating the name of the column to take the abundance values
    from the evidence.txt file. Usually "MS/MS count" or "Intensity".
  - do_ratio_per_protein; A boolean indicating if the conversion ratio per protein is performed  
  - do_remove_contaminants; A boolean indicating if the protein contaminants are removed before 
    the calculation
  - do_pyro; A boolean indicating if the PTM pattern belongs to a 'pyro speial case'. 
    See the 'ptm_rate_calculator()' function to understands the implications of calculating
    the PTM ratios under the 'pyro special case'
  - do_bootstrap; A boolean indicating if a 1000 bootstrap replicate is perform within each calculation
  """
  parser = argparse.ArgumentParser(description="This script takes an evidence.txt MaxQuant file and calculates a PTM conversion ratio between the unmodified / modified version of the amino acid")
  parser.add_argument("--evidence-file", "-e", dest="evidence_file", type=str, help="The path to the MaxQuant evidence.txt", required=True, nargs=1)
  parser.add_argument("--ptm", "-p", dest="ptm", type=str, help="The amino acid + PTM pattern, indicated in the way it's wrote in the 'Modified sequence' column e.g. Q(Deamidation (NQ)) or (Glu->pyro-Glu)E", required=True, nargs=1)
  parser.add_argument("--output-path", "-o", dest="path", type=str, help="The path to locate the folder with the result (default: evidence.txt containing folder)", required=False, nargs=1)
  parser.add_argument("--abundance-column", "-a", dest="abundance_column", type=str, help=" (default: 'MS/MS count')", default=["MS/MS count"], required=False, nargs=1)
  parser.add_argument("--per-protein", dest="per_protein", action=argparse.BooleanOptionalAction, help="Indicate if the ratio is also computed per protein", default=True, required=False)
  parser.add_argument("--remove-contaminants", dest="remove_contaminants", action=argparse.BooleanOptionalAction, help="Indicate if the contaminant and reverse proteins are removed to calculate the ratio", default=True, required=False)
  parser.add_argument("--bootstrap", dest="bootstrap", action=argparse.BooleanOptionalAction, help="Indicate if a 1000 bootstrap replicate is perform within each calculation", default=True, required=False)

  # Recovering the arguments 
  args = parser.parse_args()
  evidence_path = args.evidence_file[0]
  ptm_pattern = args.ptm[0]
  amino_acid, modification_mark = ptm_mark_split_aa_and_modification(ptm_pattern)

  # Use the PTM + the name of the evidence file to create the name of the folder
  output_folder_name = "{}_{}".format(ptm_pattern.replace(" ", "_"), os.path.splitext(os.path.basename(evidence_path))[0])
  
  # If the ouput path is not indicated, the results would be output in the same folder that contains the evidence.txt file
  if args.path:
    output_folder_path = args.path[0] + "/" + output_folder_name
  else:
    output_folder_path = os.path.dirname(evidence_path) + "/" + output_folder_name

  abundance_column = args.abundance_column[0]

  do_ratio_per_protein = args.per_protein
  do_remove_contaminants = args.remove_contaminants
  do_bootstrap = args.bootstrap
  if ptm_pattern.startswith("("):
    do_pyro = True
  elif not ptm_pattern.startswith("("):
    do_pyro = False

  return evidence_path, output_folder_path, ptm_pattern, amino_acid, modification_mark, abundance_column, do_ratio_per_protein, do_remove_contaminants, do_pyro, do_bootstrap

def ptm_mark_split_aa_and_modification(ptm_pattern):
    """
    This function takes a PTM pattern and retrieves the amino acid and modification mark
    *INPUT
    - ptm_pattern; A string indicating the amino acid + PTM pattern, indicated in the way it's wrote 
      in the 'Modified sequence' column e.g. "Q(Deamidation (NQ))" or "(Glu->pyro-Glu)E"
    *OUTPUT
    - amino_acid; A string indicating the unmodified aminoacid e.g. "Q" or "E"
    - modification_mark; A string indicating the PTM modification mark e.g. "(Deamidation (NQ))" 
      or "(Glu->pyro-Glu)"
    """
    if ptm_pattern.startswith("("):
      amino_acid = ptm_pattern[-1]
      modification_mark = ptm_pattern[:-1] 

    elif not ptm_pattern.startswith("("):
      amino_acid = ptm_pattern[0]
      modification_mark = ptm_pattern[1:] 

    return amino_acid, modification_mark

def format_errors(evidence_df, abundance_column, modification_mark, evidence_path):
  """
  This function stops the code and prints an error if some conditions are fulfill
  *INPUT
  - evidence_df; The evidence dataframe inputed from the MaxQuant evidence.txt file
  - abundance_column; A string indicating the name of the column to take the abundance values
    from the evidence.txt file. Usually "MS/MS count" or "Intensity"
  - modification_mark; A string indicating the PTM modification mark e.g. "(Deamidation (NQ))" 
    or "(Glu->pyro-Glu)"
  """
  # The abundance column does not exist in the evidence.txt fil
  try:
    if abundance_column not in evidence_df.columns:
      raise ValueError("The '{}' column does not exist in this evidence.txt file".format(abundance_column))
  except ValueError as ve:
    print("ERROR:", ve)
    sys.exit(1)

  # The way the PTM was inputed is not valid 
  try:
    if modification_mark.count("(") != modification_mark.count(")"):
      raise ValueError("The '{}' format is not valid. e.g. 'M(Oxidation (M)) or '(Glu->pyro-Glu)E'".format(modification_mark))
  except ValueError as ve:
    print("ERROR:", ve)
    sys.exit(1)

  # The evidence.txt file was not indicated as a path
  try:
    if not evidence_path.startswith((".", "/")):
      raise ValueError("The '{}' path format is not valid. e.g. './evidence.txt' or /Users/user1/Documents/evidence.txt".format(evidence_path))
  except ValueError as ve:
    print("ERROR:", ve)
    sys.exit(1)
  
  
def filtering(evidence_df, amino_acid, abundance_column, do_remove_contaminants):
  """
  This function filters the dataset by removing Reverse and Contaminant peptides (if required) and
  eliminating the rows without an abundance value (like rows without Intensity due to MS/MS errors).
  It also raises an error if the targeted amino acid is not found in the filtered dataset.
  *Input
  - evidence_df; The evidence dataframe inputed from the MaxQuant evidence.txt file
  - amino_acid; A string indicating the unmodified aminoacid e.g. "Q" or "E"
  - abundance_column; A string indicating the name of the column to take the abundance values
    from the evidence.txt file. Usually "MS/MS count" or "Intensity"
  - do_remove_contaminants; A boolean indicating if the protein contaminants are removed before 
    the calculation
  *Output
  - evidence_df; The filtered version of the dataframe
  """
  if do_remove_contaminants:
    evidence_df = evidence_df[~evidence_df["Leading razor protein"].str.startswith("REV_")]
    evidence_df = evidence_df[~evidence_df["Leading razor protein"].str.startswith("CON_")]
  
  # Remove all the rows without an abundance value
  evidence_df.dropna(subset=[abundance_column], inplace=True)

  # Check if the targeted amino acid is present in the database sequences
  try:
    if not evidence_df["Sequence"].str.contains(amino_acid).any():
      raise ValueError("Dataset unfeasible to calculate the conversion ratio. The '{}' amino acid is not present in the dataset\n".format(amino_acid))
  except ValueError as ve:
    print("ERROR:", ve)
    sys.exit(1)
  
  return evidence_df

def ptm_rate_calculator(evidence_df, amino_acid, modification_mark, ptm_pattern, abundance_column, do_pyro, do_bootstrap, do_per_protein=False):
  """
  This function takes a evidence.txt MaxQuant file (inputed as a pandas dataframe) and calculates a 
  PTM conversion ratio per raw file and protein. The way to normalize the data can be PSM or Intensity
  based, inspired by the David Lyon deamidation script (https://github.com/dblyon/deamidation)

  To calculate the conversion ratio, the total number of unmodified and modified amino acids are count 
  for each peptide. Then, this values are divided (modified / total) and multiplied by the "peptide 
  spectrum matches" (PSM) or the "Intensity" abundance values of the peptide. Afterwards, the data is 
  splitted by raw file (and protein) and the normalized ratio of the peptides is summed up and divided by 
  the total abundance value within each group.

  *PYRO SPECIAL CASE: If the PTM pattern looks like (Glu->pyro-Glu)E or (Gln->pyro-Glu)Q, the calculations
  are only based on the first amino acid of each peptide. This decition is made because the pyro-Glu conversion
  only occurs at the N-terminal amino acid. 

  *INPUT:
  - evidence_df; The evidence dataframe inputed from the MaxQuant evidence.txt file
  - amino_acid; A string indicating the unmodified aminoacid e.g. "Q" or "E" 
  - modification_mark; A string indicating the PTM modification mark e.g. "(Deamidation (NQ))" 
    or "(Glu->pyro-Glu)"
  - ptm_pattern; A string indicating the amino acid + PTM pattern, indicated in the way it's wrote 
    in the 'Modified sequence' column e.g. "Q(Deamidation (NQ))" or "(Glu->pyro-Glu)E"
  - abundance_column; A string indicating the name of the column to take the abundance values
    from the evidence.txt file. Usually "MS/MS count" or "Intensity"
  - do_pyro; A boolean indicating if the PTM pattern belongs to a 'pyro speial case'
  - do_bootstrap; A boolean indicating if a 1000 bootstrap replicate is perform within each calculation
  - do_per_protein; A boolean indicating if the conversion ratio per protein is performed
  *OUTPUT:
  - conversion_rate_df; A dataframe containing the following calculation results:
    a. Amino_acid_count: Number of amino acids employed to calculate the conversion ratio
    b. Conversion_count: Number of modified amino acids observed in the dataset
    c. Abundance_count: Number of "Peptide Spectrum Matches" employed to calculate the 
      conversion ratio. NOTE; the number of unmodified amino acids can can be higher than 
      the number of PSMs because a PSM (peptide) can contain the desired amino acid more than once.
    d. Conversion_ratio: Conversion ratio normalized by the abundance values
    *If bootstrapping:
    e'. Boostrap_mean_converion_ratio; The conversion ratio mean obtained from the bootstrap replicates
    f'. Q95_low; The value that delimits the bottom of the 95% bootstrap replicates distribution
    g'. Q95_up; The value that delimits the top of the 95% bootstrap replicates distribution
      +The Q95_low and Q95_up values delimit the 95% of the values from the bootstrap replicates distribution
  - bootstrap_replicates_df; A dataframe containing all the conversion ratio results per each
    bootstrap replicate, numerating the replicates from 0 to 999 per Raw file (and protein).
  """

  # Change some column names to make the format more parseable
  evidence_df.rename(columns={"Raw file": "Raw_file", 
                              "Leading razor protein": "Leading_razor_protein"},
                              inplace=True)

  # Remove the sequences without the targeted aminoacid and calculate the conversion ratio per peptide
  if do_pyro:
    evidence_df = evidence_df[evidence_df["Sequence"].str[0] == amino_acid].copy()
    evidence_df["Amino acid count"] = np.where(evidence_df["Sequence"].str[0] == amino_acid, 1, 0)
  elif not do_pyro:
    evidence_df = evidence_df[evidence_df["Sequence"].str.contains(amino_acid)].copy()
    evidence_df["Amino acid count"] = evidence_df["Sequence"].str.count(amino_acid)

  evidence_df["Conversion count"] = evidence_df["Modified sequence"].str.count(re.escape(ptm_pattern)) #re.escape avoids the interpretation of the string as a RegEx
  evidence_df["Conversion ration norm"] = (evidence_df["Conversion count"] / evidence_df["Amino acid count"]) * evidence_df[abundance_column]
  
  # Create the output dtaframe (conversion_rate_df) and calculate the bulk conversion ratio per RawFile (and protein)
  conversion_rate_df = pd.DataFrame()

  groupby_column = ["Raw_file"]
  if do_per_protein:
    groupby_column.append("Leading_razor_protein")

  conversion_rate_df["Amino_acid_count"] = evidence_df.groupby(groupby_column)["Amino acid count"].apply(lambda x: x.sum())
  conversion_rate_df["Conversion_count"] = evidence_df.groupby(groupby_column)["Conversion count"].apply(lambda x: x.sum())
  conversion_rate_df["Abundance_count"] = evidence_df.groupby(groupby_column)[abundance_column].apply(lambda x: x.sum())
  conversion_rate_df["Conversion_ratio"] = evidence_df.groupby(groupby_column)["Conversion ration norm"].apply(lambda x: x.sum()) / conversion_rate_df["Abundance_count"]
  
  conversion_rate_df["Amino_acid"] = amino_acid
  conversion_rate_df["PTM"] = modification_mark[1:-1]

  # Order the dataframe so the PTM column goes before the count columns
  column_list = conversion_rate_df.columns.tolist()
  column_list.insert(column_list.index("Amino_acid_count"), column_list.pop(column_list.index("Amino_acid")))
  column_list.insert(column_list.index("Amino_acid_count"), column_list.pop(column_list.index("PTM")))
  conversion_rate_df = conversion_rate_df[column_list]

  # Bootstrap
  if do_bootstrap:
    bootstrap_replicates_df, bootstrap_statistics_df = ptm_bootstrap_calculator(evidence_df, amino_acid, modification_mark, abundance_column, do_per_protein)
    conversion_rate_df = conversion_rate_df.merge(bootstrap_statistics_df, on=groupby_column)
    return conversion_rate_df, bootstrap_replicates_df

  elif not do_bootstrap:
    return conversion_rate_df

def ptm_bootstrap_calculator(conversion_ratio_evidence_df, amino_acid, modification_mark, abundance_column, do_per_protein):
  """
  This function takes as an input a evidence.txt MaxQuant file (inputed as a pandas dataframe) 
  with a column indicating the PTM conversion ratio per row, computed as part of the 
  "ptm_rate_calculator()" function. Then, it randomly samples the dataset to generate 1000 replicates
  of the dataset. In each replicate, the peptides can be over or underrepresented compared to the original
  dataset. Then, for each replicate, the conversion ratio is computed as in the "ptm_rate_calculator()" 
  function. The mean conversion ratio for all the bootstrap replicates and the 95% quartil up/down 
  delimiter values are computed.

  *INPUT:
  - conversion_ratio_evidence_df; The evidence dataframe with the conversion ratio calculated per row
  - amino_acid; A string indicating the unmodified aminoacid e.g. "Q" or "E" 
  - modification_mark; A string indicating the PTM modification mark e.g. "(Deamidation (NQ))" or 
    "(Glu->pyro-Glu)"
  - abundance_column; A string indicating the name of the column to take the abundance values
    from the evidence.txt file. Usually "MS/MS count" or "Intensity"
  - do_per_protein; A boolean indicating if the conversion ratio per protein is performed
  *OUTPUT:
  - bootstrap_replicates_df; A dataframe containing all the conversion ratio results per each
    bootstrap replicate, numerating the replicates from 0 to 999 per Raw file (and protein).
  - bootstrap_statistics_df; A dataframe including the mean conversion ratio and 95% quartil
    delimiter values for all the bootstrap replicates per Raw file (and protein).
  """
  groupby_column = ["Raw_file"]
  if do_per_protein:
    groupby_column.append("Leading_razor_protein")

  col_names = groupby_column + ["Replicate", "Amino_acid", "PTM", "Conversion_ratio"]
  bootstrap_replicates_df = pd.DataFrame(columns=col_names)

  # Bootstraping
  for tag, group_df in conversion_ratio_evidence_df.groupby(groupby_column):
    for replicate in range(1000):
      # Resample the data with replacement
      random_sample_evidence_df = group_df.sample(n=len(group_df), replace=True)  
      
      # Compute the calculation
      abundance_sum = random_sample_evidence_df[abundance_column].sum()
      conversion_sum = random_sample_evidence_df["Conversion ration norm"].sum()
      ptm_conversion_ratio = conversion_sum / abundance_sum

      # Writte the result of each iteration in a dataframe
      if not do_per_protein:
        new_row = [tag] + [replicate, amino_acid, modification_mark[1:-1], ptm_conversion_ratio]
      elif do_per_protein:
        new_row = list(tag) + [replicate, amino_acid, modification_mark[1:-1], ptm_conversion_ratio]
      
      new_row_serie = pd.Series(new_row, index=bootstrap_replicates_df.columns)
      bootstrap_replicates_df = pd.concat([bootstrap_replicates_df, pd.DataFrame([new_row_serie], columns=bootstrap_replicates_df.columns)])

  # Create and write the a df including some bootstrap statistics 
  bootstrap_statistics_df = pd.DataFrame()
  bootstrap_statistics_df["Boostrap_mean_converion_ratio"] = bootstrap_replicates_df.groupby(groupby_column)["Conversion_ratio"].mean()
  bootstrap_statistics_df["Q95_low"] = bootstrap_replicates_df.groupby(groupby_column)["Conversion_ratio"].apply(lambda x: x.quantile(0.025))
  bootstrap_statistics_df["Q95_up"] = bootstrap_replicates_df.groupby(groupby_column)["Conversion_ratio"].apply(lambda x: x.quantile(0.925))

  return bootstrap_replicates_df, bootstrap_statistics_df

if __name__ == "__main__":
  main()
"""
This script takes an evidence.txt MaxQuant file and calculates a PTM conversion ratio
between the modified / unmodified version of the amino acid. The normalization of the
process is performed using PSM information and the ratio can be computed per raw file and 
per protein. 

See the "ptm_rate_calculator()" function to know more about the output and the way 
the ratio is calculated.

The script has been coded using an evidence.txt file from MaxQuant v2.0.2 and tested with
the following PTMs: [M(Oxidation (M)), P(Oxidation (P)), (Gln->pyro-Glu)Q, (Glu->pyro-Glu)E, 
N(Deamidation (NQ)), Q(Deamidation (NQ)), R(Arg->Orn), T(Phospho (ST)), S(Phospho (ST))]
"""

# Script information - Wrote in Python 3.9.12 - December 2023
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
    evidence_path, output_folder_path, ptm_pattern, amino_acid, modification_mark, do_ratio_per_protein, do_remove_contaminants, do_pyro = parser()
    print("Computing the {} conversion ratio".format(ptm_pattern))
    
    # Upload the evidence table, remove contaminants and check if there are sequences with the desire amino acid
    evidence_df = pd.read_csv(evidence_path, sep=None, engine="python")
    if do_remove_contaminants:
      evidence_df = evidence_df[~evidence_df["Leading razor protein"].str.startswith("REV_")]
      evidence_df = evidence_df[~evidence_df["Leading razor protein"].str.startswith("CON_")]
    amino_acid_presence_check(evidence_df, amino_acid)

    # Create a folder to store the results
    if not os.path.exists(output_folder_path):
      os.makedirs(output_folder_path)

    # Calculate bulk (and per protein) conversion ratio
    if do_pyro:
       print("  PYRO SPECIAL CASE: Only the N-terminal amino acids are employed to calculate the {} ratio".format(ptm_pattern))  
    
    bulk_conversion_ratio_df = ptm_rate_calculator(evidence_df, amino_acid, modification_mark, ptm_pattern, do_pyro)
    bulk_conversion_ratio_df.to_csv("{}/bulk_conversion_ratio.csv".format(output_folder_path))
    
    if do_ratio_per_protein:
      per_proten_conversion_ratio_df = ptm_rate_calculator(evidence_df, amino_acid, modification_mark, ptm_pattern, do_pyro, do_per_protein=True)
      per_proten_conversion_ratio_df.to_csv("{}/per_protein_conversion_ratio.csv".format(output_folder_path))

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
    - do_ratio_per_protein; A boolean indicating if the conversion ratio per protein is performed  
    - do_remove_contaminants; A boolean indicating if the protein contaminants are removed before 
      the calculation
    - do_pyro; A boolean indicating if the PTM pattern belongs to a 'pyro speial case'. 
      See the 'ptm_rate_calculator()' function to understands the implications of calculating
      the PTM ratios under the 'pyro special case'
    """
    parser = argparse.ArgumentParser(description="This script takes an evidence.txt MaxQuant file and calculates a PTM conversion ratio between the unmodified / modified version of the amino acid")
    parser.add_argument("--evidence-file", "-e", dest="evidence_file", type=str, help="The path to the MaxQuant evidence.txt", required=True, nargs=1)
    parser.add_argument("--ptm", "-p", dest="ptm", type=str, help="The amino acid + PTM pattern, indicated in the way it's wrote in the 'Modified sequence' column e.g. Q(Deamidation (NQ)) or (Glu->pyro-Glu)E", required=True, nargs=1)
    parser.add_argument("--output-path", "-o", dest="path", type=str, help="The path to locate the folder with the result (default: evidence.txt containing folder)", required=False, nargs=1)
    parser.add_argument("--per-protein", dest="per_protein", action=argparse.BooleanOptionalAction, help="Indicate if the ratio is also computed per protein", default=True, required=False)
    parser.add_argument("--remove-contaminants", dest="remove_contaminants", action=argparse.BooleanOptionalAction, help="Indicate if the contaminant and reverse proteins are removed to calculate the ratio", default=True, required=False)

    #Recovering the arguments 
    args = parser.parse_args()
    evidence_path = args.evidence_file[0]
    ptm_pattern = args.ptm[0]
    amino_acid, modification_mark = ptm_mark_split_aa_and_modification(ptm_pattern)

    ##Use the PTM + the name of the evidence file to create the name of the folder 
    output_folder_name = "{}_{}".format(ptm_pattern.replace(" ", "_"), os.path.splitext(os.path.basename(evidence_path))[0])
    
    ##If the ouput path is not indicated, the results would be output in the same folther that contains the evidence.txt file
    if args.path:
      output_folder_path = args.path[0] + "/" + output_folder_name
    else:
      output_folder_path = os.path.dirname(evidence_path) + "/" + output_folder_name

    do_ratio_per_protein = args.per_protein
    do_remove_contaminants = args.remove_contaminants
    if ptm_pattern.startswith("("):
      do_pyro = True
    elif not ptm_pattern.startswith("("):
      do_pyro = False

    #Possible format error
    if modification_mark.count("(") != modification_mark.count(")"):
      print("ERROR: Check PTM format; e.g. 'M(Oxidation (M))'\n")
      sys.exit(1)

    return evidence_path, output_folder_path, ptm_pattern, amino_acid, modification_mark, do_ratio_per_protein, do_remove_contaminants, do_pyro

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

def amino_acid_presence_check(evidence_df, amino_acid):
    """
    This function takes the evidence dataframe and searchs if the desired amino acid is
    present in the dataset. If not, the code stops and an Error message is printed.
    *INPUT
    - evidence_df; The evidence dataframe inputed from the MaxQuant evidence.txt file
    - amino_acid; A string indicating the unmodified aminoacid
    *OUTPUT
    - An error message if the condition is fulfill
    """
    if not evidence_df["Sequence"].str.contains(amino_acid).any():
      print("ERROR: Dataset unfeasible to calculate the conversion ratio. The '{}' amino acid is not present in the dataset\n".format(amino_acid))
      sys.exit(1)

def ptm_rate_calculator(evidence_df, amino_acid, modification_mark, ptm_pattern, do_pyro, do_per_protein=False):
    """
    This function takes a evidence.txt MaxQuant file (inputed as a panda dataframe) and calculates a 
    PTM conversion ratio per raw file and protein. The way to normalize the data is PSM-based, inspired 
    by the David Lyon deamidation script (https://github.com/dblyon/deamidation)

    To calculate the conversion ratio, the total number of unmodified and modified amino acids are count 
    for each peptide. Then, this values are divided (modified / total) and multiplied by the "peptide 
    spectrum matches" (PSM) of the peptide. Afterwards, the data is slitted by raw file (and protein) and 
    the normalized ratio of the peptides is summed up and divided by the total PSM within each group.

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
    - do_pyro; A boolean indicating if the PTM pattern belongs to a 'pyro speial case'
    - do_per_protein; A boolean indicating if the conversion ratio per protein is performed
    *OUTPUT:
    - A csv file containing the calculation results:
        a. Amino acid count: Number of amino acids employed to calculate the conversion ratio
        b. Conversion count: Number of modified amino acids observed in the dataset
        c. Total PSM count: Number of "Peptide Spectrum Matches" employed to calculate the 
           conversion ratio. NOTE; the number of unmodified amino acids can can be higher than 
           the number of PSMs because a PSM (peptide) can contain the desired amino acid more than once.
        d. PTM norm ratio: Conversion ratio normalized by the PSM counts
    """

    #Change some column names to make the format more parseable
    evidence_df.rename(columns={"Raw file": "Raw_file", 
                               "Leading razor protein": "Leading_razor_protein"},
                               inplace=True)

    #Remove the sequences without the aminoacid and calculate the conversion ratio per peptide
    if do_pyro:
      evidence_df = evidence_df[evidence_df["Sequence"].str[0] == amino_acid].copy()
      evidence_df["Amino acid count"] = np.where(evidence_df["Sequence"].str[0] == amino_acid, 1, 0)
    elif not do_pyro:
      evidence_df = evidence_df[evidence_df["Sequence"].str.contains(amino_acid)].copy()
      evidence_df["Amino acid count"] = evidence_df["Sequence"].str.count(amino_acid)

    evidence_df["Conversion count"] = evidence_df["Modified sequence"].str.count(re.escape(ptm_pattern)) #re.escape avoids the interpretation of the string as a RegEx
    evidence_df["Conversion ration norm"] = (evidence_df["Conversion count"] / evidence_df["Amino acid count"]) * evidence_df["MS/MS count"]
    
    #Create the conversion dataframe and calculate the bulk conversion ratio per RawFile (and protein)
    conversion_rate_df = pd.DataFrame()

    groupby_column = ["Raw_file"]
    if do_per_protein:
      groupby_column.append("Leading_razor_protein")

    conversion_rate_df["Amino_acid_count"] = evidence_df.groupby(groupby_column)["Amino acid count"].apply(lambda x: x.sum())
    conversion_rate_df["Conversion_count"] = evidence_df.groupby(groupby_column)["Conversion count"].apply(lambda x: x.sum())
    conversion_rate_df["Total_PSM_count"] = evidence_df.groupby(groupby_column)["MS/MS count"].apply(lambda x: x.sum())
    conversion_rate_df["PTM_norm_ratio"] = evidence_df.groupby(groupby_column)["Conversion ration norm"].apply(lambda x: x.sum()) / conversion_rate_df["Total_PSM_count"]
    conversion_rate_df["Amino_acid"] = amino_acid
    conversion_rate_df["PTM"] = modification_mark[1:-1]

    #Order the dataframe so the PTM column goes before the count columns
    column_list = conversion_rate_df.columns.tolist()
    column_list.insert(column_list.index("Amino_acid_count"), column_list.pop(column_list.index("Amino_acid")))
    column_list.insert(column_list.index("Amino_acid_count"), column_list.pop(column_list.index("PTM")))
    conversion_rate_df = conversion_rate_df[column_list]

    return conversion_rate_df

if __name__ == "__main__":
  main()
#!/usr/local/bin/bash

## DISCLAIMER; the script folder needs to be in the same directory as the .sh script
## DISCLAIMER; 

# Input
ptm_list_path=/Users/guillermocarrillomartin/IBE_TMB/7_software/ptm_conversion_ratio/ptm.txt
evidence_path=/Users/guillermocarrillomartin/IBE_TMB/neandertal_project/analysis/merged_evidence.txt
directory_script_path=/Users/guillermocarrillomartin/IBE_TMB/7_software/ptm_conversion_ratio
output_directory=$(dirname "$evidence_path")

# Script
## Calculate the ratios using the ptm_ratio.py script
echo "#CALCULATING THE CONVERSION RATIO"

cd $directory_script_path
sed 's/UPI[A-Z0-9]\{10\}_//g' $evidence_path > ./ptm_conversion

mapfile -t ptm_array < $ptm_list_path
for ptm in "${ptm_array[@]}"
do 
    python3 script/ptm_ratio.py -e ./ptm_conversion -p "$ptm" -o $output_directory -a "Intensity"
done
rm ptm_conversion

## Plot the results
echo "#PLOTTING THE RESULTS"
directory_path=$(dirname "$evidence_path")
find $directory_path -type f -name "bulk_conversion_ratio.csv" | while IFS= read -r file 
do
    Rscript $directory_script_path/script/bulk_plot.R $file #bulk conversion ratio
done

find $directory_path -type f -name "per_protein_conversion_ratio.csv" | while IFS= read -r file 
do
    Rscript $directory_script_path/script/per_protein_plot.R $file #per protein conversion ratio
done

echo "#DONE"

# PTM ratio calculator

The **ptm_ratio.py** script takes an evidence.txt MaxQuant file and calculates a PTM conversion ratio between the modified / unmodified version of the amino acid. The normalization of the process is performed using PSM information and the ratio can be computed per raw file and per protein.

**Disclaimer**: The script has been tested only using an evidence.txt file from *MaxQuant v2.0.2*

## How do I use it?

Just type `python3 ptm_ratio.py -h` or `python3 ptm_ratio.py --help` to see the execution parameters. In detail, the option parameters are:

**Mandatory**

+----------------------+--------------------------------------------------------------------------------------------------+
| \--evidence-file; -e | The path to the MaxQuant evidence.txt file                                                       |
+----------------------+--------------------------------------------------------------------------------------------------+
| \--ptm; -p           | The amino acid + PTM pattern, indicated in the way it's wrote in the 'Modified sequence' column. |
|                      |                                                                                                  |
|                      | e.g. "Q(Deamidation (NQ))" or "(Glu-\>pyro-Glu)E"                                                |
+----------------------+--------------------------------------------------------------------------------------------------+

**Optional**

+---------------------------------------------------+-----------------------------------------------------------------------------------------------------+
| ---output-path; -o                                | The path to locate the folder with the results (default: evidence.txt containing folder)            |
+---------------------------------------------------+-----------------------------------------------------------------------------------------------------+
| ---per-protein, ---no-per-protein                 | Indicate if the ratio is also computed per protein (default: True)                                  |
+---------------------------------------------------+-----------------------------------------------------------------------------------------------------+
| ---remove-contaminants, ---no-remove-contaminants | Indicate if the contaminant and reverse proteins are removed to calculate the ratio (default: True) |
+---------------------------------------------------+-----------------------------------------------------------------------------------------------------+

The following code line exemplifies the way to execute the script to calculate the Q deamidation ratio:

```{python}
python3 ptm_ratio.py --evidence-file /example/path/evidence.txt --ptm "Q(Deamidation (NQ))" --no-remove-contaminants 
```

## How does it works?

To calculate the conversion ratio, the total number of unmodified and modified amino acids are count for each peptide. Then, this values are divided (modified / total) and multiplied by the "peptide spectrum matches" (PSM) of the peptide. Afterwards, the data is slitted by raw file (and protein) and the normalized ratio of the peptides is summed up and divided by the total PSM within each group. A more detailed explanation (including an example) can be seen in [David Lyon's deamidation script](https://github.com/dblyon/deamidation)

If the PTM pattern is "(Glu-\>pyro-Glu)E" or "(Gln-\>pyro-Glu)Q", the calculations are only based on the first amino acid of each peptide. This decision is made because the pyro-Glu conversion only occurs at the N-terminal amino acid.

## Testing and examples

An three-PTM example driven by two 10Ma *Deinotherium* samples can be found in the [**Example folder**](example)

The script has been coded using an evidence.txt file from MaxQuant v2.0.2 and tested with the following PTMs: [M(Oxidation (M)), P(Oxidation (P)), (Gln-\>pyro-Glu)Q, (Glu-\>pyro-Glu)E, N(Deamidation (NQ)), Q(Deamidation (NQ)), R(Arg-\>Orn), T(Phospho (ST)), S(Phospho (ST))]

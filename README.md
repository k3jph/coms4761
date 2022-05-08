
# Computational Genomics Final Project

This repository contains a final project for COMS 4761 -- Computational
Genomics at Columbia University in the City of New York.

* GitHub repo: <https://github.com/k3jph/coms4761.git>
* Free software: MIT

Our project replicates aspects of the
[GRADIS](https://github.com/MonaRazaghi/GRADIS) for the supervised
learning of gene regulatory networks based on graph distance profile
of transcriptomics data.  For this, we have replicated the basic
analysis in R, and used methods other than SVM.  In addition, to
simplify some of the processing, we have exported the data from
Excel to CSVs.

## Features

All referenced files can be found in `src`, which is organized:

* `data` folder containing *E. coli* and *S. cerevisiae* inputs in the correct format for GRADIS use. Also includes the raw DREAM5 data for both organisms in `raw` and the pre-processing script used to transform to the proper input form/files.
* `gnn` folder containing the GRGNN implementation from <https://github.com/juexinwang/GRGNN>, as well as the output results used in our analysis. Includes a separate README file with quickstart instructions.
* `R` folder contains the majority of our GRADIS implementation scripts, as translated from the original MATLAB implementation here: <https://github.com/MonaRazaghi/GRADIS>.

## Quickstart

**Data Transformation**

The GRADIS scripts expect three input files: 
* *Genes.csv* with the list of gene names for the organism
* *Network.csv* listing labeled TF-gene interactions (where 1 indicates a positive interaction)
* *Expression.txt* with collected gene expression levels. Each row represents a sample. 

Raw DREAM5 data for *E. coli* (Network 3) and *S. cerevisiae* (Network 4) is available in `data\raw`, as pulled from <https://www.synapse.org/#!Synapse:syn2787211>. To transform these files into the expected format for our GRADIS implementation (discussed below), navigate to the `data` directory, and run the following command:
```
python data_processing.py <network_id>
```
where network_id = 3 for *E. coli* and network_id = 4 for *S. cerevisiae*.

## Usage

### GNN

### R

To test three different methods for differentiating 
between positive and negative, for *E. coli* use:

    Rscript src/R/glmmer-ec.R

For *S. cerevisiae*, use:

    Rscript src/R/glmmer-sc.R

Both scripts include a vanilla GLM, random forest, and 
naïve Bayes as classifiers.  The code is built on top 
of the R library [caret](https://topepo.github.io/caret/),
so switching classifiers is trivial.

To identify and train negatives over the *E. coli* data,
use:

    Rscript src/R/gradis-neg-ec.R

To identify and train negatives over the *S. cerevisiae* data,
use:

    Rscript src/R/gradis-neg-sc.R

Both scripts use random forest by default, but also
use caret, so switching is, again trivial.  However,
even rapidly training models will take multiple hours.
Using random forest requires up to 24 hours on new 
Apple silicon.

## For more information

* "[Supervised learning of gene-regulatory networks based on graph distance profiles of transcriptomics data
](https://www.nature.com/articles/s41540-020-0140-1)"
* James P. Howard, II <<jh@jameshoward.us>>

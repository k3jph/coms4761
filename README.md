
# Computational Genomics Final Project


<div align="center">

[![PyPI - Version](https://img.shields.io/pypi/v/coms4761.svg)](https://pypi.python.org/pypi/coms4761)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/coms4761.svg)](https://pypi.python.org/pypi/coms4761)
[![Tests](https://github.com/howardjp/coms4761/workflows/tests/badge.svg)](https://github.com/howardjp/coms4761/actions?workflow=tests)
[![Codecov](https://codecov.io/gh/howardjp/coms4761/branch/main/graph/badge.svg)](https://codecov.io/gh/howardjp/coms4761)
[![Read the Docs](https://readthedocs.org/projects/coms4761/badge/)](https://coms4761.readthedocs.io/)
[![PyPI - License](https://img.shields.io/pypi/l/coms4761.svg)](https://pypi.python.org/pypi/coms4761)

[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)


</div>


 Final project for COMS 4761


* GitHub repo: <https://github.com/k3jph/coms4761.git>
* Free software: MIT


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

**Implementation**

## Credits

This package was created with [Cookiecutter][cookiecutter] and the [fedejaure/cookiecutter-modern-pypackage][cookiecutter-modern-pypackage] project template.

[cookiecutter]: https://github.com/cookiecutter/cookiecutter
[cookiecutter-modern-pypackage]: https://github.com/fedejaure/cookiecutter-modern-pypackage

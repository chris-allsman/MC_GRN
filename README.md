# MetaCell-based Ensemble Construction of Gene Regulatory Networks

This is supplemental code for our COMS4761 project at UC Berkeley. All code here is written by Chris Allsman, Andrey Zaznaev, and Jae Shin, with packages utilized cited at the end of the README.

Our code has two components: a collection of R scripts used to generate MetaCell samples, and a collection of Python functions for generating and analyzying GRN from those graphs. It is possible to use these components "end-to-end", but for testing purposes we recommend considering them separately.

## Environment
Due to package dependencies, the R scripts are only compatible with MacOS or Linux systems. We ran into difficulties running the scripts on Ubuntu (or more specifically, WSL2), so we recommend running on MacOS if possible. In particular, code was run on MacOS 12.3.

There are no such operating system requirements for the Python module, although it was only tested on Windows 11.

The language versions we used were:
  - R version 4.2.1
  - Python version 3.10.7

There are some package dependencies that need to be installed for Python, and are given in `python/requirements.txt` and listed below. All requirements can be installed via `pip install -r /path/to/requirements.txt` from the base directory. 

- GReNaDIne version 0.0.21
- matplotlib version 3.6.2
- networkx version 2.8.8
- numpy version 1.23.4
- pandas version 1.5.1
- scipy version 1.10.1

## Running Files
In order to generate MetaCell assignments from R scripts, they must be run interactively and sequentially. The diagram below gives the workflow that must be follows.

![image](https://user-images.githubusercontent.com/19377828/236118265-7c1cf613-ced9-4649-939f-d1bfeca3c5c4.png)

We have provided example output from the R scripts in the Python module, so the Python component can be run entirely independently from the R scripts. The python module can be executed by simply executing `main.py` from the root directory of the module, that is, `./python`.

It is worth noting that the input files have been shrunken considerably (by removing a random sample of non-transcription factor genes from the counts matrices) and the parameterization for the GENIE3 graph construction algorithm are set to construct a small number of trees with a low depth. Therefore, the results of this script will differ considerably from those presented in our paper. Replacing the input count files with the original data sources, the parameterizations with those presented in our paper, and the bulk/single cell CSV files with assignments generated from R scripts should give comporable results.

Although the count files were reduced in size, they are still somewhat large (on the order of 50 MB) in order for some amount of signal to be captured, therefore running the Python module will not be instantaneous. On our hardware (12 GB RAM, GPU not utilized) it took about 160 seconds to execute.

## Original Data Sources
The mouse counts were derived from https://singlecell.broadinstitute.org/single_cell/study/SCP1711/mouse-colon-stroma-inflammation#study-summary and the yeast counts and gold standard network were derived from https://zenodo.org/record/5272314.
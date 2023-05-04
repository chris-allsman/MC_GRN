# MetaCell-based Ensemble Construction of Gene Regulatory Networks

This is supplemental code for our COMS4761 project at Columbia University. All code here is written by Chris Allsman, Andrey Zaznaev, and Jae Shin, with packages utilized cited at the end of the README.

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

There are some package dependencies that need to be installed for R, which are listed below. All requirements can be installed via running the installation lines of code at the top of each R script needing these libraries. 

- Seurat version 4.3.0
- Matrix version 1.5.4
- Metacell version 0.3.7
- dplyr version 1.1.2


## Running Files
### R Workflow
In order to generate MetaCell assignments from R scripts, they must be run interactively and sequentially. The diagram below gives the workflow that must be follows.

![image](https://user-images.githubusercontent.com/19377828/236118265-7c1cf613-ced9-4649-939f-d1bfeca3c5c4.png)

### Demo Files for Quick Testing of Metacells
We provide two counts files with reduced geneset for quick testing of the metacell generation script (`generate_metacells.R`), located in `R`>`demo_data`. `mouse_counts_demo.txt` and `yeast_counts_demo.txt` are reduced counts matrices for mouse and yeast scRNAseq datasets, respectively. These files can be directly used with `generate_metacells.R` to generate the mapping csv file. To accomplish that, import counts matrix filename should be changed to either `mouse_counts_demo.txt` or `yeast_counts_demo.txt`.

### Running Python

We have provided example output from the R scripts in the Python module, so the Python component can be run entirely independently from the R scripts. The python module can be executed by simply executing `main.py` from the root directory of the module, that is, `./python`.

It is worth noting that the input files have been shrunken considerably (by removing a random sample of non-transcription factor genes from the counts matrices) and the parameterization for the GENIE3 graph construction algorithm are set to construct a small number of trees with a low depth. Therefore, the results of this script will differ considerably from those presented in our paper. Replacing the input count files in `./python/resources` with the original data sources, GENIE3 parameters with those presented in our paper, and the bulk/single cell CSV files in `./python/resources` with assignments generated from R scripts should give comporable results. 

Although the count files were reduced in size, they are still somewhat large (on the order of 50 MB) in order for some amount of signal to be captured, therefore running the Python module will not be instantaneous. On our hardware (12 GB RAM, GPU not utilized) it took about 160 seconds to execute.

## Original Data Sources
The mouse counts were derived from https://singlecell.broadinstitute.org/single_cell/study/SCP1711/mouse-colon-stroma-inflammation#study-summary and the yeast counts and gold standard network were derived from https://zenodo.org/record/5272314.

## Package citations

1. Harris, C. R. et al. Array programming with NumPy. Nature 585, 357–362 (2020).
2. McKinney, W. Data Structures for Statistical Computing in Python. in 56–61 (2010). doi:10.25080/Majora-92bf1922-00a.
3. Hagberg, A., Swart, P. & S Chult, D. Exploring network structure, dynamics, and function using networkx. https://www.osti.gov/biblio/960616 (2008).
4. Schmitt, P. et al. GReNaDIne: A Data-Driven Python Library to Infer Gene Regulatory Networks from Gene Expression Data. Genes 14, 269 (2023).
5. Hunter, J. D. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering 9, 90–95 (2007).
6. Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods 17, 261–272 (2020).
7. Hao, Y. et al. Integrated analysis of multimodal single-cell data. Cell 184, 3573-3587.e29 (2021).
8. Matrix: Sparse and Dense Matrix Classes and Methods. https://cran.r-project.org/web/packages/Matrix/index.html
9. Wickham H, François R, Henry L, Müller K, Vaughan D (2023). dplyr: A Grammar of Data Manipulation. https://dplyr.tidyverse.org, https://github.com/tidyverse/dplyr.
10. Baran, Y. et al. MetaCell: analysis of single-cell RNA-seq data using K-nn graph partitions. Genome Biol 20, 206 (2019).

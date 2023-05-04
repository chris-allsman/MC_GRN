Scripts used to generate inputs for the Python module. See the README in the root directory for more information about workflow.

`demo_data` contains small count files that can be used to generate MetaCells.

`sample_output` contains two CSVs corresponding to MetaCell assignments: one for mouse samples using a lower threshold and corresponding to a "single-cell-like" sample, and one for yeast using a higher threshold and corresponding to a "Pseudo-bulk-like" sample.

`generate_metacells.R` is responsible for generating the actual MetaCell assignments, while `prepare_mouse_scRNAseq.R` and `prepare_yeast_scRNAseq.R` are used to transform the counts matrices for use in `generate_metacells`.
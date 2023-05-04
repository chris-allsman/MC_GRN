Files specific to the Python module, which are used to generate and analyze GRNs. 

`resources` contains sources for constructing graphs, including:
- `*_bulk.csv` and `*_sc.csv`, which are the output files of the R component which assign MetaCell clusters to samples at two levels of granularity
- `*_counts.txt`, which are the original scRNA-seq counts matrices
- `*_tfs.txt`, which is a newline separated list of gene names correspondign to transcription factors
- `yeast_gs.tsv`, which is a "gold standard" GRN for the yeast network used for evaluation

`output` contains a few examples of GRNs output by the algorithm

`utils` contains utility functions used to run the analyses and visualization.
- `setup.py` for loading in counts dataframes and metacell assignments
- `eval.py` for getting the Jaccard index across ranked graphs, constructing and getting properties from a networkx representation of the GRN, and getting the GS evaluation result
-  `graph.py` for constructing GRN score matrices
-  `visualize.py` for visualizing networkx representations of the graph

`main.py` has three functions corresponding to the different results in our paper:
- `gold_standard_evaluation` for constructing a yeast GRN and comparing to gold standard
- `jaccard_index` for generating graphs from different metacell assignments and comparing similarity
- `visualize_graph` for constructing networkx graphs, printing properties, and rendering visualizations

# Filename: generate_metacells.R
# Authors: Chris Allsman, Jaeweon Shin, Andrey Zaznaev
# Last edit date: Apr 2023
# Purpose: calculates metacells from an input of single cell RNA-seq count 
#   matrix using Metacell package in R. Outputs csv file with metacell mappings,
#   as well as some accompanying figures.
# Code source: this code was adapted from the Metacell package tutorial:
#   https://tanaylab.github.io/metacell/articles/a-basic_pbmc8k.html

# Install and load necessary packages
if (!require("BiocManager")) install.packages('BiocManager') 
BiocManager::install("tanaylab/metacell")
install.packages("dplyr")

library("metacell")
library(dplyr)

# Initialize the database and figures directories
if(!dir.exists("db")) dir.create("db/")
scdb_init("db/", force_reinit=T)
if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")

# Load the scRNAseq counts matrix (tab-delimited format)
#   for mouse, use counts_matrix_H2O_Distal_all_genes.txt
#   for yeast, use yeast_counts.txt
mcell_import_scmat_tsv("sc_mouse_colon",
                       fn = "./counts_matrix_H2O_Distal_all_genes.txt",
                       dset_nm = "dataset_name")
mat = scdb_mat("sc_mouse_colon")

# Figure: plot UMIs per cell
mcell_plot_umis_per_cell("sc_mouse_colon")


mcell_add_gene_stat(gstat_id="sc_mouse_colon", mat_id="sc_mouse_colon", force=T)

mcell_gset_filter_varmean(gset_id="sc_mouse_colon_feats", 
                          gstat_id="sc_mouse_colon", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "sc_mouse_colon_feats", 
                      gstat_id="sc_mouse_colon", T_tot=100, T_top3=2)

mcell_plot_gstats(gstat_id="sc_mouse_colon", gset_id="sc_mouse_colon_feats")

mcell_add_cgraph_from_mat_bknn(mat_id="sc_mouse_colon", 
                               gset_id = "sc_mouse_colon_feats", 
                               graph_id="sc_mouse_colon_graph",
                               K=100,
                               dsamp=T)

# Change min_mc_size parameter to achieve desired number of metacells
mcell_coclust_from_graph_resamp(
  coc_id="sc_mouse_colon_coc500", 
  graph_id="sc_mouse_colon_graph",
  min_mc_size=100, 
  p_resamp=0.75, n_resamp=500)

# === metacell numbers ===
# yeast
  # min_mc_size=20 => 434 mcs
  # mic_mc_size=200 => 22 mcs
# mouse
  # min_mc_size=20 => 132 mcs
  # min_mc_size=40 => 105 mcs
  # min_mc_size=100 => 52 mcs
  # min_mc_size=150 => 27 mcs
  # min_mc_size=200 => 9 mcs
# ===                 ===

mcell_mc_from_coclust_balanced(
  coc_id="sc_mouse_colon_coc500", 
  mat_id= "sc_mouse_colon",
  mc_id= "sc_mouse_colon_mc", 
  K=30, min_mc_size=30, alpha=2)

# Ensure all metacells are homogeneous
mcell_plot_outlier_heatmap(mc_id="sc_mouse_colon_mc", 
                           mat_id = "sc_mouse_colon", T_lfc=3)

# Split metacells
mcell_mc_split_filt(new_mc_id="sc_mouse_colon_mc_f", 
                    mc_id="sc_mouse_colon_mc", 
                    mat_id="sc_mouse_colon",
                    T_lfc=3, plot_mats=F)

mcell_gset_from_mc_markers(gset_id="sc_mouse_colon_markers", 
                           mc_id="sc_mouse_colon_mc")

# Get mc object
mc = scdb_mc("sc_mouse_colon_mc")
print(max(mc@mc)) # print number of metacells

# Export the csv mapping file
mapping = mc@mc
write.csv(mapping,file="mapping.csv",row.names=T)

# Additional figures
mcell_mc_plot_marks(mc_id="sc_mouse_colon_mc", gset_id="sc_mouse_colon_markers", 
                    mat_id="sc_mouse_colon")
mcell_mc2d_force_knn(mc2d_id="sc_mouse_colon_2dproj",mc_id="sc_mouse_colon_mc", 
                     graph_id="sc_mouse_colon_graph")
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="sc_mouse_colon_2dproj")

# Visualizing the MC confusion matrix
mc_hc = mcell_mc_hclust_confu(mc_id="sc_mouse_colon_mc", 
                      graph_id="sc_mouse_colon_graph")
mc_sup = mcell_mc_hierarchy(mc_id="sc_mouse_colon_mc",
                            mc_hc=mc_hc, T_gap=0.04)
mcell_mc_plot_hierarchy(mc_id="sc_mouse_colon_mc", 
                        graph_id="sc_mouse_colon_graph", 
                        mc_order=mc_hc$order, 
                        sup_mc = mc_sup, 
                        width=2800, height=2000, min_nmc=2)
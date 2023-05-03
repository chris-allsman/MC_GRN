from grenadine.Inference.inference import score_links
from grenadine.Inference.regression_predictors import GENIE3
from grenadine.Preprocessing.standard_preprocessing import z_score 
import pandas as pd
import numpy as np
from scipy.stats import beta
from utils.setup import *

"""
Utilities for creating and combining graphs
"""

def run_graph(counts_df, assignment_df, tf_file_name, n_estimators, depth, save=False, name=""):
  """
    Creates a co-expression score matrix for the provided counts matrix and transcription factors, 
    normalizing the counts matrix and running the GENIE3 algorithm with the provided parameters.
  """
  grouped_counts_df = get_grouped_counts_df(counts_df, assignment_df)
  tfs = get_tfs(grouped_counts_df, tf_file_name)
  tf = pd.read_csv(tf_file_name, header=None)
  tf = tf[tf[0].isin(tfs)][0]


  GENIE3_params = {"n_estimators": n_estimators,
                  'max_depth': depth}

  X = z_score(grouped_counts_df,axis=1)
  score_matrix=score_links(X, GENIE3, tf, **GENIE3_params)
  if save:
    score_matrix.to_csv(name)
  return score_matrix

def ensemble_score_links(score_links_matrices, weight_functions=None):
    """
        Given matrices to combine and a function to apply to the given weights,
        produces a merged matrix according to the procedure described in the paper.
    """

    if weight_functions == None:
       weight_functions = [lambda x: x for _ in score_links_matrices]

    assert len(score_links_matrices) == 2, "Can only support averaging two GRNs"

    single_cell_grn, bulk_grn = score_links_matrices[1], score_links_matrices[0]

    score_links = weight_functions[1](np.clip(single_cell_grn - bulk_grn, 0, 1))
    score_links += weight_functions[0](bulk_grn)

    return score_links

default_weight_bulk = lambda x: .5 * beta.cdf((x * (x > .01)), .05, .5)
default_weight_sc = lambda x: x


import pandas as pd

"""
Utilities for preparing files for further processing given MetaCell assignments
"""

def get_cluster_df(name):
  """
    Given the name of a file assigning samples to MetaCells, prepare a Pandas DataFrame for
    further use
  """
  cluster_df = pd.read_csv(name)
  cluster_df = cluster_df.rename(columns={"Unnamed: 0": "Name", "x": "Cluster"})
  return cluster_df

def get_grouped_counts_df(counts_df, cluster_df, bootstrap=False, random_state=None):
  """
    Given a DataFrame with raw counts and a DataFrame of MetaCell assignments, produce
    a new DataFrame aggregating counts across samples in each MetaCell.

    Optionally supports bootstrap assignment.
  """
  grouped_counts_df = pd.DataFrame(index=counts_df.index)
  if bootstrap:
    sample_df = cluster_df.sample(frac=1, replace=True, random_state=random_state)
  else:
    sample_df = cluster_df
  for index, row in sample_df.iterrows():
      cluster = row["Cluster"]
      name = row["Name"]
      if cluster in grouped_counts_df.columns:
          grouped_counts_df[cluster] = grouped_counts_df[cluster] + counts_df[name]
      else:
          grouped_counts_df[cluster] = counts_df[name]
  return grouped_counts_df

def get_tfs(grouped_counts_df, tf_file_name):
  """
    Given the name of a flat file containing transcription factor names and a DataFrame
    associating counts with MetaCells, produces a list of transcription factor present
    in the DataFrame
  """
  with open(tf_file_name, 'r') as tf_file:
    tfs = tf_file.read().splitlines()
  tf_set = set(tfs)
  tfs = list(tf_set.intersection(set(grouped_counts_df.index)))
  return tfs
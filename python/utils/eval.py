import pandas as pd
import networkx as nx
from grenadine.Evaluation.evaluation import rank_GRN, jaccard_similarity, rank_to_networkx

"""
Methods for network evaluation
"""

def prepare_gs_result(gs_file_name):
    gs_df = pd.read_csv(gs_file_name, index_col=0, sep="\t")
    gs_rank = rank_GRN(gs_df)
    gs_edge_list = gs_rank.drop(["rank"], axis=1).rename(columns={"score": "IS_REGULATED"})
    gs_edge_list["IS_REGULATED"] = gs_edge_list["IS_REGULATED"].abs().astype("int")
    return gs_edge_list

def get_jaccard_index(dataframes, descriptors, n):
    """
        Given a list of rank dataframes (output from rank_GRN) and a list of sample names,
        returns the Jaccard index for the top n edges across the networks.
    """
    assert len(dataframes) >= 2, "Need at least 2 dataframes"
    processed_dfs = [df[["rank"]].rename(columns={"rank": d}) for df, d in zip(dataframes, descriptors)]
    joined = processed_dfs[0]
    for df in processed_dfs[1:]:
        joined = joined.join(df, how="outer")
    return jaccard_similarity(joined, k=n)

def make_graph_from_score_df(score_dataframe, n_edges):
    GRN = rank_GRN(score_dataframe)
    G = rank_to_networkx(GRN, top_n = n_edges, to_undirected=True)
    return G

def get_sccs(G, top=100):
    return sorted(nx.connected_components(G), key=lambda x: -len(x))[:top]

def get_communities(G, threshold=1, random_state=1):
    communities = nx.community.louvain_communities(G, seed=random_state)
    return list(set.union(*list(filter(lambda x: len(x) >= threshold, communities))))
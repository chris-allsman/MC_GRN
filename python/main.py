import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from utils.setup import *
from utils.eval import *
from utils.graph import *
from utils.visualize import *
from grenadine.Evaluation.evaluation import evaluate_result, rank_GRN
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

import time

start_time = time.time()


################################
# Do evaluation for GS dataset #
################################
def gold_standard_evaluation():
    np.random.seed(0)
    print("\nConstructing yeast graphs")
    yeast_counts = pd.read_csv("resources/yeast_counts.txt", sep="\t")
    yeast_tf_filename = "resources/yeast_tf.tsv"
    yeast_bulk_assignment = get_cluster_df("resources/yeast_bulk.csv")
    yeast_sc_assignment = get_cluster_df("resources/yeast_sc.csv")
    yeast_bulk_scores = run_graph(yeast_counts, yeast_bulk_assignment, yeast_tf_filename, 3, 2)
    yeast_sc_scores = run_graph(yeast_counts, yeast_sc_assignment, yeast_tf_filename, 3, 2)
    ensemble_scores = ensemble_score_links([yeast_bulk_scores.copy(), yeast_sc_scores.copy()], 
                                        [default_weight_bulk, default_weight_sc])

    print("\nComparing to gold standard network:")
    yeast_gs_scores = prepare_gs_result("resources/yeast_gs.tsv")

    print("\nMetrics for combined graph:")
    print(evaluate_result(ensemble_scores, yeast_gs_scores, n_links=1000))

    print("\nMetrics for SC graph with no weighting or averaging:")
    print(evaluate_result(yeast_sc_scores, yeast_gs_scores, n_links=1000))

    print("\nMetrics for bulk graph only with nonlinear transformation:")
    print(evaluate_result(pd.DataFrame(default_weight_bulk(yeast_bulk_scores.copy()), columns=yeast_bulk_scores.columns, index=yeast_bulk_scores.index), yeast_gs_scores, n_links=1000))

def jaccard_index():
    np.random.seed(44)
    print("\nConstructing mouse graphs")
    mouse_counts = pd.read_csv("resources/mouse_counts.txt", sep="\t").set_index('GeneSymbol')
    mouse_tf_filename = "resources/mouse_tfs.txt"
    mouse_bulk_1 = get_cluster_df("resources/mouse_bulk_1.csv")
    mouse_bulk_2 = get_cluster_df("resources/mouse_bulk_2.csv")
    mouse_sc_1 = get_cluster_df("resources/mouse_sc_1.csv")
    mouse_sc_2 = get_cluster_df("resources/mouse_sc_2.csv")
    ensemble_1_scores = ensemble_score_links([run_graph(mouse_counts, mouse_bulk_1.copy(), mouse_tf_filename, 3, 2), run_graph(mouse_counts, mouse_sc_1.copy(), mouse_tf_filename, 3, 2)], 
                                        [default_weight_bulk, default_weight_sc])
    ensemble_2_scores = ensemble_score_links([run_graph(mouse_counts, mouse_bulk_2.copy(), mouse_tf_filename, 3, 2), run_graph(mouse_counts, mouse_sc_2.copy(), mouse_tf_filename, 3, 2)], 
                                        [default_weight_bulk, default_weight_sc])
    ensemble_3_scores = ensemble_score_links([run_graph(mouse_counts, mouse_bulk_1.copy(), mouse_tf_filename, 3, 2), run_graph(mouse_counts, mouse_sc_2.copy(), mouse_tf_filename, 3, 2)], 
                                        [default_weight_bulk, default_weight_sc])
    print("Jaccard Index:\n")
    print(get_jaccard_index([rank_GRN(ensemble_1_scores), rank_GRN(ensemble_2_scores), rank_GRN(ensemble_3_scores)], ["ensemble 1", "ensemble 2", "ensemble 3"], 100))

def visualize_graph():
    np.random.seed(44)
    print("\nConstructing mouse graphs")
    mouse_counts = pd.read_csv("resources/mouse_counts.txt", sep="\t").set_index('GeneSymbol')
    mouse_tf_filename = "resources/mouse_tfs.txt"
    mouse_bulk = get_cluster_df("resources/mouse_bulk_1.csv")
    mouse_sc = get_cluster_df("resources/mouse_sc_1.csv")
    bulk_scores = run_graph(mouse_counts, mouse_bulk, mouse_tf_filename, 3, 2)
    sc_scores = run_graph(mouse_counts, mouse_sc, mouse_tf_filename, 3, 2)
    ensemble_scores = ensemble_score_links([bulk_scores.copy(), sc_scores.copy()], 
                                        [default_weight_bulk, default_weight_sc])
    
    G_ensemble = make_graph_from_score_df(ensemble_scores, 500)
    ensemble_communities = get_communities(G_ensemble)
    print("Maximum community size for ensemble:", max([len(c) for c in ensemble_communities]))
    print("Average community size for ensemble:", sum([len(c) for c in ensemble_communities]) / len(ensemble_communities))
    plot_closeness_centrality(G_ensemble, .1)

    d = {
     "Col1a1": ["Fibroblast", "Endothelial cellthelial cell of lymphatic vessel", "Interstitial cell", "Glial Cell", "Pericyte cell"], 
     "Col3a1": ["Fibroblast", "Endothelial cellthelial cell of lymphatic vessel", "Interstitial cell", "Endothelial cell", "Glial Cell", "Pericyte cell"], 
     "Dcn": ["Endothelial cell", "Glial Cell", "Pericyte cell"], 
     "Acta2": ["Interstitial cell"], 
     "Cnn1": ["Smooth Muscle Cell"], 
     "Mylk": ["Smooth Muscle Cell", "Interstitial cell"],
     "Fn1": ["Interstitial cell", "Fibroblast"],
     "Dan": ["Fibroblast"],
     "Loxl1": ["Pericyte cell"],
     "Mmp2": ["Fibroblast"],
     "Loxl2": ["Fibroblast"],
     "Pgdfra": ["Fibroblast"],
     "Acta2": ["Smooth Muscle Cell", "Endothelial cell", "Interstitial cell", "Glial Cell"],
     "Pecam1": ["Endothelial cell"],
     "Kdr": ["Endothelial cell"],
     "Lyve1": ["Endothelial cell"],
     "Kit": ["Interstitial cell"],
     "Ano1": ["Interstitial cell"],
     "Res": ["Interstitial cell"],
     "Pecam1": ["Endothelial cell"],
     "Plvap": ["Endothelial cell"],
     "Mmp2": ["Glial Cell"],
     "S100b": ["Glial Cell"],
     "Rgs5": ["Pericyte cell"]
     }

    cmap = {
        "Smooth Muscle Cell": "red",
        "Fibroblast": "green",
        "Endothelial cell": "orange",
        "Endothelial cellthelial cell of lymphatic vessel": "orange",
        "Interstitial cell": "purple",
        "Glial Cell": "blue",
        "Pericyte cell": "pink"
    }

    plot_marker_genes(G_ensemble, d, cmap)

    G_sc = make_graph_from_score_df(sc_scores, 500)
    sc_communities = get_communities(G_sc)
    print("Maximum community size for baseline graph:", max([len(c) for c in sc_communities]))
    print("Average community size for baseline graph:", sum([len(c) for c in sc_communities]) / len(sc_communities))


gold_standard_evaluation()
jaccard_index()
visualize_graph()

print("--- %s seconds ---" % (time.time() - start_time))



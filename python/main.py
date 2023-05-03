from utils.setup import *
from utils.eval import *
from utils.graph import *
from utils.visualize import *
from grenadine.Evaluation.evaluation import evaluate_result


################################
# Do evaluation for GS dataset #
################################
print("Constructing yeast graphs:")
yeast_counts = pd.read_csv("resources/yeast_count.txt")
yeast_tf_filename = "resources/yeast_tf.txt"
yeast_bulk_assignment = pd.read_csv("resources/yeast_bulk_assignment")
yeast_sc_assignment = pd.read_csv("resources/yeast_sc_assignment")
yeast_bulk_scores = run_graph(yeast_counts, yeast_bulk_assignment, yeast_tf_filename, 5, 2)
yeast_sc_scores = run_graph(yeast_counts, yeast_sc_assignment, yeast_tf_filename, 5, 2)
ensemble_scores = ensemble_score_links([yeast_bulk_scores, yeast_sc_scores], 
                                       [default_weight_bulk, default_weight_sc])

print("Comparing to gold standard network:")
yeast_gs_scores = prepare_gs_result("resources/yeast_gs_network")

print("Metrics for combined graph:")
print(evaluate_result(ensemble_scores, yeast_gs_scores, n_links=1000))

print("Metrics for SC graph with no weighting or averaging:")
print(evaluate_result(ensemble_scores, yeast_sc_scores, n_links=1000))



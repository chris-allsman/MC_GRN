from grenadine.Visualisation import visualisation
from matplotlib import pyplot as plt
import networkx as nx

def plot_closeness_centrality(G, size):
    closeness = visualisation.closeness_centrality(G)
    visualisation.draw_by(G, closeness, size)
    plt.title("Network visualized by closeness centrality")
    plt.show()


def plot_marker_genes(G, marker_gene_dict, cmap, display_threshold=3):
    myNet = G.copy()
    myNet.remove_nodes_from(list(set.union(*list(filter(lambda x: len(x) < display_threshold, nx.connected_components(G))))))
    colnodes = list()
    sizenodes = list()
    for node in myNet.nodes:
        if node in marker_gene_dict:
            colnodes.append(cmap[marker_gene_dict[node][0]])
            sizenodes.append(50)
        else:
            colnodes.append("black")
            sizenodes.append(2)
    nx.draw(myNet, pos = nx.nx_agraph.graphviz_layout(G), node_color = colnodes, node_size = sizenodes)
    pairs = [[k, v] for k, v in cmap.items()]
    for p in pairs:
        plt.scatter([],[], c=[p[1]], label=p[0])

    plt.legend()
    plt.title("Marker genes in graph")
    plt.show()


import networkx as nx
import random
def network_mote_carlo_points(G, weightID, number_of_nodes):
    edge_weights = nx.get_edge_attributes(G, weightID)
    #edge_oneway = nx.get_edge_attributes(G, "oneway")
    L = 0
    for edge in G.edges:
        L = L + edge_weights[edge]/2
    edge_weights_percentage = edge_weights/L
    random_draws = random.choices(G.edges, weights = edge_weights_percentage, k = number_of_nodes)
    for 
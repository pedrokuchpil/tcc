import networkx as nx
import numpy as np
import random

# Step 1: Generate a random network
def generate_network(num_nodes, edge_prob):
    # Create a random Erdős-Rényi graph with num_nodes and edge probability edge_prob
    G = nx.erdos_renyi_graph(num_nodes, edge_prob)
    return G

# Step 2: Generate attributes based on distributions (exponential, gaussian, gamma)
def generate_node_attributes(G, distribution='gaussian'):
    attributes = {}
    if distribution == 'gaussian':
        attributes = {node: np.random.normal(0, 1) for node in G.nodes()}
    elif distribution == 'gamma':
        attributes = {node: np.random.gamma(2, 2) for node in G.nodes()}
    elif distribution == 'exponential':
        attributes = {node: np.random.exponential(1) for node in G.nodes()}
    nx.set_node_attributes(G, attributes, 'attribute')
    return attributes

# Step 3: Assign edge weights based on a distribution
def generate_edge_weights(G, distribution='gaussian'):
    weights = {}
    if distribution == 'gaussian':
        weights = {(u, v): np.random.normal(0, 1) for u, v in G.edges()}
    elif distribution == 'gamma':
        weights = {(u, v): np.random.gamma(2, 2) for u, v in G.edges()}
    elif distribution == 'exponential':
        weights = {(u, v): np.random.exponential(1) for u, v in G.edges()}
    nx.set_edge_attributes(G, weights, 'weight')
    return weights

# Step 4: Generate ground truth (community assignment)
def assign_communities(G, num_communities):
    community_assignment = {}
    for node in G.nodes():
        community_assignment[node] = random.randint(0, num_communities-1)  # Assign random community
    nx.set_node_attributes(G, community_assignment, 'community')
    return community_assignment

# Step 5: Write the network and attributes to a file
def write_to_file(G, filename):
    with open(filename, 'w') as f:
        # Write node attributes and communities
        f.write("# Node attributes (id,attribute,community)\n")
        for node, data in G.nodes(data=True):
            f.write(f"NODE,{node},{data['attribute']},{data['community']}\n")
        
        # Write edge weights
        f.write("\n# Edge weights (source,target,weight)\n")
        for u, v, data in G.edges(data=True):
            f.write(f"EDGE,{u},{v},{data['weight']}\n")


# Main function
def main():
    num_nodes = 10            # Number of nodes in the network
    edge_prob = 0.3           # Probability of edge creation
    num_communities = 3       # Number of communities

    # Step 1: Generate the network
    G = generate_network(num_nodes, edge_prob)

    # Step 2: Generate node attributes using a specified distribution
    node_attributes = generate_node_attributes(G, distribution='gaussian')

    # Step 3: Generate edge weights using a specified distribution
    edge_weights = generate_edge_weights(G, distribution='gaussian')

    # Step 4: Assign communities (ground truth)
    communities = assign_communities(G, num_communities)

    # Step 5: Write everything to a file
    write_to_file(G, "network_output.txt")

if __name__ == "__main__":
    main()

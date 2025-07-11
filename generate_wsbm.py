import numpy as np
import networkx as nx

class WSBMGenerator:
    def __init__(self, weight_distribution, weight_variance, min_, max_):
        self.n_clusters = n_clusters = 3   # Number of clusters (communities)
        nodes_per_cluster = 5
        n = int(n_clusters*nodes_per_cluster)
        pout = 0.05
        pin = 0.1
        p = (pin - pout) * np.eye( n_clusters ) + pout * np.ones( (n_clusters, n_clusters) )
        sizes = int(n//n_clusters)
        self.probability_matrix=p # Block connectivity probabilities
        self.communities_sizes = [sizes]*n_clusters
        self.weight_distribution = weight_distribution  # Weight distribution function
        self.weight_variance = weight_variance  # Variance for weight distribution
        self.min_ = min_  # Minimum weight for distribution
        self.max_ = max_  # Maximum weight for distribution
        self.weight_centers = None  # To store weight centers

    # Helper function to calculate weight parameters
    def get_w_params(self, weight_centers, weight_variance, n_clusters):
        params = np.zeros((n_clusters, n_clusters, 2))  # Assuming mean and variance for each block pair
        for i in range(n_clusters):
            for j in range(n_clusters):
                params[i, j] = [weight_centers[i, j], weight_variance]
        return params

    # Function to generate WSBM
    def generate_WSBM(self, complete_graph=False):
        N = np.sum(self.communities_sizes)

        # Generate binary connectivity matrix using the stochastic block model
        if complete_graph:
            G = nx.stochastic_block_model(self.communities_sizes, np.ones((self.n_clusters, self.n_clusters)), directed=False, seed=42)
        else:
            G = nx.stochastic_block_model(self.communities_sizes, self.probability_matrix, directed=False, seed=42)

        # Draw the means of the weight distributions for each pair of community interaction
        if self.weight_centers is None:
            self.weight_centers = np.zeros((self.n_clusters, self.n_clusters))
            self.weight_centers[np.triu_indices(self.n_clusters, k=0)] = \
                np.linspace(self.min_, self.max_, num=int(self.n_clusters * (self.n_clusters + 1) / 2))
            self.weight_centers = self.weight_centers + self.weight_centers.T - np.diag(np.diag(self.weight_centers))

        # Get weight parameters for the weight distributions
        params = self.get_w_params(self.weight_centers, self.weight_variance, self.n_clusters)

        # Assign weights to edges based on the weight distribution and the block matrix
        for i, j in G.edges():
            block_i = G.nodes[i]["block"]  # Community of node i
            block_j = G.nodes[j]["block"]  # Community of node j
            weight_params = params[block_i][block_j]  # Get the parameters for the weight distribution
            G[i][j]['weight'] = self.weight_distribution(weight_params[0], weight_params[1])  # Apply weight distribution

        return G

    # Function to write the graph to a file in the NODE/EDGE format
    def write_to_file(self, G, filename):
        with open(filename, 'w') as f:
            # Write node attributes
            f.write("# Node attributes (id,attribute,community)\n")
            for node, data in G.nodes(data=True):
                f.write(f"NODE,{node},{np.random.normal(0, 1)},{data['block']}\n")

            # Write edge weights
            f.write("\n# Edge weights (source,target,weight)\n")
            for u, v, data in G.edges(data=True):
                f.write(f"EDGE,{u},{v},{data['weight']}\n")

            # Write ground truth (communities)
            f.write("\n# Ground Truth (id,community)\n")
            for node, data in G.nodes(data=True):
                f.write(f"GROUND_TRUTH,{node},{data['block']}\n")

# Example of usage
def main():
    # Example parameters
    weight_variance = 0.1  # Small variance for weights
    min_ = 1.0  # Minimum possible weight
    max_ = 5.0  # Maximum possible weight

    # Weight distribution: Gaussian with given mean and variance
    def gaussian_weight(mean, variance):
        return np.random.normal(mean, np.sqrt(variance))

    # Create WSBM generator
    generator = WSBMGenerator(gaussian_weight, weight_variance, min_, max_)

    # Generate the WSBM graph
    G = generator.generate_WSBM()

    # Write to file
    generator.write_to_file(G, "network_output.txt")

if __name__ == "__main__":
    main()

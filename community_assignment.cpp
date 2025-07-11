#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include <cmath>
#include <cassert>
#include <tuple>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <random>
#include <algorithm>
#include <numeric>
#include <set>

using namespace std;

struct node_t {
    int id;
    double weight;
    node_t* next;
};

enum class DistributionType { GAUSSIAN, EXPONENTIAL, POISSON };

double diverg_phi(double y, double nu, DistributionType distribution) {
    int theta = 1;

    switch (distribution) {
        case DistributionType::GAUSSIAN: // Gaussian
            return ( pow (theta , 2) / 2 )  * pow (y-nu, 2);
        case DistributionType::POISSON: // Gamma -- IMPLEMENT
            return 0;
        case DistributionType::EXPONENTIAL: // Exponential
            return y / nu - log (y/nu) - 1;
        default:
            return 0;
    }
}

double diverg_psi(double weight, double mu, DistributionType distribution) { // wrong
    switch (distribution) {
        case DistributionType::GAUSSIAN :
            return ( pow (weight - mu, 2) / 2 );
        case DistributionType::POISSON: // IMPLEMENT
            return 0;
        case DistributionType::EXPONENTIAL: // IMPLEMENT
            return 0;
        default:
            return 0;
    }
}

double diverg_KL(double p, double q) {
    if (q == 0 && p == 1) {
        return numeric_limits<double>::infinity();
    }
    else if (q == 0 && p == 0) {
        return 0;
    }
    else if (q == 1 && p == 1) {
        return 0;
    }
    else if (q == 1 && p == 0) {
        return numeric_limits<double>::infinity();
    }
    else if (p == 0 && q > 0) {
        return - log( (1-q) );
    }
    else if (p == 1 && q > 0) {
        return - log ( q );
    }
    else if (q == 0 && p > 0) { //HANDLE DIVISION BY ZERO
        return numeric_limits<double>::infinity();
    }
    else if (q == 1 && p > 0) { //HANDLE DIVISION BY ZERO
        return numeric_limits<double>::infinity();;
    }

    return p * log ( p/q ) + (1-p) * log( (1-p) / (1-q) );
}

void compute_maximum_assignment(
    vector<node_t*> A,           // adjacency list
    vector<double> Y,            // node attributes
    vector<unsigned int>& z,     // community assignment vector
    vector<vector<unsigned int>> edge_count,  // edge count matrix
    vector<vector<double>> edge_weight,       // edge weight matrix
    vector<unsigned int> node_count,               // node count per community
    vector<double> node_attribute,                 // node attribute per community
    vector<vector<unsigned int>> neighbor_count,  // neighbors count matrix
    vector<vector<double>> p,                 // edge probability matrix
    vector<vector<double>> mu,                // average edge weight matrix
    vector<double> nu,                             // average node attribute per community
    vector<unsigned int>& z_new,                   // new community assignment vector
    DistributionType distribution              // distribution type
) {
    int n = A.size();
    int k = p.size();

    // Loop through all nodes
    for (int i = 0; i < n; ++i) {
        double max_L_ia = -numeric_limits<double>::infinity();
        int max_z = k+1;

        // Loop through all communities for node i
        for (int a = 0; a < k; ++a) {
            // Compute the value of the likelihood function
            double L_ia = -2 * diverg_phi(Y[i], nu[a], distribution);
            // Go over the neighbors of node i to update psi
            node_t* v = A[i];
            while (v != nullptr) {
                L_ia -= diverg_psi(v->weight, mu[a][z[v->id]], distribution);
                v = v->next;
            }

            // Update KL of the edges and non-edges of i
            for (int b = 0; b < k; ++b) {
                // Update divergence with present edges in b
                if (neighbor_count[i][b] > 0) {
                    L_ia -= neighbor_count[i][b] * diverg_KL(1, p[a][b]);
                } else {
                    L_ia += 0; // Handle division by zero case
                }
                // Update divergence with missing edges in b
                L_ia -= (node_count[b] - neighbor_count[i][b]) * diverg_KL(0, p[a][b]);
            }

            // Divide total divergence by two
            L_ia /= 2;

            // Update the maximum value and community
            if (max_z == k+1 || L_ia > max_L_ia) {
                max_z = a;
                max_L_ia = L_ia;
            }
        }

        // Update new community for i
        z_new[i] = max_z;
    }

    // Update assignment vector for all nodes
    for (int i = 0; i < n; ++i) {
        z[i] = z_new[i];
    }
}


void update_parameters(
    vector<node_t*> A,           // adjacency list
    vector<double> Y,            // node attributes
    vector<unsigned int>& z,     // community assignment vector
    vector<vector<unsigned int>>& edge_count,  // edge count matrix
    vector<vector<double>>& edge_weight,       // edge weight matrix
    vector<unsigned int>& node_count,               // node count per community
    vector<double>& node_attribute,                 // node attribute per community
    vector<vector<unsigned int>>& neighbor_count,  // neighbors count matrix
    vector<vector<double>>& p,                 // edge probability matrix
    vector<vector<double>>& mu,                // average edge weight matrix
    vector<double>& nu                              // average node attribute per community
) {
    int n = A.size();
    int k = p.size();

    // Initialize auxiliary variables
    for (int a = 0; a < k; ++a) {
        node_count[a] = 0;
        node_attribute[a] = 0;
        for (int b = 0; b < k; ++b) {
            edge_count[a][b] = 0;
            edge_weight[a][b] = 0;
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int a = 0; a < k; ++a) {
            neighbor_count[i][a] = 0;
        }
    }
    // Main loop over all nodes of the graph
    for (int i = 0; i < n; ++i) {
        // Update number of nodes and total attribute
        node_count[z[i]] += 1;
        node_attribute[z[i]] += Y[i];
        // Go over the neighbors of node i
        node_t* v = A[i];
        while (v != nullptr) {
            neighbor_count[i][z[v->id]] += 1; // Update neighbor count
            // Consider each edge only once (skip edges already seen)
            if (v->id < i) {
                v = v->next;
                continue;
            }

            // Update edge count and weight
            edge_count[z[i]][z[v->id]] += 1;
            edge_weight[z[i]][z[v->id]] += v->weight;
            v = v->next;
        }
    }

    // Update model parameters
    for (int a = 0; a < k; ++a) {
        if (node_count[a] > 0) {
            nu[a] = node_attribute[a] / node_count[a];
        } else {
            nu[a] = 0.0; // Handle division by zero case
        }

        for (int b = 0; b < k; ++b) {
            // Add both counts from [a][b] and [b][a] since matrix is symmetric
            if (edge_count[a][b] + edge_count[b][a] > 0) {
                mu[a][b] = (edge_weight[a][b] + edge_weight[b][a]) / (edge_count[a][b] + edge_count[b][a]);
                p[a][b] = static_cast<double>(edge_count[a][b] + edge_count[b][a]) / (node_count[a] * node_count[b]);
            } else {  // Handle division by zero case
                mu[a][b] = 0.0;
                p[a][b] = 0.0; 
            }
        }
    }
}


vector<unsigned int> bregman_clustering(
    vector<node_t*> A,           // adjacency list
    vector<double> Y,            // node attributes
    vector<unsigned int> z,      // initial community assignment vector
    vector<vector<unsigned int>>& edge_count,  // edge count matrix
    vector<vector<double>>& edge_weight,       // edge weight matrix
    vector<unsigned int>& node_count,               // node count per community
    vector<double>& node_attribute,                 // node attribute per community
    vector<vector<unsigned int>>& neighbor_count,  // neighbors count matrix
    vector<vector<double>>& p,                 // edge probability matrix
    vector<vector<double>>& mu,                // average edge weight matrix
    vector<double>& nu,                             // average node attribute per community
    DistributionType distribution                            // distribution type
) {
    vector<unsigned int> z_new(z.size(), 0);  // Initialize z_new with zeros
    int max_iterations = 20;  // Maximum number of iterations
    int iterations = 0;

    bool converged = false;
    while (!converged) {
        // 3: Update parameters
        update_parameters(A, Y, z, edge_count, edge_weight, node_count, node_attribute, neighbor_count, p, mu, nu);

        // 5-9: Compute maximum assignment and update community memberships
        compute_maximum_assignment(A, Y, z, edge_count, edge_weight, node_count, node_attribute, neighbor_count, p, mu, nu, z_new, distribution);

        // 10: Check for convergence
        converged = (max_iterations <= iterations); // maximum number of iterations OR if there is no change OR number of changes is smaller than C we stop

        // Update z with z_new
        z = z_new;

        // Increment the number of iterations
        iterations++;
    }

    return z;
}
void add_edge(vector<node_t*>& A, int u, int v, double weight) {
    // Create a new node_t for v and add it to u's adjacency list
    node_t* new_node = new node_t{v, weight, nullptr};
    if (A[u] == nullptr) {
        A[u] = new_node;
    } else {
        node_t* current = A[u];
        while (current->next != nullptr) {
            current = current->next;
        }
        current->next = new_node;
    }
}

// Helper function to generate random numbers based on a Gaussian distribution
double generate_weight(double mean, double variance) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(mean, sqrt(variance));  // Gaussian distribution
    return d(gen);
}
    // Function to generate a random value based on a Gaussian distribution
double generate_gaussian(double mean, double variance) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(mean, sqrt(variance));
    return d(gen);
}

// Function to generate a random value based on an exponential distribution
double generate_exponential(double lambda) {
    random_device rd;
    mt19937 gen(rd());
    exponential_distribution<> d(lambda);
    return d(gen);
}

double generate_poisson(double mean) {
    random_device rd;
    mt19937 gen(rd());
    poisson_distribution<int> d(mean);
    return static_cast<double>(d(gen));
}

// Function to print the adjacency list of the graph
void print_graph(const vector<node_t*>& graph) {
    for (int i = 0; i < graph.size(); ++i) {
        cout << "Node " << i << " edges: ";
        node_t* current = graph[i];
        while (current != nullptr) {
            cout << " -> (Node " << current->id << ", Weight: " << current->weight << ")";
            current = current->next;
        }
        cout << endl;
    }
}

class WSBMGenerator {
public:
    WSBMGenerator(int n, int k, DistributionType attr_dist, DistributionType edge_dist, DistributionType weight_dist,
                  double attr_param1, double attr_param2, double edge_param1, double edge_param2,
                  double weight_param1, double weight_param2, double min_, double max_)
        : attribute_dist(attr_dist), edge_dist(edge_dist), weight_dist(weight_dist),
          attr_param1(attr_param1), attr_param2(attr_param2), edge_param1(edge_param1), edge_param2(edge_param2),
          weight_param1(weight_param1), weight_param2(weight_param2), min_(min_), max_(max_) {
        n_clusters = k;  // Number of clusters (communities)
        N = n;

        nodes_per_cluster = n / k;

        // Define intra- and inter-cluster probabilities
        pout = 0.15;
        pin = 0.3;

        // Generate the block probability matrix
        generate_probability_matrix();
        generate_weight_centers();
    }

    // Function to generate WSBM
    tuple<vector<node_t*>, vector<double>, vector<int>> generate_WSBM(bool complete_graph = false) {
        vector<node_t*> graph(N, nullptr);
        node_attributes.resize(N);
        vector<int> Z_true(N);

        // Assign node attributes based on the chosen distribution
        for (int i = 0; i < N; ++i) {
            node_attributes[i] = generate_value(attribute_dist, attr_param1, attr_param2);
            Z_true[i] = i / nodes_per_cluster;
        }

        // Simulate the stochastic block model
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                int block_i = i / nodes_per_cluster;
                int block_j = j / nodes_per_cluster;

                double base_prob = (complete_graph) ? 1.0 : probability_matrix[block_i][block_j];
                double modified_prob = base_prob * generate_value(edge_dist, edge_param1, edge_param2);

                if (random_prob() <= modified_prob) {
                    // Assign weight based on the block interactions
                    double mean = weight_centers[block_i][block_j];
                    double weight = generate_value(weight_dist, weight_param1, weight_param2);

                    // Add the edge to the adjacency list (undirected graph, so add both directions)
                    add_edge(graph, i, j, weight);
                    add_edge(graph, j, i, weight);
                }
            }
        }

        return make_tuple(graph, node_attributes, Z_true);
    }

    // Function to print node attributes
    void print_node_attributes() {
        for (int i = 0; i < node_attributes.size(); ++i) {
            cout << "Node " << i << ": Attribute = " << node_attributes[i] << endl;
        }
    }

private:
    int n_clusters;                    // Number of clusters
    int nodes_per_cluster;             // Nodes per cluster
    int N;                             // Total number of nodes
    double pout, pin;                  // Inter- and intra-cluster probabilities
    double attr_param1, attr_param2;   // Parameters for attribute distribution
    double edge_param1, edge_param2;   // Parameters for edge distribution
    double weight_param1, weight_param2; // Parameters for weight distribution
    double min_, max_;                 // Minimum and maximum weights for distributions
    vector<double> node_attributes;  // Node attributes
    vector<vector<double>> probability_matrix;  // Block probability matrix
    vector<vector<double>> weight_centers;      // Weight centers for each block
    DistributionType attribute_dist;   // Distribution type for node attributes
    DistributionType edge_dist;        // Distribution type for edge probabilities
    DistributionType weight_dist;      // Distribution type for edge weights

    // Helper function to generate the block probability matrix
    void generate_probability_matrix() {
        probability_matrix.resize(n_clusters, vector<double>(n_clusters, pout));
        for (int i = 0; i < n_clusters; ++i) {
            probability_matrix[i][i] = pin;
        }
    }

    // Helper function to generate the weight centers for the block interactions
    void generate_weight_centers() {
        weight_centers.resize(n_clusters, vector<double>(n_clusters, 0.0));
        int index = 0;
        for (int i = 0; i < n_clusters; ++i) {
            for (int j = i; j < n_clusters; ++j) {
                double value = min_ + (max_ - min_) * (index / double(n_clusters * (n_clusters + 1) / 2));
                weight_centers[i][j] = value;
                weight_centers[j][i] = value;  // Symmetric
                index++;
            }
        }
    }

    // Helper function to generate a random value based on the selected distribution
    double generate_value(DistributionType dist, double param1, double param2) {
        switch (dist) {
            case DistributionType::GAUSSIAN:
                return generate_gaussian(param1, param2);
            case DistributionType::EXPONENTIAL:
                return generate_exponential(param1);
            case DistributionType::POISSON:
                return generate_poisson(param1);
            default:
                return 0.0;
        }
    }

    // Helper function to generate a random probability between 0 and 1
    double random_prob() {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0.0, 1.0);
        return dis(gen);
    }
};

// Function to initialize the Z vector with random community assignments
vector<unsigned int> initialize_Z(int num_nodes, int num_clusters) {
    vector<unsigned int> Z(num_nodes);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, num_clusters - 1);

    for (int i = 0; i < num_nodes; ++i) {
        Z[i] = dis(gen);  // Assign a random cluster to node i
    }

    return Z;
}

double compute_cluster_accuracy(const std::vector<unsigned int>& Z_pred,
                                const std::vector<int>& Z_true) {
    if (Z_pred.size() != Z_true.size()) {
        std::cerr << "Error: Z_pred and Z_true must be the same length!" << std::endl;
        return 0.0;
    }

    int n = Z_pred.size();

    // Get unique labels
    std::set<int> true_label_set(Z_true.begin(), Z_true.end());
    std::set<int> pred_label_set(Z_pred.begin(), Z_pred.end());

    if (true_label_set.size() != pred_label_set.size()) {
        std::cerr << "Error: Number of clusters in Z_true and Z_pred must match!" << std::endl;
        return 0.0;
    }

    std::vector<int> true_labels(true_label_set.begin(), true_label_set.end());
    std::vector<int> pred_labels(pred_label_set.begin(), pred_label_set.end());

    std::vector<int> perm = pred_labels;
    double best_accuracy = 0.0;

    do {
        // Map predicted label to true label based on current permutation
        std::unordered_map<int, int> label_map;
        for (size_t i = 0; i < perm.size(); ++i) {
            label_map[perm[i]] = true_labels[i];
        }

        int correct = 0;
        for (int i = 0; i < n; ++i) {
            if (label_map[Z_pred[i]] == Z_true[i]) {
                ++correct;
            }
        }

        double accuracy = static_cast<double>(correct) / n;
        best_accuracy = std::max(best_accuracy, accuracy);

    } while (std::next_permutation(perm.begin(), perm.end()));

    return best_accuracy;
}


int main() {
    // 1-2: Initialize adjacency list (A) and attribute vector Y
    int n = 300; // number of nodes
    int k = 3;  // number of communities

    vector<node_t*> A(n, nullptr);
    vector<double> Y(n);
    vector<int> Z_true(n);

    // 3: Initialize initial community assignments (z)
    vector<unsigned int> z = initialize_Z(n, k);  // Example initial clustering

    // 4: Initialize necessary matrices and vectors for the bregman_clustering algorithm

    // Edge count matrix (K x K)
    vector<vector<unsigned int>> edge_count(k, vector<unsigned int>(k, 0));

    // Edge weight matrix (K x K)
    vector<vector<double>> edge_weight(k, vector<double>(k, 0.0));

    // Node count per community (K)
    vector<unsigned int> node_count(k, 0);

    // Node attribute sum per community (K)
    vector<double> node_attribute(k, 0.0);

    // Neighbor count matrix (N x K)
    vector<vector<unsigned int>> neighbor_count(n, vector<unsigned int>(k, 0));

    // Edge probability matrix (K x K)
    vector<vector<double>> p(k, vector<double>(k, 0.0));

    // Average edge weight matrix (K x K)
    vector<vector<double>> mu(k, vector<double>(k, 0.0));

    // Average node attribute per community (K)
    vector<double> nu(k, 0.0);

    //Distribution type
    DistributionType distribution = DistributionType::GAUSSIAN; // 1: Gaussian, 2: Poisson, 3: Exponential

    // Example parameters for weight generation and attribute distributions
    DistributionType attr_dist = DistributionType::GAUSSIAN;    // Node attributes distribution
    DistributionType edge_dist = DistributionType::GAUSSIAN;    // Edge probabilities distribution
    DistributionType weight_dist = DistributionType::GAUSSIAN;   // Weight distribution
    double attr_param1 = 0.0;     // Mean for Gaussian distribution (attributes)
    double attr_param2 = 1.0;     // Stddev for Gaussian distribution (attributes)
    double edge_param1 = 0.0;     // Mean for Gaussian distribution (edges)
    double edge_param2 = 1.0;    // Stddev for Gaussian distribution (edges)
    double weight_param1 = 0.0;   // Mean for Poisson distribution (weights)
    double weight_param2 = 1.0;   // Not used for Poisson distribution
    double min_ = 1.0;            // Minimum weight
    double max_ = 5.0;            // Maximum weight

    // Create the WSBM generator
    WSBMGenerator generator(n, k, attr_dist, edge_dist, weight_dist, attr_param1, attr_param2, edge_param1, edge_param2, weight_param1, weight_param2, min_, max_);

    // Generate the WSBM
    tie(A, Y, Z_true) = generator.generate_WSBM();

    // Print the generated graph
    print_graph(A);

    // 5: Call the bregman_clustering function
    vector<unsigned int> final_assignment = bregman_clustering(A, Y, z, edge_count, edge_weight, node_count, node_attribute, neighbor_count, p, mu, nu, distribution);

    // Output the final community assignment
    cout << "Final community assignment:" << endl;
    for (int i = 0; i < final_assignment.size(); ++i) {
        cout << "Node " << i << ": Community " << final_assignment[i] << endl;
    }

    cout << "Accuracy" << compute_cluster_accuracy(final_assignment, Z_true);

    // Cleanup dynamically allocated memory (if any)
    for (int i = 0; i < n; ++i) {
        node_t* current = A[i];
        while (current != nullptr) {
            node_t* to_delete = current;
            current = current->next;
            delete to_delete;
        }
    }

    return 0;
}
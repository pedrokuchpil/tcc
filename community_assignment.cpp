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
#include <chrono>

#include "wsbm_io_bin.hpp"   // node_t, load_wsbm_bin(...), free_graph(...)
#include "WSBMGenerator.hpp" // only for DistributionType in the clustering signature (if needed)

using namespace std;

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
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();      
    vector<unsigned int> z_new(z.size(), 0);  // Initialize z_new with zeros
    int max_iterations = 20;  // Maximum number of iterations
    int iterations = 0;

    bool converged = false;
    while (!converged) {
        update_parameters(A, Y, z, edge_count, edge_weight, node_count, node_attribute, neighbor_count, p, mu, nu);

        compute_maximum_assignment(A, Y, z, edge_count, edge_weight, node_count, node_attribute, neighbor_count, p, mu, nu, z_new, distribution);

        converged = (max_iterations <= iterations); // maximum number of iterations OR if there is no change OR number of changes is smaller than C we stop

        z = z_new;

        iterations++;
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Time difference = " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
    return z;
}
void add_edge(vector<node_t*>& A, int u, int v, double weight) {
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

double compute_cluster_accuracy(const vector<unsigned int>& Z_pred,
                                const vector<unsigned int>& Z_true) {
    if (Z_pred.size() != Z_true.size()) {
        cerr << "Error: Z_pred and Z_true must be the same length!" << endl;
        return 0.0;
    }

    int n = Z_pred.size();

    // Get unique labels
    set<int> true_label_set(Z_true.begin(), Z_true.end());
    set<int> pred_label_set(Z_pred.begin(), Z_pred.end());

    if (true_label_set.size() != pred_label_set.size()) {
        cerr << "Error: Number of clusters in Z_true and Z_pred must match!" << endl;
        cout << "Z_true" << true_label_set.size() << "Z_pred" << pred_label_set.size() << endl;
        return 0.0;
    }

    vector<int> true_labels(true_label_set.begin(), true_label_set.end());
    vector<int> pred_labels(pred_label_set.begin(), pred_label_set.end());

    vector<int> perm = pred_labels;
    double best_accuracy = 0.0;

    do {
        // Map predicted label to true label based on current permutation
        unordered_map<int, int> label_map;
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
        best_accuracy = max(best_accuracy, accuracy);

    } while (next_permutation(perm.begin(), perm.end()));

    return best_accuracy;
}

double log2_safe(double x) {
    return x > 0 ? log(x) / log(2) : 0;
}

double entropy(const unordered_map<int, int>& label_counts, int total) {
    double H = 0.0;
    for (const auto& [label, count] : label_counts) {
        double p = static_cast<double>(count) / total;
        H -= p * log2_safe(p);
    }
    return H;
}

double normalized_mutual_information(const vector<unsigned int>& Z_pred,
                                     const vector<unsigned int>& Z_true) {
    if (Z_pred.size() != Z_true.size() || Z_pred.empty())
        throw invalid_argument("Vectors must be of the same non-zero length");

    const int N = Z_pred.size();

    unordered_map<int, int> true_counts;
    unordered_map<int, int> pred_counts;
    unordered_map<int, unordered_map<int, int>> joint_counts;

    for (int i = 0; i < N; ++i) {
        int true_label = Z_true[i];
        int pred_label = Z_pred[i];
        true_counts[true_label]++;
        pred_counts[pred_label]++;
        joint_counts[true_label][pred_label]++;
    }

    double I = 0.0;
    for (const auto& [true_label, pred_map] : joint_counts) {
        for (const auto& [pred_label, joint_count] : pred_map) {
            double p_ij = static_cast<double>(joint_count) / N;
            double p_i = static_cast<double>(true_counts[true_label]) / N;
            double p_j = static_cast<double>(pred_counts[pred_label]) / N;
            I += p_ij * log2_safe(p_ij / (p_i * p_j));
        }
    }

    double H_true = entropy(true_counts, N);
    double H_pred = entropy(pred_counts, N);

    if (H_true == 0 || H_pred == 0)
        return 0.0;

    return I / sqrt(H_true * H_pred);
}

inline double comb2(int n) {
    return n < 2 ? 0.0 : (n * (n - 1)) / 2.0;
}

double adjusted_rand_index(const std::vector<unsigned int>& Z_pred,
                           const std::vector<unsigned int>& Z_true) {
    if (Z_pred.size() != Z_true.size() || Z_pred.empty())
        throw std::invalid_argument("Vectors must be of the same non-zero length");

    const int N = Z_pred.size();

    // Count true and predicted labels
    std::unordered_map<int, int> true_counts;
    std::unordered_map<int, int> pred_counts;
    std::unordered_map<int, std::unordered_map<int, int>> contingency;

    for (int i = 0; i < N; ++i) {
        int true_label = Z_true[i];
        int pred_label = Z_pred[i];
        true_counts[true_label]++;
        pred_counts[pred_label]++;
        contingency[true_label][pred_label]++;
    }

    double sum_comb_c = 0.0;
    for (const auto& [true_label, pred_map] : contingency) {
        for (const auto& [pred_label, count] : pred_map) {
            sum_comb_c += comb2(count);
        }
    }

    double sum_comb_true = 0.0;
    for (const auto& [label, count] : true_counts) {
        sum_comb_true += comb2(count);
    }

    double sum_comb_pred = 0.0;
    for (const auto& [label, count] : pred_counts) {
        sum_comb_pred += comb2(count);
    }

    double comb_n = comb2(N);
    if (comb_n == 0.0) return 1.0; // Perfect match for single elements

    double expected_index = (sum_comb_true * sum_comb_pred) / comb_n;
    double max_index = 0.5 * (sum_comb_true + sum_comb_pred);
    double ARI = (sum_comb_c - expected_index) / (max_index - expected_index);

    return ARI;
}

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " in_file n k\n";
    std::cerr << "Example: " << prog << " wsbm_data.bin 100 4\n";
}

int main(int argc, char** argv) {
    if (argc != 4) { usage(argv[0]); return 1; }

    std::string in_path = argv[1];
    int n = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);

    if (n <= 0 || k <= 0) {
        usage(argv[0]); return 1;
    }

    vector<node_t*> A(n, nullptr);
    vector<double> Y(n);
    vector<unsigned int> Z_true(n);

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

    try {
        load_wsbm_bin(in_path, A, Y, Z_true);
    } catch (const std::exception& e) {
        std::cerr << "Failed to load '" << in_path << "': " << e.what() << "\n";
        return 1;
    }


    // 5: Call the bregman_clustering function
    vector<unsigned int> final_assignment = bregman_clustering(A, Y, z, edge_count, edge_weight, node_count, node_attribute, neighbor_count, p, mu, nu, distribution);

    cout << "Accuracy: " << compute_cluster_accuracy(final_assignment, Z_true) << endl;
    cout << "NMI: " << normalized_mutual_information(final_assignment, Z_true) << endl;
    cout << "ARI: " << adjusted_rand_index(final_assignment, Z_true) << endl;


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
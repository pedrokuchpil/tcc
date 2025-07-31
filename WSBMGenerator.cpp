// WSBMGenerator.cpp
#include "WSBMGenerator.hpp"
#include <random>
#include <cmath>
#include <stdexcept>
using namespace std;

// File-local RNG shared by helper functions
namespace {
    std::mt19937& rng() {
        static thread_local std::mt19937 gen{ std::random_device{}() };
        return gen;
    }

    // File-local helper to add an edge to the linked-list adjacency.
    // Matches your node_t {id, weight, next}.
    inline void add_edge(std::vector<node_t*>& G, int u, int v, double w) {
        // Allocate if needed? For linked-list head, the vector holds head pointers (can be nullptr).
        node_t* e = new node_t{ v, w, G[u] };
        G[u] = e;
    }
}

// ---- Constructor ----
WSBMGenerator::WSBMGenerator(int n, int k,
                             DistributionType attr_dist, DistributionType edge_dist, DistributionType weight_dist,
                             double attr_p1, double attr_p2,
                             double edge_p1, double edge_p2,
                             double weight_p1, double weight_p2,
                             double minv, double maxv,
                             double p, double q)
    : attribute_dist(attr_dist), edge_dist(edge_dist), weight_dist(weight_dist),
      attr_param1(attr_p1), attr_param2(attr_p2),
      edge_param1(edge_p1), edge_param2(edge_p2),
      weight_param1(weight_p1), weight_param2(weight_p2),
      min_(minv), max_(maxv), pin(p), pout(q)
{
    n_clusters = k;
    N = n;
    if (k <= 0 || n <= 0) throw std::invalid_argument("WSBMGenerator: n and k must be > 0");
    nodes_per_cluster = n / k;
    if (nodes_per_cluster <= 0) throw std::invalid_argument("WSBMGenerator: n/k must be >= 1");

    generate_probability_matrix();
    generate_weight_centers();
}

// ---- Generate ----
std::tuple<std::vector<node_t*>, std::vector<double>, std::vector<unsigned int>>
WSBMGenerator::generate_WSBM(bool complete_graph)
{
    std::vector<node_t*> graph(N, nullptr);
    node_attributes.resize(N);
    std::vector<unsigned int> Z_true(N);

    // Node attributes and ground-truth labels
    for (int i = 0; i < N; ++i) {
        node_attributes[i] = generate_value(attribute_dist, attr_param1, attr_param2);
        Z_true[i] = static_cast<unsigned int>(i / nodes_per_cluster);
    }

    // Edges
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            int block_i = i / nodes_per_cluster;
            int block_j = j / nodes_per_cluster;

            double base_prob = complete_graph ? 1.0 : probability_matrix[block_i][block_j];
            double mod_prob  = base_prob * generate_value(edge_dist, edge_param1, edge_param2);

            if (random_prob() <= mod_prob) {
                double mean = weight_centers[block_i][block_j];
                double w = generate_value(weight_dist, weight_param1, weight_param2);

                // Undirected: add both directions
                add_edge(graph, i, j, w);
                add_edge(graph, j, i, w);
            }
        }
    }

    return std::make_tuple(std::move(graph), std::move(node_attributes), std::move(Z_true));
}

// ---- Diagnostics ----
void WSBMGenerator::print_node_attributes() {
    for (int i = 0; i < static_cast<int>(node_attributes.size()); ++i) {
        std::cout << "Node " << i << ": Attribute = " << node_attributes[i] << std::endl;
    }
}

// ---- Helpers ----
void WSBMGenerator::generate_probability_matrix() {
    probability_matrix.assign(n_clusters, std::vector<double>(n_clusters, pout));
    for (int i = 0; i < n_clusters; ++i) probability_matrix[i][i] = pin;
}

void WSBMGenerator::generate_weight_centers() {
    weight_centers.assign(n_clusters, std::vector<double>(n_clusters, 0.0));
    int index = 0;
    const int total = n_clusters * (n_clusters + 1) / 2;
    for (int i = 0; i < n_clusters; ++i) {
        for (int j = i; j < n_clusters; ++j) {
            double value = min_ + (max_ - min_) * (static_cast<double>(index) / total);
            weight_centers[i][j] = value;
            weight_centers[j][i] = value;
            ++index;
        }
    }
}

double WSBMGenerator::generate_value(DistributionType dist, double p1, double p2) {
    switch (dist) {
        case DistributionType::GAUSSIAN:     return generate_gaussian(p1, p2); // mean, stddev
        case DistributionType::EXPONENTIAL:  return generate_exponential(p1);   // lambda
        case DistributionType::POISSON:      return generate_poisson(p1);       // mean
        default: return 0.0;
    }
}

double WSBMGenerator::generate_gaussian(double mean, double stddev) {
    std::normal_distribution<double> d(mean, std::max(stddev, 1e-12));
    return d(rng());
}

double WSBMGenerator::generate_exponential(double lambda) {
    std::exponential_distribution<double> d(std::max(lambda, 1e-12));
    return d(rng());
}

double WSBMGenerator::generate_poisson(double mean) {
    std::poisson_distribution<int> d(std::max(mean, 1e-12));
    return static_cast<double>(d(rng()));
}

double WSBMGenerator::random_prob() {
    std::uniform_real_distribution<double> d(0.0, 1.0);
    return d(rng());
}

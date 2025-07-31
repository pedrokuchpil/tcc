// WSBMGenerator.hpp
#pragma once
#include <vector>
#include <tuple>
#include <iostream>

// Bring in node_t (linked-list) definition used across the project
#include "wsbm_io_bin.hpp"   // defines: struct node_t { int id; double weight; node_t* next; }

// ---- Distributions used by the generator ----
enum class DistributionType {
    GAUSSIAN,
    EXPONENTIAL,
    POISSON
};

class WSBMGenerator {
public:
    WSBMGenerator(int n, int k,
                  DistributionType attr_dist, DistributionType edge_dist, DistributionType weight_dist,
                  double attr_param1, double attr_param2,
                  double edge_param1, double edge_param2,
                  double weight_param1, double weight_param2,
                  double min_, double max_,
                  double p, double q);

    // Generate the WSBM: adjacency (linked lists), attributes, and ground-truth labels
    std::tuple<std::vector<node_t*>, std::vector<double>, std::vector<unsigned int>>
    generate_WSBM(bool complete_graph = false);

    void print_node_attributes();

private:
    // ---- Parameters / state ----
    int n_clusters;
    int nodes_per_cluster;
    int N;
    double pout, pin;
    double attr_param1, attr_param2;
    double edge_param1, edge_param2;
    double weight_param1, weight_param2;
    double min_, max_;

    std::vector<double> node_attributes;
    std::vector<std::vector<double>> probability_matrix;
    std::vector<std::vector<double>> weight_centers;
    DistributionType attribute_dist;
    DistributionType edge_dist;
    DistributionType weight_dist;

    // ---- Helpers (implemented in .cpp) ----
    void generate_probability_matrix();
    void generate_weight_centers();

    double generate_value(DistributionType dist, double param1, double param2);
    double generate_gaussian(double mean, double stddev);
    double generate_exponential(double lambda);
    double generate_poisson(double mean);

    double random_prob();
};

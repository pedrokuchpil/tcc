// wsbm_write.cpp
#include <iostream>
#include <tuple>
#include <vector>

#include "WSBMGenerator.hpp"  // DistributionType + WSBMGenerator
#include "wsbm_io_bin.hpp"    // save_wsbm_bin, free_graph

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " in_file n k pin pout\n";
    std::cerr << "Example: " << prog << " wsbm_data.bin 100 4 0.6 0.1\n";
}

int main(int argc, char** argv) {
    if (argc < 5 || argc > 6) { usage(argv[0]); return 1; }

    int n = std::atoi(argv[1]);
    int k = std::atoi(argv[2]);
    double pin = std::atof(argv[3]);
    double pout = std::atof(argv[4]);
    std::string out_path = (argc == 6) ? argv[5] : "wsbm_data.bin";

    DistributionType attr_dist   = DistributionType::GAUSSIAN;
    DistributionType edge_dist   = DistributionType::GAUSSIAN;
    DistributionType weight_dist = DistributionType::GAUSSIAN;

    double attr_param1   = 0.0, attr_param2    = 1.0; // mean, stddev for GAUSSIAN
    double edge_param1   = 1.0, edge_param2    = 0.2; // used by generate_value; edge_param2 unused for EXPONENTIAL/POISSON
    double weight_param1 = 1.0, weight_param2  = 0.3; // mean/stddev or lambda, depending on dist
    double min_w = 0.1, max_w = 2.0;

    try {
        WSBMGenerator gen(n, k,
                          attr_dist, edge_dist, weight_dist,
                          attr_param1, attr_param2,
                          edge_param1, edge_param2,
                          weight_param1, weight_param2,
                          min_w, max_w,
                          pin, pout);

        std::vector<node_t*> A;
        std::vector<double> Y;
        std::vector<unsigned int> Z_true;

        std::tie(A, Y, Z_true) = gen.generate_WSBM(false);

        save_wsbm_bin(out_path, A, Y, Z_true);
        std::cout << "Saved " << out_path
                  << " (N=" << A.size()
                  << ", |Y|=" << Y.size()
                  << ", |Z|=" << Z_true.size() << ")\n";

        free_graph(A); // avoid leaks
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

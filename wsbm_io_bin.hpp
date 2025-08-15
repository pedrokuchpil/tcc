#pragma once
#include <vector>
#include <tuple>
#include <fstream>
#include <stdexcept>
#include <cstdint>
#include <type_traits>

// ===== Graph node type (single source of truth) =====
struct node_t {
    int id;           // neighbor id
    double weight;    // edge weight
    node_t* next;     // next edge in adjacency list
};

// ====== Binary format (little-endian on most platforms) ======
// File layout:
// int N
// int M
// repeated M times: { int u, int v, double w }   // each undirected edge once, with u < v
// double Y[N]
// unsigned int Z_true[N]

// Compile-time sanity checks (optional, but helpful for portability)
static_assert(sizeof(int) == 4, "This binary format assumes 32-bit int.");
static_assert(sizeof(unsigned int) == 4, "This binary format assumes 32-bit unsigned int.");
static_assert(sizeof(double) == 8, "This binary format assumes 64-bit double.");

// ===== Save {A, Y, Z_true} to binary file =====
inline void save_wsbm_bin(
    const std::string& filename,
    const std::vector<node_t*>& A,
    const std::vector<double>& Y,
    const std::vector<unsigned int>& Z_true)
{
    const int N = static_cast<int>(A.size());
    if (static_cast<int>(Y.size()) != N || static_cast<int>(Z_true.size()) != N) {
        throw std::runtime_error("save_wsbm_bin: size mismatch among A, Y, Z_true");
    }

    // Collect undirected edges once (u < v) to avoid duplicates
    std::vector<std::tuple<int,int,double>> edges;
    edges.reserve(N * 4); // heuristic
    for (int u = 0; u < N; ++u) {
        for (const node_t* e = A[u]; e; e = e->next) {
            const int v = e->id;
            if (u < v) {
                edges.emplace_back(u, v, e->weight);
            }
        }
    }
    const int M = static_cast<int>(edges.size());

    std::ofstream out(filename, std::ios::binary);
    if (!out) throw std::runtime_error("save_wsbm_bin: cannot open file for writing: " + filename);

    // Header
    out.write(reinterpret_cast<const char*>(&N), sizeof(int));
    out.write(reinterpret_cast<const char*>(&M), sizeof(int));

    // Edges
    for (const auto& t : edges) {
        const int u = std::get<0>(t);
        const int v = std::get<1>(t);
        const double w = std::get<2>(t);
        out.write(reinterpret_cast<const char*>(&u), sizeof(int));
        out.write(reinterpret_cast<const char*>(&v), sizeof(int));
        out.write(reinterpret_cast<const char*>(&w), sizeof(double));
    }

    // Y
    out.write(reinterpret_cast<const char*>(Y.data()), static_cast<std::streamsize>(N * sizeof(double)));

    // Z_true
    out.write(reinterpret_cast<const char*>(Z_true.data()), static_cast<std::streamsize>(N * sizeof(unsigned int)));

    if (!out.good()) throw std::runtime_error("save_wsbm_bin: write error");
}

inline void save_wsbm_txt(
    const std::string& filename,
    const std::vector<node_t*>& A,
    const std::vector<double>& Y,
    const std::vector<unsigned int>& Z_true)
{
    const int N = static_cast<int>(A.size());
    if (static_cast<int>(Y.size()) != N || static_cast<int>(Z_true.size()) != N) {
        throw std::runtime_error("save_wsbm_txt: size mismatch among A, Y, Z_true");
    }

    // Collect undirected edges once (u < v)
    std::vector<std::tuple<int,int,double>> edges;
    edges.reserve(N * 4); // heuristic
    for (int u = 0; u < N; ++u) {
        for (const node_t* e = A[u]; e; e = e->next) {
            const int v = e->id;
            if (u < v) {
                edges.emplace_back(u, v, e->weight);
            }
        }
    }
    const int M = static_cast<int>(edges.size());

    std::ofstream out(filename);
    if (!out) throw std::runtime_error("save_wsbm_txt: cannot open file for writing: " + filename);

    // Write header
    out << N << " " << M << "\n";

    // Write edges
    for (const auto& t : edges) {
        out << std::get<0>(t) << " "
            << std::get<1>(t) << " "
            << std::get<2>(t) << "\n";
    }

    // Write Y
    out << "# Y\n";
    for (double y : Y) {
        out << y << "\n";
    }

    // Write Z_true
    out << "# Z_true\n";
    for (unsigned int z : Z_true) {
        out << z << "\n";
    }

    if (!out.good()) throw std::runtime_error("save_wsbm_txt: write error");
}


// ===== Load {A, Y, Z_true} from binary file =====
inline void load_wsbm_bin(
    const std::string& filename,
    std::vector<node_t*>& A,
    std::vector<double>& Y,
    std::vector<unsigned int>& Z_true)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in) throw std::runtime_error("load_wsbm_bin: cannot open file for reading: " + filename);

    int N = 0, M = 0;
    in.read(reinterpret_cast<char*>(&N), sizeof(int));
    in.read(reinterpret_cast<char*>(&M), sizeof(int));
    if (!in.good() || N < 0 || M < 0) {
        throw std::runtime_error("load_wsbm_bin: bad header");
    }

    A.assign(N, nullptr);

    auto add_edge_one_dir = [&](int src, int dst, double w) {
        node_t* e = new node_t{dst, w, A[src]};
        A[src] = e;
    };

    // Edges (undirected)
    for (int i = 0; i < M; ++i) {
        int u, v; double w;
        in.read(reinterpret_cast<char*>(&u), sizeof(int));
        in.read(reinterpret_cast<char*>(&v), sizeof(int));
        in.read(reinterpret_cast<char*>(&w), sizeof(double));
        if (!in.good() || u < 0 || v < 0 || u >= N || v >= N) {
            throw std::runtime_error("load_wsbm_bin: bad edge at index " + std::to_string(i));
        }
        add_edge_one_dir(u, v, w);
        add_edge_one_dir(v, u, w);
    }

    // Y
    Y.resize(N);
    in.read(reinterpret_cast<char*>(Y.data()), static_cast<std::streamsize>(N * sizeof(double)));
    if (!in.good()) throw std::runtime_error("load_wsbm_bin: failed reading Y");

    // Z_true
    Z_true.resize(N);
    in.read(reinterpret_cast<char*>(Z_true.data()), static_cast<std::streamsize>(N * sizeof(unsigned int)));
    if (!in.good()) throw std::runtime_error("load_wsbm_bin: failed reading Z_true");
}

// ===== Utility: free adjacency lists (use in your programs to avoid leaks) =====
inline void free_graph(std::vector<node_t*>& A) {
    for (node_t*& head : A) {
        while (head) {
            node_t* tmp = head;
            head = head->next;
            delete tmp;
        }
        head = nullptr;
    }
}
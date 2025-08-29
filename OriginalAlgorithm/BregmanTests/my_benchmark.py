import struct
import numpy as np
import networkx as nx
from distributions import *
from BregmanClustering.models import BregmanNodeEdgeAttributeGraphClusteringSoft as softBreg
from BregmanKernel.kernel_models import BregmanKernelClustering
from scipy.sparse import csr_matrix, csc_matrix, csr_array
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, confusion_matrix, normalized_mutual_info_score
from sklearn.datasets import make_classification
import time
import pandas as pd
import subprocess
from tqdm import tqdm
import os
from copy import deepcopy
import csv
import argparse

def read_wsbm_txt(filename):
    """
    Reads a WSBM graph saved in human-readable text format:
      - header: N M
      - M edges: u v w
      - Y (node attributes)
      - Z_true (true labels)
    Returns: A (NxN), Y (Nx1), Z_true (N,)
    """
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        lines = [line.strip() for line in f if line.strip()]

    # --- Header ---
    N, M = map(int, lines[0].split())

    # --- Adjacency ---
    A = np.zeros((N, N))
    for i in range(1, M + 1):
        u, v, w = lines[i].split()
        u, v = int(u), int(v)
        w = float(w)
        if u >= N or v >= N:
            raise ValueError(f"Edge ({u},{v}) exceeds N={N}")
        A[u, v] = w
        A[v, u] = w  # undirected

    # --- Find Y and Z_true sections ---
    try:
        y_idx = lines.index("# Y") + 1
        z_idx = lines.index("# Z_true")
    except ValueError:
        raise ValueError("Missing '# Y' or '# Z_true' sections in file")

    # --- Node attributes ---
    Y_lines = lines[y_idx:z_idx]
    if len(Y_lines) != N:
        raise ValueError(f"Number of Y entries ({len(Y_lines)}) does not match N ({N})")
    Y = np.array([float(y) for y in Y_lines], dtype=np.float64).reshape(N, 1)

    # --- True labels ---
    Z_lines = lines[z_idx + 1 :]
    if len(Z_lines) != N:
        raise ValueError(f"Number of Z_true entries ({len(Z_lines)}) does not match N ({N})")
    Z_true = np.array([int(z) for z in Z_lines], dtype=np.uint32)

    return A, Y, Z_true

def prepare_matrices(filename):
    """
    Reads a WSBM file and returns matrices aligned by N nodes.
    Returns:
        X: adjacency matrix (N,N)
        Y: node attributes (N,1)
        Z_true: true labels (N,)
    """
    X, Y, Z_true = read_wsbm_txt(filename)
    N = X.shape[0]

    # Sanity checks
    if Y.shape[0] != N or Z_true.shape[0] != N:
        raise ValueError("Mismatch in number of nodes between X, Y, and Z_true")

    return X, Y, Z_true

def clustering_accuracy(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    r, c = linear_sum_assignment(cm.max() - cm)
    return cm[r, c].sum() / cm.sum()

def run_from_wsbm_file(wsbm_path, runs=10, csv_path="python_benchmark_results.csv", n_iters=10, **model_kwargs):
    X, Y, z_true = prepare_matrices(wsbm_path)
    print("X:", X.shape)
    print("Y:", Y.shape)
    print("Z_true:", z_true.shape)

    N = int(X.shape[0])
    K = int(len(np.unique(z_true)))

    edge_index = np.array(X.nonzero())
    E = X[edge_index[0], edge_index[1]][:, None]

    needs_header = not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0
    with open(csv_path, "a", newline="") as f:
        w = csv.writer(f)
        if needs_header:
            w.writerow(["Run", "N", "K", "Time(ms)", "Accuracy", "NMI", "ARI"])

        for run_id in range(1, runs + 1):
            model = softBreg(n_clusters=K, n_iters=n_iters, **model_kwargs)

            t0 = time.perf_counter()
            model.initialize(edge_index, E, Y)
            model.assignInitialLabels(None, None)
            z_init = np.zeros(N).reshape(N, 1)
            for i in range(0,N,2):
                z_init[i] = 1
            z_pred = model.fit(edge_index, E, Y).predict(X, Y)
            dt_ms = (time.perf_counter() - t0) * 1000.0

            acc = clustering_accuracy(z_true, z_pred)
            nmi = normalized_mutual_info_score(z_true, z_pred)
            ari = adjusted_rand_score(z_true, z_pred)

            w.writerow([run_id, N, K, f"{dt_ms:.3f}", f"{acc:.6f}", f"{nmi:.6f}", f"{ari:.6f}"])

    print(f"✅ Done. Appended {runs} runs to {csv_path} for N={N}, K={K}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run WSBM model from a file and record results.")
    parser.add_argument("--wsbm_path", required=True, help="Path to the WSBM binary file.")
    parser.add_argument("--runs", type=int, default=10, help="Number of runs.")
    parser.add_argument("--csv_path", required=True, help="Path to CSV output file.")
    parser.add_argument("--n_iters", type=int, default=1, help="Number of iterations for the model.")

    args = parser.parse_args()

    run_from_wsbm_file(
        wsbm_path=args.wsbm_path,
        runs=args.runs,
        csv_path=args.csv_path,
        n_iters=args.n_iters,
        attributeDistribution="gaussian",
        edgeDistribution="gaussian",
        weightDistribution="gaussian",
        reduce_by=None,
        divergence_precomputed=False,
        initializer="random",
        use_random_init=True,
    )

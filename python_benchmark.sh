#!/bin/bash
# run_full_benchmark.sh
# Generates WSBM graphs and runs Python benchmark, storing results in CSV

set -e  # exit on error

# Paths
ROOT_DIR=$(pwd)
WSBM_WRITE="$ROOT_DIR/wsbm_write"
PY_DIR="$ROOT_DIR/OriginalWork/BregmanTests"
PY_RUN="$PY_DIR/my_benchmark.py"

# CSV file
CSV_FILE="$ROOT_DIR/benchmark_results.csv"
rm -f "$CSV_FILE"

# Parameters
P_IN=0.8
P_OUT=0.2
RUNS=10
N_LIST=(100 500)
K_LIST=(2 4)

# Loop over N and K
for N in "${N_LIST[@]}"; do
    for K in "${K_LIST[@]}"; do
        echo "=== Generating WSBM for N=$N, K=$K ==="
        $WSBM_WRITE "$N" "$K" "$P_IN" "$P_OUT" "$ROOT_DIR/wsbm_data.bin" "$ROOT_DIR/wsbm_data.txt" 

        echo "=== Running Python benchmark for N=$N, K=$K ==="
        PYTHONPATH="$ROOT_DIR/OriginalWork" python3 "$PY_RUN" --wsbm_path "$ROOT_DIR/wsbm_data.txt" --runs "$RUNS" --csv_path "$CSV_FILE"
    done
done

echo "✅ All benchmarks complete. Results saved to $CSV_FILE."

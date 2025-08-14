#!/bin/bash

WSBM_EXEC="./wsbm_write"
MAIN_EXEC="./main"
OUTPUT_FILE="augmenting_k_results.csv"
BIN_FILE="wsbm_data.bin"
RUNS=10

# Parameter values
N_values=(10000)
K_values=(2 4 8)
P_values=(0.02)


# Write CSV header
echo "Run,N,K,Time(ms),Accuracy,NMI,ARI,P" > "$OUTPUT_FILE"

# Loop over N and K
for N in "${N_values[@]}"; do
  for K in "${K_values[@]}"; do
    for P in "${P_values[@]}"; do
        echo "Generating WSBM with N=$N, K=$K..."

        # Generate new binary file with default pin/pout
        $WSBM_EXEC "$N" "$K" "$P" 0.004 "$BIN_FILE"

        # Check if generation succeeded
        if [ $? -ne 0 ]; then
        echo "WSBM generation failed for N=$N, K=$K"
        continue
        fi

        # Run main 10 times and log results
        for i in $(seq 1 $RUNS); do
        echo "Run $i for N=$N, K=$K"

        OUTPUT=$($MAIN_EXEC "$BIN_FILE" "$N" "$K")

        TIME=$(echo "$OUTPUT" | grep "Time difference" | awk -F= '{print $2}' | sed 's/\[ms\]//g' | xargs)
        ACC=$(echo "$OUTPUT" | grep "Accuracy:" | awk '{print $2}')
        NMI=$(echo "$OUTPUT" | grep "NMI:" | awk '{print $2}')
        ARI=$(echo "$OUTPUT" | grep "ARI:" | awk '{print $2}')
        PIN=$P

        echo "$i,$N,$K,$TIME,$ACC,$NMI,$ARI,$PIN" >> "$OUTPUT_FILE"
        done
    done
  done
done

echo "Benchmark complete. Results saved to $OUTPUT_FILE"

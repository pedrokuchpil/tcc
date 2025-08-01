import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load CSV
df = pd.read_csv("benchmark_results.csv")

# Convert columns to proper numeric types (handle errors silently)
df['Time(ms)'] = pd.to_numeric(df['Time(ms)'], errors='coerce')
df['Accuracy'] = pd.to_numeric(df['Accuracy'], errors='coerce')
df['NMI'] = pd.to_numeric(df['NMI'], errors='coerce')
df['N'] = pd.to_numeric(df['N'], errors='coerce')
df['K'] = pd.to_numeric(df['K'], errors='coerce')

# Drop rows with missing (NaN) values in any of the relevant columns
df = df.dropna(subset=['N', 'K', 'Time(ms)', 'Accuracy', 'NMI'])

# Group and average
agg_df = df.groupby(['K', 'N']).agg({
    'Time(ms)': 'mean',
    'Accuracy': 'mean',
    'NMI': 'mean'
}).reset_index()

# Plot settings
sns.set(style="whitegrid")
k_values = sorted(agg_df['K'].unique())

# Create plots
for k in k_values:
    data_k = agg_df[agg_df['K'] == k]

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot Time on left Y-axis
    color = 'tab:red'
    ax1.set_xlabel("Number of Nodes (N)")
    ax1.set_ylabel("Time (ms)", color=color)
    ax1.plot(data_k['N'], data_k['Time(ms)'], marker='o', color=color, label='Time (ms)')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_xscale('log')

    # Create second y-axis for Accuracy and NMI
    ax2 = ax1.twinx()  # share same x-axis
    color_acc = 'tab:blue'
    color_nmi = 'tab:green'
    ax2.set_ylabel("Accuracy / NMI")
    ax2.plot(data_k['N'], data_k['Accuracy'], marker='s', color=color_acc, label='Accuracy')
    ax2.plot(data_k['N'], data_k['NMI'], marker='^', color=color_nmi, label='NMI')
    ax2.tick_params(axis='y')

    # Title and layout
    plt.title(f"Benchmark Metrics for K = {k}")
    fig.tight_layout()

    # Create combined legend
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax2.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left')

    # Save and close
    plt.savefig(f"benchmark_k{k}.png")
    plt.close()

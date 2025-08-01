import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load CSV
df = pd.read_csv("augmenting_p_results.csv")

# Ensure correct types
df['P'] = pd.to_numeric(df['P'], errors='coerce')
df['Time(ms)'] = pd.to_numeric(df['Time(ms)'], errors='coerce')
df['Accuracy'] = pd.to_numeric(df['Accuracy'], errors='coerce')
df['NMI'] = pd.to_numeric(df['NMI'], errors='coerce')

# Drop invalid rows
df = df.dropna(subset=['P', 'Time(ms)', 'Accuracy', 'NMI'])

# Aggregate by P
agg_df = df.groupby('P').agg({
    'Time(ms)': 'mean',
    'Accuracy': 'mean',
    'NMI': 'mean'
}).reset_index()

# Begin plot
sns.set(style="whitegrid")
fig, ax1 = plt.subplots(figsize=(10, 6))

# Left Y-axis: Execution Time
color_time = 'tab:red'
ax1.set_xlabel("P")
ax1.set_ylabel("Execution Time (ms)", color=color_time)
ax1.plot(agg_df['P'], agg_df['Time(ms)'], marker='o', color=color_time, label='Time (ms)')
ax1.tick_params(axis='y', labelcolor=color_time)

# Right Y-axis: Accuracy and NMI
ax2 = ax1.twinx()
color_acc = 'tab:blue'
color_nmi = 'tab:green'
ax2.set_ylabel("Accuracy / NMI")
ax2.plot(agg_df['P'], agg_df['Accuracy'], marker='s', color=color_acc, label='Accuracy')
ax2.plot(agg_df['P'], agg_df['NMI'], marker='^', color=color_nmi, label='NMI')
ax2.tick_params(axis='y')

# Combine legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc='center right')

# Title and layout
plt.title("Execution Time, Accuracy, and NMI vs P (N=1000, K=4, Q=0.05)")
fig.tight_layout()
plt.savefig("time_accuracy_nmi_vs_p.png")
plt.close()

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import io
import re

csv_data = """Configuration,Total Cores,Computation (s),Communication (s),Total Time (s)
"1 Node, 1 Task, 8 Threads",8,37.33,0.00,37.33
"1 Node, 1 Task, 32 Threads",32,18.22,0.00,18.22
"1 Node, 2 Tasks, 16 Threads",32,11.88,0.45,12.33
"1 Node, 4 Tasks, 8 Threads",32,9.23,0.31,9.54
"2 Nodes, 2 Tasks, 1 Thread",2,88.79,0.33,89.12
"2 Nodes, 2 Tasks, 8 Threads",16,16.79,0.73,17.52
"2 Nodes, 2 Tasks, 16 Threads",32,11.58,0.43,12.01
"2 Nodes, 2 Tasks, 32 Threads",64,9.14,0.34,9.48
"2 Nodes, 4 Tasks, 16 Threads",64,5.84,0.44,6.28
"2 Nodes, 8 Tasks, 8 Threads",64,4.49,0.50,4.99
"4 Nodes, 4 Tasks, 16 Threads",64,5.84,0.32,6.16
"4 Nodes, 4 Tasks, 32 Threads",128,5.01,0.57,5.58
"4 Nodes, 16 Tasks, 8 Threads",128,2.36,0.40,2.76
"4 Nodes, 64 Tasks, 2 Threads",128,1.99,0.69,2.68
"8 Nodes, 8 Tasks, 32 Threads",256,2.75,0.73,3.48
"8 Nodes, 32 Tasks, 8 Threads",256,1.23,2.64,3.87
"8 Nodes, 128 Tasks, 2 Threads",256,1.22,0.92,2.14
"16 Nodes, 16 Tasks, 32 Threads",512,1.63,0.97,2.60
"16 Nodes, 64 Tasks, 8 Threads",512,0.64,3.04,3.68
"16 Nodes, 128 Tasks, 4 Threads",512,0.59,3.62,4.21
"""

# Read data
df = pd.read_csv(io.StringIO(csv_data))

# Function to create short labels and extract Node count for sorting
def parse_config(row):
    conf = row['Configuration']
    # Regex to find N, T (Tasks/Ranks), Th (Threads)
    # Pattern: "X Node(s), Y Task(s), Z Thread(s)"
    match = re.search(r'(\d+)\s+Nodes?,\s+(\d+)\s+Tasks?,\s+(\d+)\s+Threads?', conf)
    if match:
        n = int(match.group(1))
        r = int(match.group(2))
        t = int(match.group(3))
        label = f"{n}N {r}R {t}T"
        return pd.Series([n, r, t, label])
    return pd.Series([0, 0, 0, conf])

df[['Nodes', 'Ranks', 'Threads', 'ShortLabel']] = df.apply(parse_config, axis=1)

# Sort by Nodes, then by Total Time (descending) to show improvement? 
# Or by Ranks? 
# User asked "In increasing order of nodes only".
# Inside nodes, I'll keep the order provided in the CSV which seemed logical (often increasing cores/performance).
# Actually, let's explicit sort: Nodes (asc), Total Cores (asc), Total Time (desc).
df = df.sort_values(by=['Nodes', 'Total Cores', 'Total Time (s)'], ascending=[True, True, False])

# Set index for plotting
df_plot = df.set_index('ShortLabel')

# Plot
fig, ax = plt.subplots(figsize=(15, 8))
colors = ['#4c72b0', '#c44e52'] # Blue, Red

df_plot[['Computation (s)', 'Communication (s)']].plot(
    kind='bar', stacked=True, color=colors, ax=ax, width=0.7, logy=True
)

ax.set_title('Combined Performance Summary: 1 to 16 Nodes\nComputation vs Communication', fontsize=16)
ax.set_ylabel('Time (Seconds) - Log Scale', fontsize=14)
ax.set_xlabel('Configuration (N=Nodes, R=MPI Ranks, T=OpenMP Threads)', fontsize=14)
plt.xticks(rotation=45, ha='right')

# Formatting Y-axis to scalar
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())

# Add labels
for i, v in enumerate(df_plot['Total Time (s)']):
    # Add * for the last one if it was the unstable one (1.17s)
    label = f"{v:.2f}s"
    if v == 1.17:
        label += "*"
        # Add red text for failure
        ax.text(i, v * 1.1, label, ha='center', va='bottom', fontsize=10, fontweight='bold', color='red')
    else:
        ax.text(i, v * 1.1, label, ha='center', va='bottom', fontsize=10, fontweight='bold')

# Grid
ax.grid(axis='y', which='major', linestyle='-', alpha=0.5)
ax.grid(axis='y', which='minor', linestyle=':', alpha=0.3)

# Add vertical lines to separate Node groups
# Find indices where Node count changes
node_counts = df_plot['Nodes'].values
for i in range(len(node_counts) - 1):
    if node_counts[i] != node_counts[i+1]:
        ax.axvline(x=i + 0.5, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)

plt.tight_layout()
plt.savefig('combined_summary_plot.png')

print("Plot created.")
print(df[['ShortLabel', 'Total Cores', 'Computation (s)', 'Communication (s)', 'Total Time (s)']].to_markdown(index=False))
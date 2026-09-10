import pandas as pd
import matplotlib.pyplot as plt
import io

data_csv = """Configuration,Total Cores,Computation (s),Communication (s),Total Time (s)
"1 Node, 1 Task, 8 Threads",8,37.33,0.00,37.33
"1 Node, 1 Task, 32 Threads",32,18.22,0.00,18.22
"1 Node, 2 Tasks, 16 Threads",32,11.88,0.45,12.33
"1 Node, 4 Tasks, 8 Threads",32,9.23,0.31,9.54"""

df = pd.read_csv(io.StringIO(data_csv))

# Setup Plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot Stacked Bar
# We need 'Configuration' as index for x-axis labels
df.set_index('Configuration', inplace=True)

# Select only the components to stack
plot_data = df[['Computation (s)', 'Communication (s)']]

# Plot
plot_data.plot(kind='bar', stacked=True, ax=ax, color=['#4c72b0', '#c44e52'], width=0.6)

# Formatting
ax.set_title('Execution Time Breakdown: Computation vs Communication', fontsize=14)
ax.set_ylabel('Time (Seconds)', fontsize=12)
ax.set_xlabel('Configuration', fontsize=12)
plt.xticks(rotation=15, ha='right')

# Add Labels on top of bars
totals = df['Total Time (s)']
for i, v in enumerate(totals):
    ax.text(i, v + 0.5, f"{v:.2f}s", ha='center', va='bottom', fontsize=10, fontweight='bold')

# Add grid
ax.grid(axis='y', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('stacked_histogram_summary.png')
print("Plot saved.")
#!/usr/bin/env python3
"""
Script to visualize multi-GPU scaling results from atmospheric model logs.
Creates a stacked bar chart showing timing breakdown by number of GPUs.
"""

import re
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def parse_log_file(filepath):
    """Parse a log file and extract timing statistics and GPU count."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Extract number of MPI tasks (= number of GPUs)
    mpi_match = re.search(r'Number of MPI tasks:\s+(\d+)', content)
    if not mpi_match:
        return None
    
    n_gpus = int(mpi_match.group(1))
    
    # Extract number of nodes
    nodes_match = re.search(r'Nodes = (\d+)', content)
    n_nodes = int(nodes_match.group(1)) if nodes_match else 1
    
    # Extract problem size
    nx_match = re.search(r'nx_global\s*:\s*(\d+)', content)
    nz_match = re.search(r'nz\s*:\s*(\d+)', content)
    nx_global = int(nx_match.group(1)) if nx_match else None
    nz = int(nz_match.group(1)) if nz_match else None
    
    # Extract timing statistics
    timing_section = re.search(
        r'PARALLEL TIMING STATISTICS.*?-{100}(.*?)-{100}',
        content,
        re.DOTALL
    )
    
    if not timing_section:
        return None
    
    computation_time = 0.0
    communication_time = 0.0
    
    timing_lines = timing_section.group(1).strip().split('\n')
    
    for line in timing_lines:
        # Skip header and separator lines
        if 'Function' in line or line.strip().startswith('*') or not line.strip():
            continue
            
        # Parse timing lines: Function name, Max Total, Max Excl, Avg Total, Calls
        # Match lines with timing data
        match = re.match(r'\s*(.+?)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s*$', line)
        if match:
            func_name = match.group(1).strip()
            max_excl = float(match.group(3))  # Max Excl in microseconds
            
            # Categorize as computation or communication
            if 'Communication' in func_name or 'MPI' in func_name:
                communication_time += max_excl
            else:
                computation_time += max_excl
    
    return {
        'n_gpus': n_gpus,
        'n_nodes': n_nodes,
        'computation': computation_time,
        'communication': communication_time,
        'nx_global': nx_global,
        'nz': nz
    }

def main():
    # Find all log files
    log_dir = Path('/home/raiong/Downloads/logs_leo/logs')
    log_files = sorted(log_dir.glob('multigpu_scaling_*.out'))
    
    # Parse all log files
    results = []
    for log_file in log_files:
        data = parse_log_file(log_file)
        if data:
            results.append(data)
    
    # Sort by number of GPUs
    results.sort(key=lambda x: x['n_gpus'])
    
    if not results:
        print("No valid log files found!")
        return
    
    # Extract data for plotting
    n_gpus_list = [r['n_gpus'] for r in results]
    n_nodes_list = [r['n_nodes'] for r in results]
    computation_times = np.array([r['computation'] / 1_000_000 for r in results])  # Convert to seconds
    communication_times = np.array([r['communication'] / 1_000_000 for r in results])  # Convert to seconds
    
    # Create labels with node and GPU info
    labels = [f"{nodes}N×{gpus//nodes}G\n({gpus} GPUs)" for nodes, gpus in zip(n_nodes_list, n_gpus_list)]
    
    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(12, 8))
    
    x = np.arange(len(n_gpus_list))
    width = 0.6
    
    # Colors
    comp_color = '#2E86AB'  # Blue for computation
    comm_color = '#A23B72'  # Purple for communication
    
    # Create stacked bars
    p1 = ax.bar(x, computation_times, width, label='Computation', color=comp_color)
    p2 = ax.bar(x, communication_times, width, bottom=computation_times, label='Communication', color=comm_color)
    
    # Get problem size from first result
    nx_global = results[0].get('nx_global')
    nz = results[0].get('nz')
    problem_size = f" (nx={nx_global}, nz={nz})" if nx_global and nz else ""
    
    # Customize plot
    ax.set_xlabel('Server Configuration', fontsize=14, fontweight='bold')
    ax.set_ylabel('Time (seconds)', fontsize=14, fontweight='bold')
    ax.set_title(f'Multi-GPU Scaling: Computation vs Communication Time{problem_size}\n(Max Exclusive Times)', 
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.legend(loc='upper right', fontsize=12, framealpha=0.95)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Determine appropriate y-axis scale
    max_time = max(computation_times + communication_times)
    if max_time < 100:
        ax.set_ylim(0, max_time * 1.15)
    else:
        ax.set_ylim(0, max_time * 1.1)
    
    # Add total time labels on top of bars
    totals = computation_times + communication_times
    for i, (gpu_count, total, comp, comm) in enumerate(zip(n_gpus_list, totals, computation_times, communication_times)):
        # Total time on top
        ax.text(i, total + max_time * 0.02, f'{total:.1f}s', 
                ha='center', va='bottom', fontweight='bold', fontsize=10)
        
        # Percentage labels inside bars
        if comp > max_time * 0.05:  # Only show if bar is large enough
            ax.text(i, comp / 2, f'{(comp/total)*100:.0f}%', 
                    ha='center', va='center', color='white', fontweight='bold', fontsize=9)
        if comm > max_time * 0.05:
            ax.text(i, comp + comm / 2, f'{(comm/total)*100:.0f}%', 
                    ha='center', va='center', color='white', fontweight='bold', fontsize=9)
    
    plt.tight_layout()
    
    # Save figure
    output_file = log_dir / 'scaling_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    
    # Print summary statistics
    print("\n" + "="*80)
    print("SCALING SUMMARY - Max Exclusive Times")
    print("="*80)
    print(f"{'Config':<15} {'GPUs':<8} {'Computation (s)':<18} {'Communication (s)':<18} {'Total (s)':<12}")
    print("-"*80)
    
    for nodes, gpus, comp, comm, total in zip(n_nodes_list, n_gpus_list, computation_times, communication_times, totals):
        config = f"{nodes}N×{gpus//nodes}G"
        print(f"{config:<15} {gpus:<8} {comp:<18.2f} {comm:<18.2f} {total:<12.2f}")
    
    print("="*80)
    
    # Speedup and efficiency
    baseline_time = totals[0]
    print(f"\n{'Config':<15} {'Speedup':<12} {'Efficiency (%)':<15} {'Comm Overhead (%)':<20}")
    print("-"*80)
    for nodes, gpus, total, comm in zip(n_nodes_list, n_gpus_list, totals, communication_times):
        config = f"{nodes}N×{gpus//nodes}G"
        speedup = baseline_time / total
        efficiency = (speedup / gpus) * 100
        comm_overhead = (comm / total) * 100
        print(f"{config:<15} {speedup:<12.2f} {efficiency:<15.1f} {comm_overhead:<20.1f}")
    print("="*80)
    
    plt.show()

if __name__ == '__main__':
    main()

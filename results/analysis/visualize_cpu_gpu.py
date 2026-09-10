#!/usr/bin/env python3
"""
Script to compare CPU vs GPU performance from atmospheric model logs.
Creates visualizations comparing computation and communication times.
"""

import re
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def parse_log_file(filepath):
    """Parse a log file and extract timing statistics."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Extract execution mode
    is_gpu = 'OpenACC: ENABLED' in content
    is_cpu = 'OpenACC: DISABLED' in content
    
    # Extract number of MPI tasks
    mpi_match = re.search(r'Number of MPI tasks:\s+(\d+)', content)
    if not mpi_match:
        return None
    
    n_tasks = int(mpi_match.group(1))
    
    # Extract OpenMP threads
    omp_match = re.search(r'Number of OpenMP threads:\s+(\d+)', content)
    n_threads = int(omp_match.group(1)) if omp_match else 1
    
    # Extract number of nodes
    nodes_match = re.search(r'Nodes = (\d+)', content)
    n_nodes = int(nodes_match.group(1)) if nodes_match else 1
    
    # Extract problem size
    nx_match = re.search(r'nx_global\s*:\s*(\d+)', content)
    nz_match = re.search(r'nz\s*:\s*(\d+)', content)
    nx_global = int(nx_match.group(1)) if nx_match else None
    nz = int(nz_match.group(1)) if nz_match else None
    
    # Extract simulation time
    sim_time_match = re.search(r'final time\s*:\s*([\d.]+)', content)
    sim_time = float(sim_time_match.group(1)) if sim_time_match else None
    
    # Extract wall clock time
    wall_start = re.search(r'Wall Clock Start:\s+([\d:\.]+)', content)
    wall_end = re.search(r'Wall Clock End:\s+([\d:\.]+)', content)
    
    # Extract CPU time used
    cpu_time_match = re.search(r'USED CPU TIME:\s+([\d.]+)', content)
    total_cpu_time = float(cpu_time_match.group(1)) if cpu_time_match else None
    
    # Extract timing statistics
    timing_section = re.search(
        r'PARALLEL TIMING STATISTICS.*?-{100}(.*?)-{100}',
        content,
        re.DOTALL
    )
    
    if not timing_section:
        return None
    
    timings = {
        'thermal': 0.0,
        'hydrostatic': 0.0,
        'rungekutta': 0.0,
        'step': 0.0,
        'communication': 0.0,
        'mass_energy': 0.0,
        'init': 0.0
    }
    
    timing_lines = timing_section.group(1).strip().split('\n')
    
    for line in timing_lines:
        if 'Function' in line or line.strip().startswith('*') or not line.strip():
            continue
            
        match = re.match(r'\s*(.+?)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s*$', line)
        if match:
            func_name = match.group(1).strip()
            max_excl = float(match.group(3)) / 1_000_000  # Convert to seconds
            
            if 'thermal' in func_name.lower():
                timings['thermal'] = max_excl
            elif 'hydrostatic' in func_name.lower():
                timings['hydrostatic'] = max_excl
            elif 'rungekutta' in func_name.lower():
                timings['rungekutta'] = max_excl
            elif 'step' in func_name.lower():
                timings['step'] = max_excl
            elif 'Communication' in func_name or 'MPI' in func_name:
                timings['communication'] = max_excl
            elif 'mass_energy' in func_name.lower():
                timings['mass_energy'] = max_excl
            elif 'INIT' in func_name:
                timings['init'] = max_excl
    
    return {
        'mode': 'GPU' if is_gpu else 'CPU',
        'n_tasks': n_tasks,
        'n_threads': n_threads,
        'n_nodes': n_nodes,
        'nx_global': nx_global,
        'nz': nz,
        'sim_time': sim_time,
        'total_cpu_time': total_cpu_time,
        'timings': timings
    }

def create_comparison_plot(cpu_data, gpu_data):
    """Create comparison visualization between CPU and GPU runs."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # === Plot 1: Stacked Bar Chart ===
    categories = ['Thermal', 'Hydrostatic', 'Step', 'Communication', 'Init']
    
    cpu_times = [
        cpu_data['timings']['thermal'],
        cpu_data['timings']['hydrostatic'],
        cpu_data['timings']['step'],
        cpu_data['timings']['communication'],
        cpu_data['timings']['init']
    ]
    
    gpu_times = [
        gpu_data['timings']['thermal'],
        gpu_data['timings']['hydrostatic'],
        gpu_data['timings']['step'],
        gpu_data['timings']['communication'],
        gpu_data['timings']['init']
    ]
    
    x = np.arange(len(categories))
    width = 0.35
    
    colors_cpu = '#E63946'  # Red for CPU
    colors_gpu = '#06FFA5'  # Green for GPU
    
    bars1 = ax1.bar(x - width/2, cpu_times, width, label='CPU', color=colors_cpu, alpha=0.8)
    bars2 = ax1.bar(x + width/2, gpu_times, width, label='GPU', color=colors_gpu, alpha=0.8)
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0.1:  # Only show label if bar is visible
                ax1.text(bar.get_x() + bar.get_width()/2., height,
                        f'{height:.1f}s',
                        ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax1.set_xlabel('Computation Component', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Time (seconds)', fontsize=12, fontweight='bold')
    ax1.set_title('CPU vs GPU: Component-wise Timing Breakdown\n(Max Exclusive Times)', 
                  fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(categories, rotation=15, ha='right')
    ax1.legend(fontsize=11)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    
    # === Plot 2: Speedup Analysis ===
    speedups = []
    labels = []
    
    for cat, cpu_t, gpu_t in zip(categories, cpu_times, gpu_times):
        if gpu_t > 0:
            speedup = cpu_t / gpu_t
            speedups.append(speedup)
            labels.append(cat)
    
    # Add total time speedup
    cpu_total = cpu_data['total_cpu_time']
    gpu_total = gpu_data['total_cpu_time']
    if cpu_total and gpu_total:
        total_speedup = cpu_total / gpu_total
        speedups.append(total_speedup)
        labels.append('TOTAL')
    
    colors = ['#2E86AB'] * (len(speedups) - 1) + ['#A23B72']  # Different color for total
    
    bars = ax2.barh(labels, speedups, color=colors, alpha=0.8)
    
    # Add speedup values on bars
    for i, (bar, speedup) in enumerate(zip(bars, speedups)):
        width = bar.get_width()
        ax2.text(width + 0.1, bar.get_y() + bar.get_height()/2.,
                f'{speedup:.2f}x',
                ha='left', va='center', fontsize=11, fontweight='bold')
    
    # Add reference line at 1x
    ax2.axvline(x=1, color='red', linestyle='--', linewidth=2, alpha=0.7, label='No speedup')
    
    ax2.set_xlabel('Speedup Factor (CPU time / GPU time)', fontsize=12, fontweight='bold')
    ax2.set_title('Speedup Analysis: CPU → GPU\n(Higher is Better)', 
                  fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(axis='x', alpha=0.3, linestyle='--')
    
    # Add configuration info
    cpu_config = f"CPU: {cpu_data['n_tasks']} MPI × {cpu_data['n_threads']} OMP = {cpu_data['n_tasks'] * cpu_data['n_threads']} cores"
    gpu_config = f"GPU: {gpu_data['n_tasks']} MPI tasks × {gpu_data['n_nodes']} nodes"
    problem_size = f"Problem: nx={cpu_data['nx_global']}, nz={cpu_data['nz']}, t={cpu_data['sim_time']}s"
    
    fig.text(0.5, 0.02, f"{cpu_config}  |  {gpu_config}  |  {problem_size}", 
             ha='center', fontsize=10, style='italic')
    
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    
    return fig

def print_detailed_comparison(cpu_data, gpu_data):
    """Print detailed comparison statistics."""
    print("\n" + "="*100)
    print("DETAILED PERFORMANCE COMPARISON: CPU vs GPU")
    print("="*100)
    
    print("\n--- CONFIGURATION ---")
    print(f"CPU Setup: {cpu_data['n_tasks']} MPI tasks × {cpu_data['n_threads']} OpenMP threads = {cpu_data['n_tasks'] * cpu_data['n_threads']} total cores")
    print(f"GPU Setup: {gpu_data['n_tasks']} MPI tasks across {gpu_data['n_nodes']} nodes")
    print(f"Problem Size: nx_global={cpu_data['nx_global']}, nz={cpu_data['nz']}")
    print(f"Simulation Time: {cpu_data['sim_time']}s (CPU), {gpu_data['sim_time']}s (GPU)")
    
    print("\n--- TIMING BREAKDOWN (Max Exclusive Times in seconds) ---")
    print(f"{'Component':<25} {'CPU Time':<15} {'GPU Time':<15} {'Speedup':<12} {'GPU %':<12}")
    print("-"*100)
    
    components = [
        ('Initialization', 'init'),
        ('Thermal Computation', 'thermal'),
        ('Hydrostatic Constant', 'hydrostatic'),
        ('Time Step', 'step'),
        ('MPI Communication', 'communication'),
        ('Mass/Energy Check', 'mass_energy')
    ]
    
    for name, key in components:
        cpu_t = cpu_data['timings'][key]
        gpu_t = gpu_data['timings'][key]
        speedup = cpu_t / gpu_t if gpu_t > 0 else 0
        gpu_pct = (gpu_t / gpu_data['total_cpu_time'] * 100) if gpu_data['total_cpu_time'] else 0
        
        print(f"{name:<25} {cpu_t:<15.2f} {gpu_t:<15.2f} {speedup:<12.2f}x {gpu_pct:<12.1f}%")
    
    print("-"*100)
    cpu_total = cpu_data['total_cpu_time']
    gpu_total = gpu_data['total_cpu_time']
    total_speedup = cpu_total / gpu_total if gpu_total else 0
    
    print(f"{'TOTAL EXECUTION TIME':<25} {cpu_total:<15.2f} {gpu_total:<15.2f} {total_speedup:<12.2f}x")
    print("="*100)
    
    print("\n--- PERFORMANCE METRICS ---")
    print(f"Overall Speedup: {total_speedup:.2f}x")
    print(f"Time Saved: {cpu_total - gpu_total:.2f} seconds ({(1 - gpu_total/cpu_total)*100:.1f}% reduction)")
    
    # Communication overhead analysis
    cpu_comm_pct = (cpu_data['timings']['communication'] / cpu_total * 100) if cpu_total else 0
    gpu_comm_pct = (gpu_data['timings']['communication'] / gpu_total * 100) if gpu_total else 0
    
    print(f"\nCommunication Overhead:")
    print(f"  CPU: {cpu_comm_pct:.1f}% of total time")
    print(f"  GPU: {gpu_comm_pct:.1f}% of total time")
    
    print("\n" + "="*100)

def main():
    # File paths
    cpu_file = Path('/home/raiong/Documents/Github/fightClub/code/results/gpu/best_cpu.out')
    gpu_file = Path('/home/raiong/Documents/Github/fightClub/code/results/gpu/multigpu_scaling_27508614.out')
    
    # Parse files
    print("Parsing CPU log file...")
    cpu_data = parse_log_file(cpu_file)
    
    print("Parsing GPU log file...")
    gpu_data = parse_log_file(gpu_file)
    
    if not cpu_data or not gpu_data:
        print("Error: Could not parse one or both log files!")
        return
    
    # Print detailed comparison
    print_detailed_comparison(cpu_data, gpu_data)
    
    # Create visualization
    print("\nGenerating comparison plots...")
    fig = create_comparison_plot(cpu_data, gpu_data)
    
    # Save figure
    output_file = Path('/home/raiong/Documents/Github/fightClub/code/results/gpu/cpu_vs_gpu_comparison.png')
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    
    plt.show()

if __name__ == '__main__':
    main()
"""
Script for performing a Design of Experiments (DoE) on GROMACS Molecular Dynamics (MD) simulations.
This script:
1. Generates a full-factorial DoE for MD parameters.
2. Runs MD simulations for all combinations.
3. Monitors GPU usage during each run.
4. Collects performance data and saves results to CSV.
"""

# =========================
# Import required libraries
# =========================
import subprocess  # Execute shell commands (GROMACS, file operations)
import os          # File and directory operations
import time        # Track script runtime and timestamps
from pyDOE3 import fullfact  # Full factorial experimental design
import numpy as np
import pandas as pd
import GPUtil       # Monitor GPU load and memory
import csv          # Save GPU monitoring data

# =========================
# Initialize logging
# =========================
start_script = time.time()
logfile = open('doe.log', 'w', encoding='utf-8')
logfile.write('Creating DoE factors and levels (1/4)\n')

# =========================
# Define MD run parameters
# =========================
gpu = 2  # Number of GPUs to use. If different from 2, change this value.

# Define possible ntmpi (MPI processes per GPU) values
ntmpi_2 = [2, 4 , 6, 8, 10]
ntmpi_3 = [3, 6, 9, 12, 15]
ntmpi_4 = [4, 8, 12, 16, 20]

# Define GPU task distribution patterns (string representations)
gputasks_2= ['01', '0011', '0001', '000111', '000011', '000001', '00001111', '00000111', '00000011', '00000001', '0000011111', '0000001111', '0000000111', '0000000011', '0000000001']
gputasks_3 = ['012', '001122', '000112', '000111222', '000011112', '000011112222', '000000111112']
gputasks_4 = ['0123', '00112233', '00011223', '000111222333', '000011112223', '0000111122223333', '0000011111222223']

# Assign ntmpi and gputasks based on number of GPUs
if gpu == 2:
    ntmpi = ntmpi_2
    gputasks = gputasks_2
elif gpu == 3:
    ntmpi = ntmpi_3
    gputasks = gputasks_3
elif gpu == 4:
    ntmpi = list(set(ntmpi_2 + ntmpi_4))
    gputasks = list(set(gputasks_2 + gputasks_4))
print("ntmpi options:", ntmpi)
print("GPU task options:", gputasks)

# =========================
# Experimental design setup
# =========================
logfile.write('Generating full-factorial design (2/4)\n')

# Define the MD parameters and their possible values
mdrun_parameters= {
    'ntmpi' : ntmpi,
    'pin': ['off', 'on'],          # CPU pinning on/off
    'gputasks': gputasks,          # GPU task distribution
    'bonded': ['cpu', 'gpu'],      # Where bonded interactions are calculated
    'update': ['cpu', 'gpu']       # Where particle positions are updated
        }

# List of factors and number of levels for full factorial design
factors = list(mdrun_parameters.keys())
levels = [len(mdrun_parameters[p]) for p in mdrun_parameters]

# Generate full factorial design (indices)
design = fullfact(levels).astype(object)

# Replace indices with actual values for each factor
for i, factor in enumerate(factors):
    design[:, i] = [mdrun_parameters[factor][int(level)] for level in design[:, i]]

# Convert design matrix to DataFrame and save
doe = pd.DataFrame(design, columns=factors)
doe.to_csv('gmx_fullfact_doe.csv', index=False)
print(doe)

# =========================
# Perform MD runs
# =========================
logfile.write('Performing Molecular Dynamics simulations (3/4)\n')

# Load the DOE file (ensure all values are strings for command construction)
doe = pd.read_csv('gmx_fullfact_doe.csv', dtype=str)
runs_num = len(doe.index)
interval = 0.1  # GPU monitoring interval in seconds

# Prepare GROMACS input
# Adjust path if .mdp, .gro, .top files are located elsewhere
subprocess.run('gmx grompp -f ../../*.mdp -o md.tpr -c ../../*.gro -p ../../*.top -v -maxwarn 2', shell=True, check=True)

# Loop over each DoE run
for run in range(runs_num):
    # Create and access the working directory
    os.mkdir(f'{run+1}')
    os.chdir(f'{run+1}')
    
    try:
        # Copy MD input files into the run directory
        subprocess.run("cp ../md* .", shell=True, check=True)
        
        # Construct GROMACS command with parameters from DOE
        cmd = f"gmx mdrun -v -notunepme -resethway -deffnm md -nb gpu -nsteps 25000 -npme 1 -ntmpi {doe.loc[run, 'ntmpi']} -pin {doe.loc[run, 'pin']} -gputasks {doe.loc[run, 'gputasks']} -pme gpu -pmefft gpu -bonded {doe.loc[run, 'bonded']} -update {doe.loc[run, 'update']}"
        logfile.write(f'Command: {cmd}\n')
        
        # Prepare GPU monitoring CSV
        gpu_usage_file = f"gpu_usage_run{run+1}.csv"
        with open(gpu_usage_file, "w", newline="") as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(["Time", "GPU ID", "Load (%)", "Memory Used (MB)", "Memory Total (MB)"])
            
            # Start GROMACS simulation
            process = subprocess.Popen(cmd, shell=True)
            
            # Monitor GPU usage while simulation is running
            while process.poll() is None:
                gpus = GPUtil.getGPUs()
                for gpu in gpus:
                    gpu_load = gpu.load * 100
                    memory_used = gpu.memoryUsed
                    memory_total = gpu.memoryTotal
                    current_time = time.strftime('%H:%M:%S')
                    csv_writer.writerow([current_time, gpu.id, gpu_load, memory_used, memory_total])
                time.sleep(interval)
            
            # Wait for process to finish
            process.wait()
    
    except Exception as e:
        print(f'Cannot execute run {run+1}/{runs_num}: {e}')
        logfile.write(f'Cannot execute run {run+1}/{runs_num}: {e}\n')
    
    # Clean up unnecessary files (keep log and GPU usage)
    file_list = os.listdir()
    for file_name in file_list:
        if file_name not in ['md.log', gpu_usage_file]:
            try:
                os.remove(file_name)
                print(f"File '{file_name}' successfully removed.")
            except Exception as e:
                print(f"Error occurred while removing file '{file_name}': {e}")
    
    # Return to parent directory for next run
    os.chdir('../')

logfile.write(f'\n{runs_num} Molecular Dynamics runs performed\n')

# =========================
# Collect performance data
# =========================
logfile.write('Collecting performances metrics (4/4)\n')
performance = np.array([])

for run in range(runs_num):
    os.chdir(f'{run+1}')
    out_cmd = ["grep", "Performance", "md.log"]
    try:
        # Extract performance from MD log using grep
        output = float(subprocess.check_output(out_cmd, universal_newlines=True).split('\n', maxsplit = 1)[0].split()[1])
        performance = np.append(performance, output)
    except subprocess.CalledProcessError:
        logfile.write(f'No matches found in md.log in directory {run+1}')
        performance = np.append(performance, 0)
    os.chdir('../')
# Moving run directories in runs/
np.set_printoptions(threshold=np.inf)
logfile.write(str(performance))
# Save performance to DOE results
doe['Performance (ns/day)'] = performance
doe.to_csv('gmx_fullfact_doe_results.csv', index=False)

# =========================
# End of script
# =========================
end_script = time.time()
logfile.write(f'\nDoE runs performed in {(end_script-start_script)/60:.3f} min')
logfile.close()
print(f'DoE runs performed in {(end_script-start_script)/60:.3f} min')

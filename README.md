# GROMODEX

**GROMODEX** is a Python script designed to automate a **Design of Experiments (DoE)** for **GROMACS** Molecular Dynamics (MD) performance optimization.  
It generates a full-factorial experimental design for MD parameters, runs the simulations, monitors GPU usage, and collects performance metrics automatically.

---

## Features

1. **Full-factorial DoE generation**  
   Generates all combinations of MD parameters, including:
   - Number of MPI processes per GPU (`ntmpi`)
   - CPU pinning configuration (`on/off`)
   - GPU task distribution patterns
   - CPU/GPU assignment for bonded and update calculations

2. **Automated MD runs**  
   Executes all parameter combinations using GROMACS, creating a dedicated directory for each simulation.

3. **GPU monitoring**  
   Tracks GPU load and memory usage in real time during each run, saving data to CSV files.

4. **Performance collection**  
   Extracts performance metrics (`ns/day`) from GROMACS log files and stores them in a summary CSV file.

5. **Detailed logging**  
   Creates a `doe.log` file with step-by-step progress, executed commands, and potential errors.

---

## Requirements

- **Python** ≥ 3.8  
- [**GROMACS**](http://www.gromacs.org/) installed and accessible via command line  
- **Python dependencies:**
  ```bash
  pip install pyDOE3 numpy pandas GPUtil
  ```
- Hardware: 2–4 GPUs (configurable in the script)

---

## Usage
1. Prepare input files
Ensure your .mdp, .gro, and .top files are located in the correct directory.
If they are elsewhere, adjust their relative paths inside the script.

2. Configure GPU settings
At the top of the script, set the number of GPUs to use:

gpu = 2  # Number of GPUs (supported: 2, 3, or 4)
3. Run the script
Execute the script with:
  
```bash
python gromodex.py
```

## Output Files
File	Description
gmx_fullfact_doe.csv	Full-factorial design of experiments
gmx_fullfact_doe_results.csv	DOE results including performance metrics
gpu_usage_runX.csv	GPU usage data recorded for each run
doe.log	Log file containing workflow progress and errors

---

## Workflow Overview
Generate all MD parameter combinations using a full-factorial design

Prepare simulation inputs with gmx grompp

For each design point:

Create a dedicated run directory

Copy necessary input files

Execute gmx mdrun with parameters defined by the DoE

Monitor GPU load and memory usage

Remove unnecessary temporary files

Extract simulation performance (ns/day) from each md.log file

Aggregate results into a summary CSV file for analysis

---

## Notes
The script supports configurations with 2 to 4 GPUs.

GPU task distributions and MPI settings are predefined for each GPU count.

The number of simulation steps (-nsteps 25000) and other GROMACS flags can be customized directly in the script.

The system must include the grep command (used for performance extraction).

Ensure that all runs have sufficient permissions to create and modify directories.

---

 ## Author
#### Marco Savioli 
Department of Mathematics, University of Rome Tor Vergata, Via della Ricerca Scientifica 1, Rome 00133, Italy.
#### Paolo Calligari
Department of Chemical Sciences and Technology, University of Rome Tor Vergata, Via della Ricerca Scientifica 1, Rome 00133, Italy.
#### Ugo Locatelli
Department of Mathematics, University of Rome Tor Vergata, Via della Ricerca Scientifica 1, Rome 00133, Italy.
#### Gianfranco Bocchinfuso 
Department of Chemical Sciences and Technology, University of Rome Tor Vergata, Via della Ricerca Scientifica 1, Rome 00133, Italy.

---

## Citation / Reference
If you use GROMODEX in your research or publication, please cite it as:

GROMODEX: Optimisation of GROMACS Performance through a Design of Experiment Approach
Marco Savioli, Paolo Calligari, Ugo Locatelli, Gianfranco Bocchinfuso
bioRxiv 2025.10.08.681202; doi: https://doi.org/10.1101/2025.10.08.681202

---

## License
This project is released without a specific license.
You are free to use, modify, and adapt it for research and educational purposes.
If used in publications or derivative works, please provide proper attribution to the author.

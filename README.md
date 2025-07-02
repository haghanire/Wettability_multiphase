# Spatial Characterization of Wetting in Porous Media Using Local Lattice‑Boltzmann Simulations

A comprehensive Python workflow for measuring in situ wettability properties (contact angles) in porous rock samples using multiphase fluid flow simulations and optimization techniques. The numerical method for interface-capturing and hydrodynamic equations is the lattice Boltzmann method (LBM).

## Overview

This workflow performs automated wettability characterization of ganglia in porous media by:

1. **Ganglion Detection**: Identifying and labeling individual ganglia (e.g., oil or water in two-phase systems) in 3D rock samples
2. **Contact Line Analysis**: Detecting three-phase contact lines (solid-water-oil interfaces)
3. **Wettability Optimization**: Using constrained optimization to determine an optimal wettability parameter (contact angle) for each ganglion
4. **Spatial Interpolation**: Filling gaps in the wettability map using interpolation techniques
5. **Visualization Export**: Converting results to VTK format for 3D visualization

## Scientific Background

Wettability is a crucial property in petroleum engineering and geoscience that describes the tendency of a fluid to spread on or adhere to a solid surface in the presence of other immiscible fluids. Various approaches were followed by different researchers to measure local contact angles in porous media to characterize wettability. Recently, a new workflow has been developed and the results have been published in the paper [*Spatial Characterization of Wetting in Porous Media Using Local Lattice‑Boltzmann Simulations*](https://link.springer.com/article/10.1007/s11242-023-02044-x). This workflow uses:

- **Lattice Boltzmann Method (LBM)**: For simulating multiphase fluid flow in porous media using the open-source package **LBPM**
- [**LBPM (Lattice Boltzmann Method for Porous Media)**](https://github.com/OPM/LBPM): an open source software framework designed to model flow processes based on digital rock physics  
- **Constrained Nelder-Mead Optimization**: For finding optimal wettability parameters
- **Spatial Interpolation**: For creating complete wettability maps

## Requirements

### Software Dependencies

```bash
# Core Python packages
numpy
scipy
matplotlib
scikit-image

# Specialized packages
constrNMPy      # Constrained Nelder-Mead optimization
hickle          # HDF5-based pickle replacement
pyevtk          # VTK export functionality
```

### External Software

- **LBPM**: Lattice Boltzmann simulation package
  - Must be compiled and accessible at: `/path/to/lbpm/bin/lbpm_color_simulator`
  - Requires MPI for parallel execution

### System Requirements

- **Memory**: Minimum 8GB RAM (16GB+ recommended for large domains)
- **Storage**: ~1GB free space for intermediate files
- **OS**: Linux/Unix (tested on Ubuntu 20.04.4 LTS)
- **CPU**: Multi-core recommended for optimization

## Installation

1. **Clone or download the workflow**:
   ```bash
   git clone <this repository url>
   cd Wettability_multiphase
   ```

2. **Install Python dependencies**:
   ```bash
   pip install numpy scipy matplotlib scikit-image hickle pyevtk
   pip install constrNMPy  # May require compilation
   ```

3. **Install LBPM**:
   ```bash
   # Follow LBPM installation instructions (https://github.com/OPM/LBPM/wiki)
   # Ensure the binary is accessible at the specified path
   ```

4. **Verify installation**:
   ```bash
   python -c "import numpy, scipy, constrNMPy, hickle, pyevtk; print('All dependencies installed successfully')"
   ```

## Usage

### Quick Start

1. **Prepare input files**:
   - Place your 3D rock sample in the working directory as `image.raw`
   - Ensure LBPM input template `input.db` is available (an example follows)
   - Modify the parameters related to the `Domain` in the input file `input.db`
   - Verify file paths `path` and `lbpm_binary` in the script match your system

2. **Run the complete workflow**:
   ```bash
   python wettability_characterization.py
   ```

3. **Monitor progress**:
   - The script will display progress for each stage
   - Intermediate results are saved automatically
   - Check log output for any errors

### Input Files

| File | Format | Description |
|------|--------|-------------|
| `image.raw` | Binary (uint8) | 3D rock sample (100×100×100 voxels) |
| `input.db` | Text | LBPM configuration template |
| `solid_fluid.raw` | Binary (uint8) | Pre-computed solid-fluid interfaces (optional) |
| `contact_line.raw` | Binary (uint8) | Pre-computed contact lines (optional) |

### Input Data Format

The main input file `image.raw` should contain (numbers as an example):
- **Dimensions**: 100×100×100 voxels (adjustable in script)
- **Data type**: 8-bit unsigned integers
- **Phase encoding**:
  - `0`: Solid (rock)
  - `1`: Water
  - `2`: Oil

The selection of the ganglion (either water or oil) can be done in the labeling part of `ndimage.measurements.label`. 

### File Structure

```
project_directory/
├── # Input files
│   ├── wettability_characterization.py   # Main workflow script
│   ├── input.db                          # LBPM configuration template
│   ├── image.raw                         # Input 3D rock sample
│   ├── README.md                         # This file
│   └── LICENSE                           # License    
├── # Generated output files
│   ├── affinity_map_matrix.npy           # Generated affinity
│   ├── affinity_map_processed_matrix.npy # Propagated affinity
│   ├── affinity_map_processed_vtk.vtr    # Affinity saved in VTK format
│   ├── results                           # Pickled optimization results
│   └── *.txt                             # Optimization data files
└── # Temporary LBPM files    
    ├── subdomain_confined.raw   
    ├── subdomain.raw            
    └── input_mod.db             
```

## Output Files

### Primary Outputs

| File | Format | Description |
|------|--------|-------------|
| `affinity_map_processed_matrix.npy` | NumPy array | Final interpolated wettability map |
| `affinity_map_processed_vtk.vtr` | VTK | 3D visualization file |
| `results` | Hickle | Complete optimization results |

### Intermediate Outputs

| File | Format | Description |
|------|--------|-------------|
| `affinity_map_matrix.npy` | NumPy array | Raw wettability map (before interpolation) |
| `affinity_vector.txt` | Text | Optimized affinity values |
| `error_vector.txt` | Text | Optimization errors |
| `ganglion_volume.txt` | Text | Ganglion volumes |
| `ganglion_count.txt` | Text | Ganglion identifiers |
| `contact_line.raw` | Binary | Three-phase contact line map |
| `solid_fluid.raw` | Binary | Solid-fluid interface map |

## Workflow Stages

### Stage 1: Wettability Optimization

- **Input**: 3D rock sample, LBPM configuration
- **Process**: 
  - Detect and label oil (default) ganglia
  - Identify three-phase contact lines
  - For each ganglion, optimize wettability using LBPM simulations
- **Output**: Raw wettability map with optimized values at contact lines

### Stage 2: Post-processing

- **Input**: Raw wettability map, solid-fluid interface map
- **Process**:
  - Identify locations needing interpolation
  - Apply linear interpolation where possible
  - Use nearest-neighbor interpolation as fallback
- **Output**: Complete wettability map with interpolated values

### Stage 3: VTK Export

- **Input**: Processed wettability map
- **Process**: Convert NumPy array to VTK structured grid format
- **Output**: VTK file for visualization in ParaView, VisIt, etc.

## Configuration

### Key Parameters

Edit these variables in the script to customize the analysis:

```python
# Domain size (must match input.db and image.raw information)
original_domain_size = np.array((100, 100, 100), dtype=int)

# File paths
path = '/path/to/your/working_directory/'
lbpm_binary = '/path/to/lbpm/bin/lbpm_color_simulator'
timeout = 120  # Increase for large ganglia

# Optimization bounds
LB = [-1.0]  # Lower bound for wettability
UB = [1.0]   # Upper bound for wettability
xtol = [0.001]  # Optimization tolerance

# Ganglion size threshold
mgv = 10  # Minimum ganglion volume in voxels to process
bf = 10  # Boundary buffer
```

### LBPM Configuration

The `input.db` file contains LBPM simulation parameters. Key settings:

```
Domain {
   Filename = "image.raw"  
   ReadType = "8bit"        // data type
   nproc = 1, 1, 1          // process grid
   n = 30, 26, 26           // sub-domain size
   N = 30, 26, 26           // size of original_image
   voxel_length = 3.58      // voxel length (in microns)
   ReadValues = 0, 1, 2     // labels within the original image
   WriteValues = 0, 1, 2    // associated labels to be used by LBPM
   offset = 0, 0, 0         // offset from the origin for the simulation region
   InletLayers = 0, 0, 0    // specify 10 layers along the z-inlet
   checkerSize = 0          // size of the checker to use
   BC = 4                   // boundary condition type (0 for periodic, 1 for P, 2 for velocity, 3 for dynamic P, 4 for flux)
}
Color {
    tauA = 0.7;             // relaxation time for fluid A (labeled as "1")       
    tauB = 0.7;             // relaxation time for fluid B (labeled as "2") 
    rhoA   = 1.0;           // density for fluid A (in lattice units)
    rhoB   = 1.0;           // density for fluid B (in lattice units)
    alpha = 1e-3;           // controls the surface tension
    beta  = 0.95;           // controls the interface width 
    F = 0, 0, 0             // controls the external force
    Restart = false         // initialize simulation from restart file?
    timestepMax = 10000     // maximum number of timesteps to perform before exit
    ComponentLabels = 0     // number of immobile component labels in the input image
    ComponentAffinity = 0.9 // wetting condition for each immobile component
    flux = 0.0              // volumetric flux at the z-inlet in voxels per timestep
}
Analysis {
    analysis_interval = 250         // Frequency to perform analysis
    visualization_interval = 1000   // Frequency to write visualization data
    restart_interval = 10000        // Frequency to write restart data
    restart_file = "Restart"        // Filename to use for restart file (will append rank)
    N_threads    = 1                // Number of threads to use for analysis
    load_balance = "independent"    // Load balance method to use: "none", "default", independent"
}
Visualization {
    write_silo = false          // write SILO databases with assigned variables
    save_8bit_raw = true        // write labeled 8-bit binary files with phase assignments
    save_phase_field = false    // save phase field within SILO database
    save_pressure = false       // save pressure field within SILO database
    save_velocity = false       // save velocity field within SILO database
}
FlowAdaptor {

}
```

## Visualization

### Using ParaView

1. **Open VTK file**:
   ```bash
   paraview affinity_map_processed_vtk.vtr
   ```

2. **Visualization tips**:
   - Use "Slice" filter to examine cross-sections
   - Apply "Threshold" filter to isolate specific wettability ranges
   - Use "Contour" filter to create isosurfaces
   - Color by "affinity_map_processed_vtk" variable

### Using Python

```python
import numpy as np
import matplotlib.pyplot as plt

# Load processed data
affinity_map = np.load('affinity_map_processed_matrix.npy')

# Create 2D slice visualization
plt.figure(figsize=(10, 8))
plt.imshow(affinity_map[50, :, :], cmap='RdBu', vmin=-1, vmax=1)
plt.colorbar(label='Wettability')
plt.title('Wettability Map - Middle Slice')
plt.show()
```

## Troubleshooting

### Common Issues

1. **LBPM simulation timeout**:
   - Increase timeout value in `subprocess.run()`
   - Check LBPM binary path and permissions
   - Verify MPI installation

2. **Memory errors**:
   - Reduce domain size or process fewer ganglia
   - Increase system RAM or add swap space
   - Process ganglia in batches

3. **Interpolation failures**:
   - Check for sufficient data points
   - Verify contact line detection
   - Adjust interpolation parameters

4. **File not found errors**:
   - Verify all input files exist
   - Check file paths in script
   - Ensure proper file permissions

### Performance Optimization

1. **Parallel processing**:
   - Use multiple MPI processes for LBPM
   - Implement multiprocessing for ganglion analysis

2. **Memory management**:
   - Process ganglia in smaller batches
   - Use memory-mapped arrays for large datasets

3. **Storage optimization**:
   - Use compressed file formats
   - Clean up temporary files regularly

## License

This project is licensed under the The GNU General Public License v3.0 - see the LICENSE file for details.

## Citation

If you use this workflow in your research, please cite:

```bibtex
@article{erfani2024spatial,
  title={Spatial characterization of wetting in porous media using local lattice-Boltzmann simulations},
  author={Erfani, Hamidreza and Haghani, Reza and McClure, James and Boek, Edo and Berg, Carl Fredrik},
  journal={Transport in Porous Media},
  volume={151},
  number={3},
  pages={429--448},
  year={2024},
  publisher={Springer}
}
```

## Support and Contributing

For questions, issues, or contributions:
- Open an issue on GitHub
- Contact: [haghani.re@gmail.com]
- Submit a pull request

---

**Note**: This workflow is designed for research purposes and may require adaptation for specific use cases.
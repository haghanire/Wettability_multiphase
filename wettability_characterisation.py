#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spatial Characterization of Wetting in Porous Media Using Local Lattice-Boltzmann Simulations
=============================================================================


This script combines three main stages:
1. Wettability optimization on oil (default) ganglia (main analysis)
2. Post-processing with interpolation to fill missing values
3. VTK export for visualization

Authors:    Hamidreza Erfani    (written)
            Carl Fredrik Berg   (supervision)
            Reza Haghani        (modified)

Utilized in the following scientific paper:

    Erfani, H., Haghani, R., McClure, J., Boek, E. and Berg, C.F., 
    2024. Spatial characterization of wetting in porous media using local lattice-Boltzmann simulations. 
    Transport in Porous Media, 151(3), pp.429-448. 

"""

# =============================================================================
# IMPORTS
# =============================================================================
from scipy import ndimage
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
import constrNMPy as cNM  # Constrained Nelder-Mead optimization
import numpy as np
import subprocess
from pyevtk.hl import gridToVTK
import hickle as hkl

# =============================================================================
# STAGE 1: MAIN WETTABILITY OPTIMIZATION
# =============================================================================

def wettability_opt(affinity):
    """
    Optimization function that runs LBPM simulation with given affinity values
    and returns error compared to experimental equilibrated image
    """

    # Load equilibrated experimental image for comparison
    filename = 'subdomain_confined.raw'
    f = open(path + filename, 'rb')
    img_eq = np.fromfile(f, dtype=np.uint8)
    img_eq = img_eq.reshape(subdomain_confined_size[0], subdomain_confined_size[1], subdomain_confined_size[2])
    f.close()
    
    # Modify LBPM input file with new affinity value
    fin = open(path + 'input_mod.db', "rt").readlines()
    fout = open(path + 'input_mod.db', "wt")
    
    for line in fin:
        if 'ComponentAffinity' in line:
            fout.write('    ComponentAffinity = ' + str.format("{0:.5f}", affinity[0]) + 
                      ' // wetting condition for each immobile component\r\n')
        else:
            fout.write(line)
    
    # print(affinity.ndim)
    fout.close()
    
    # Run LBPM simulation with timeout
    subprocess.run('timeout 120 ' + 
                  lbpm_binary + 
                  path + 'input_mod.db ; exit', shell=True)
    
    # Read simulation result
    filename = 'id_t4000.raw'
    f = open(path + filename, 'rb')
    img_sim = np.fromfile(f, dtype=np.uint8)
    
    # Check if simulation completed successfully
    if ((np.shape(img_sim)[0]) != (subdomain_confined_size[0] * subdomain_confined_size[1] * subdomain_confined_size[2])):
        img_sim = np.ones((subdomain_confined_size[0], subdomain_confined_size[1], subdomain_confined_size[2]))
    else:
        img_sim = img_sim.reshape(subdomain_confined_size[0], subdomain_confined_size[1], subdomain_confined_size[2])
    
    f.close()
    
    # Calculate error between simulation and experimental data
    error = 0
    for i in range(subdomain_confined_size[0]):
        for j in range(subdomain_confined_size[1]):
            for k in range(subdomain_confined_size[2]):
                if (img_eq[i, j, k] == 2.0) and (img_sim[i, j, k] != img_eq[i, j, k]):
                    error = error + 1   
    # print(error)
    return error


print("="*60)
print("STAGE 1: WETTABILITY OPTIMIZATION")
print("="*60)

# Setup paths and parameters
global path
global lbpm_binary
path = '/path/to/your/working_directory/'
lbpm_binary = '/path/to/lbpm/bin/lbpm_color_simulator '

original_domain_size = np.array((100, 100, 100), dtype=int)

# Load original 3D domain (0:Solid, 1:Water, 2:Oil)
filename = 'image.raw'
f = open(path + filename, 'rb')
original_domain = np.fromfile(f, dtype=np.uint8)
original_domain = original_domain.reshape(original_domain_size[0], original_domain_size[1], original_domain_size[2])
f.close()

# Detect three-phase contact lines (solid-water-oil interfaces)
print("Detecting three-phase contact lines...")
available_contact_raw = 0

if (available_contact_raw == 0):
    contact_line = np.zeros(original_domain_size, dtype=int)
    contact_detector = np.zeros((3, 3, 3), dtype=int)
    
    # Check each voxel for contact with all three phases
    for ii in range(1, original_domain_size[0]):
        for jj in range(1, original_domain_size[1]):
            for kk in range(1, original_domain_size[2]):
                if (original_domain[ii, jj, kk] == 0):  # If solid voxel
                    # Check 3x3x3 neighborhood
                    contact_detector = original_domain[(ii-1):(ii+2), (jj-1):(jj+2), (kk-1):(kk+2)]
                    # Check if both water and oil are present in neighborhood
                    if (np.sum((contact_detector == 1.0).astype(int)[:, :, :]) != 0) and \
                       (np.sum((contact_detector == 2.0).astype(int)[:, :, :]) != 0):
                        contact_line[ii, jj, kk] = 1
    
    contact_line.astype(np.uint8).tofile('contact_line.raw')
else:
    # Load pre-computed contact line data
    filename = 'contact_line.raw'
    f = open(path + filename, 'rb')
    contact_line = np.fromfile(f, dtype=np.uint8)
    contact_line = contact_line.reshape(original_domain_size[0], original_domain_size[1], original_domain_size[2])
    f.close()

# Separate and label different phases
print("Labeling oil ganglia...")
oil_skeleton = (original_domain == 2.0)
water_skeleton = (original_domain == 1.0)
rock_skeleton = (original_domain == 0.0)

print("oil_skeleton.size=", oil_skeleton.shape)
oil_skeleton = oil_skeleton.astype(int) * 2  # 2 is oil
water_skeleton = water_skeleton.astype(int) * 1  # 1 is water
rock_skeleton = rock_skeleton.astype(int)  # 1 is rock

# Label connected oil regions (ganglia)
lO, num_O = ndimage.measurements.label(oil_skeleton)

# Initialize affinity map (2 means no data)
affinity_map = np.zeros(original_domain_size)
affinity_map = affinity_map + 2

# Process each oil ganglion
print(f"Processing {num_O} oil ganglia...")
opt_counter = 0

for i in range(1, num_O):

    print(f"*** Processing ganglion {i}/{num_O}... ***")
    
    # Extract current ganglion
    ganglion = (lO == i)
    ganglion = ganglion.astype(int)
    ganglia_vol = sum(sum(sum(ganglion)))  # Volume in voxels
    
    # Find bounding box of ganglion
    indices = np.nonzero(ganglion)
    min_x = min(indices[0])
    max_x = max(indices[0])
    min_y = min(indices[1])
    max_y = max(indices[1])
    min_z = min(indices[2])
    max_z = max(indices[2])
    bf = 10 # boundary_buffer
    mgv = 10    # Minimum ganglion volume voxels to process
    # Check for contact lines in ganglion region
    subdomain_contacts = contact_line[(min_x-1):(max_x+1), (min_y-1):(max_y+1), (min_z-1):(max_z+1)]
    
    # Process ganglion if it meets size and boundary criteria
    if ((ganglia_vol > mgv) and (min_x) > 2 and (original_domain_size[0] - max_x) > 2 and 
        (min_y) > 2 and (original_domain_size[1] - max_y) > 2 and 
        (min_z) > 2 and (original_domain_size[2] - max_z) > 2 and 
        (sum(sum(sum(subdomain_contacts))) != 0)):
        
        # Create subdomain around ganglion
        if ((min_x) > bf and (original_domain_size[0] - max_x) > bf and 
            (min_y) > bf and (original_domain_size[1] - max_y) > bf and 
            (min_z) > bf and (original_domain_size[2] - max_z) > bf):
            # Use larger subdomain if space allows
            subdomain = original_domain[(min_x-bf):(max_x+bf), (min_y-bf):(max_y+bf), (min_z-bf):(max_z+bf)]
            subdomain = (subdomain != 0)
            subdomain = subdomain.astype(int)
            subdomain = subdomain + ganglion[(min_x-bf):(max_x+bf), (min_y-bf):(max_y+bf), (min_z-bf):(max_z+bf)]
        else:
            # Use minimal subdomain
            subdomain = original_domain[(min_x-1):(max_x+1), (min_y-1):(max_y+1), (min_z-1):(max_z+1)]
            subdomain = (subdomain != 0)
            subdomain = subdomain.astype(int)
            subdomain = subdomain + ganglion[(min_x-1):(max_x+1), (min_y-1):(max_y+1), (min_z-1):(max_z+1)]
        
        opt_counter = opt_counter + 1
        subdomain_size = np.shape(subdomain)
        
        # Create confined subdomain with padding
        subdomain_confined = np.zeros((subdomain_size[0] + (subdomain_size[0] % 2) + 2,
                                     subdomain_size[1] + (subdomain_size[1] % 2) + 2,
                                     subdomain_size[2] + (subdomain_size[2] % 2) + 2))
        
        for ii in range(subdomain_size[0]):
            for jj in range(subdomain_size[1]):
                for kk in range(subdomain_size[2]):
                    subdomain_confined[ii+1, jj+1, kk+1] = subdomain[ii, jj, kk]
        
        # Save subdomain data for LBPM
        global subdomain_confined_size
        subdomain_confined_size = np.shape(subdomain_confined)
        subdomain_confined.astype(np.uint8).tofile('subdomain_confined.raw')
        subdomain.astype(np.uint8).tofile('subdomain.raw')
        
        # Prepare LBPM input file with correct dimensions
        fin = open(path + 'input.db', "rt").readlines()
        fout = open(path + 'input_mod.db', "wt")
        
        for line in fin:
            if 'sub-domain' in line:
                fout.write('   n = ' + str(int(subdomain_confined_size[2])) + ', ' + 
                          str(int(subdomain_confined_size[1])) + ', ' + 
                          str(int(subdomain_confined_size[0])) + ' // sub-domain size\r\n')
            elif 'original_image' in line:
                fout.write('   N = ' + str(subdomain_confined_size[2]) + ', ' + 
                          str(subdomain_confined_size[1]) + ', ' + 
                          str(subdomain_confined_size[0]) + ' // size of original_image\r\n')
            else:
                fout.write(line)
        
        fout.close()
        
        # Run optimization to find best wettability
        print(f"Optimizing wettability for ganglion {i}...")
        if opt_counter == 1:
            x0 = [0.50]  # Initial guess
        else:
            x0[0] = (res['xopt'][0])  # Use previous result as starting point
        
        LB = [-1.0]  # Lower bound
        UB = [1.0]   # Upper bound
        
        # Run constrained optimization
        res = cNM.constrNM(wettability_opt, x0, LB, UB, ftol=[1.0], xtol=[0.001], full_output=True)
        print(f"Optimal affinity: {res['xopt']}")
        print(f"Ganglion {i} completed")
        
        # Store results
        if opt_counter == 1:
            affinity_vector = np.array(res['xopt'][0])
            error_vector = np.array(res['fopt'])
            ganglion_volume = np.array([ganglia_vol])
            ganglion_count = np.array([i])
        else:
            affinity_vector = np.append(affinity_vector, res['xopt'][0])
            error_vector = np.append(error_vector, res['fopt'])
            ganglion_count = np.append(ganglion_count, [i])
            ganglion_volume = np.append(ganglion_volume, [ganglia_vol])
            
            # Save intermediate results
            np.savetxt('affinity_vector.txt', affinity_vector, delimiter=',')
            np.savetxt('error_vector.txt', error_vector, delimiter=',')
            np.savetxt('ganglion_volume.txt', ganglion_volume, delimiter=',')
            np.savetxt('ganglion_count.txt', ganglion_count, delimiter=',')
        
        # Map optimized affinity values to contact line locations
        for ii in range(max_x - min_x + 1):
            for jj in range(max_y - min_y + 1):
                for kk in range(max_z - min_z + 1):
                    if (contact_line[(ii + min_x), (jj + min_y), (kk + min_z)] == 1):
                        affinity_map[(ii + min_x), (jj + min_y), (kk + min_z)] = (res['xopt'][0])

# Save final results
print("Saving Stage 1 results...")
np.save('affinity_map_matrix', affinity_map)
hkl.dump([ganglion_count, ganglion_volume, error_vector, affinity_vector], 'results')

print(f"Stage 1 completed. Processed {opt_counter} ganglia.")

# =============================================================================
# STAGE 2: POST-PROCESSING WITH INTERPOLATION
# =============================================================================

print("\n" + "="*60)
print("STAGE 2: POST-PROCESSING WITH INTERPOLATION")
print("="*60)

# Load the affinity map from Stage 1
affinity_map = np.load('affinity_map_matrix.npy')

# Load or create solid-fluid interface map
available_solid_fluid_raw = 1

if (available_solid_fluid_raw == 0):
    print("Creating solid-fluid interface map...")
    solid_fluid = np.zeros(original_domain_size, dtype=int)
    contact_detector = np.zeros((3, 3, 3), dtype=int)
    
    # Detect solid-fluid interfaces
    for ii in range(1, original_domain_size[0] - 1):
        for jj in range(1, original_domain_size[1] - 1):
            for kk in range(1, original_domain_size[2] - 1):
                if (rock_skeleton[ii, jj, kk] == 1):
                    contact_detector = rock_skeleton[(ii-1):(ii+2), (jj-1):(jj+2), (kk-1):(kk+2)]
                    if (np.sum(contact_detector) != 27):  # Not completely surrounded by rock
                        solid_fluid[ii, jj, kk] = 1
    
    solid_fluid.astype(np.uint8).tofile('solid_fluid.raw')
else:
    # Load pre-computed solid-fluid interface
    print("Loading solid-fluid interface map...")
    filename = 'solid_fluid.raw'
    f = open(path + filename, 'rb')
    csolid_fluid = np.fromfile(f, dtype=np.uint8)
    solid_fluid = csolid_fluid.reshape(original_domain_size[0], original_domain_size[1], original_domain_size[2])
    f.close()

# Prepare data for interpolation
print("Preparing interpolation data...")
x = []
y = []
z = []
affinity_sim = []

# Collect known affinity values (where affinity_map != 2)
for ii in range(1, original_domain_size[0] - 1):
    for jj in range(1, original_domain_size[1] - 1):
        for kk in range(1, original_domain_size[2] - 1):
            if (affinity_map[ii, jj, kk] != 2):
                x.append(ii)
                y.append(jj)
                z.append(kk)
                affinity_sim.append(affinity_map[ii, jj, kk])

# Create interpolators
print("Creating interpolators...")
interp_lin = LinearNDInterpolator(list(zip(x, y, z)), affinity_sim)
interp_near = NearestNDInterpolator(list(zip(x, y, z)), affinity_sim)

# Apply interpolation to fill missing values
print("Applying interpolation...")
affinity_map_processed = np.zeros(original_domain_size)

for ii in range(1, original_domain_size[0] - 1):
    for jj in range(1, original_domain_size[1] - 1):
        for kk in range(1, original_domain_size[2] - 1):
            if (solid_fluid[ii, jj, kk] == 1):  # Only at solid-fluid interfaces
                if (affinity_map[ii, jj, kk] != 2):
                    # Use original value if available
                    affinity_map_processed[ii, jj, kk] = affinity_map[ii, jj, kk]
                else:
                    # Use linear interpolation first
                    affinity_map_processed[ii, jj, kk] = interp_lin(ii, jj, kk)
                    # If linear interpolation fails, use nearest neighbor
                    if (np.isnan(affinity_map_processed[ii, jj, kk]) == 1):
                        affinity_map_processed[ii, jj, kk] = interp_near(ii, jj, kk)

# Save processed affinity map
print("Saving processed affinity map...")
np.save('affinity_map_processed_matrix', affinity_map_processed)
print("Stage 2 completed.")

# =============================================================================
# STAGE 3: VTK EXPORT FOR VISUALIZATION
# =============================================================================

print("\n" + "="*60)
print("STAGE 3: VTK EXPORT FOR VISUALIZATION")
print("="*60)

def NParray2VTK(npmodel, dims, scale, vtkfilename):
    """
    Convert numpy array to VTK format for visualization
    
    Parameters:
    - npmodel: 3D numpy array to export
    - dims: dimensions [nx, ny, nz]
    - scale: physical scale factor
    - vtkfilename: output filename
    """
    # Get dimensions
    nx, ny, nz = dims[0], dims[1], dims[2]
    
    # Physical lengths (assuming isotropic voxels)
    lx, ly, lz = scale, scale, scale
    
    # Voxel sizes
    dx, dy, dz = lx/nx, ly/ny, lz/nz
    
    # Create coordinate arrays
    X = np.arange(0, lx + 0.1*dx, dx, dtype='float64')
    Y = np.arange(0, ly + 0.1*dy, dy, dtype='float64')
    Z = np.arange(0, lz + 0.1*dz, dz, dtype='float64')
    
    # Initialize coordinate grids
    x = np.zeros((nx + 1, ny + 1, nz + 1))
    y = np.zeros((nx + 1, ny + 1, nz + 1))
    z = np.zeros((nx + 1, ny + 1, nz + 1))
    
    # Fill coordinate grids
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                x[i, j, k] = X[i]
                y[i, j, k] = Y[j]
                z[i, j, k] = Z[k]
    
    # Export to VTK
    vtkfilename = "./%s" % vtkfilename
    gridToVTK(vtkfilename, x, y, z, cellData={vtkfilename: npmodel})
    print(f"VTK file saved as: {vtkfilename}.vtr")

# Export processed affinity map to VTK
print("Exporting to VTK format...")
affinity_map_processed = np.load('affinity_map_processed_matrix.npy')
NParray2VTK(affinity_map_processed, [100, 100, 100], 1, 'affinity_map_processed_vtk')

print("Stage 3 completed.")
print("\n" + "="*60)
print("PIPELINE COMPLETED SUCCESSFULLY!")
print("="*60)
print("Output files:")
print("- affinity_map_matrix.npy (raw affinity map)")
print("- affinity_map_processed_matrix.npy (interpolated affinity map)")
print("- affinity_map_processed_vtk.vtr (VTK file for visualization)")
print("- results (pickled optimization results)")
print("- Various .txt files with optimization data")

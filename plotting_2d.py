import matplotlib.pyplot as plt
import numpy as np
import os
import glob 
import sys
import matplotlib
import math
from scipy.special import hermite
import shutil
import re
script_dir = os.path.dirname(os.path.abspath(__file__))
# Find all __init__.py files in the directory and its subdirectories
init_files = glob.glob(os.path.join(script_dir, '**/__init__.py'), recursive=True)
# List all found __init__.py files
# font

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'font.family': 'sans-serif'})

for file in init_files:
     dir_path = os.path.dirname(file)
     os.chdir(dir_path)
     print(f'Running avg.py in {os.system("pwd")}')
     #os.system('python3 __init__.py')
     os.system('mv output/splitop_tilted-EP1D-0.txt ./splitop_tilted-EP1D.txt')

##Now need to extract Iw from the files
# Get the directory where this script is located
base_path = os.path.dirname(os.path.abspath(__file__))
# Load the model
sys.path.append(base_path+"/Model")
exec(f"import EP1D as model")
par = model.parameters()
ω =  np.array([int(sys.argv[1])])/27.2114 #np.arange(0.0, 5.0, 0.01)/27.2114
def Ct(fname):
    dat = np.loadtxt(fname)

    dtE = par.dtE
    t =  np.arange(len(dat[:,0])) * dtE 
    return t, dat[:,0] + 1j * dat[:,1]
def Iw(t, Ct):

    dt = t[1] - t[0]
    
    I = ω * 0

    for iω, ωj in enumerate(ω): 
        #print(iω, ωj)   
        expiωt = np.exp(1j * ωj * t)
        Ω = np.pi / (2 * t[-1])
        cosωt = np.cos( Ω * t)
        I[iω] = (np.sum(Ct * expiωt * cosωt * dt)).real
    return ω, I
def extract_number(path):
    # Get the last part like 'initk_5' and extract the number
    basename = os.path.basename(path)
    return int(basename.split("_")[1])  # e.g., 'initk_5' → 5

# Loop through all initk_* folders
initk_paths = sorted(glob.glob(os.path.join(base_path, "initk_*")), key=extract_number)
init8_counts = [len(glob.glob(os.path.join(k, "init_*"))) for k in initk_paths]
max_init8 = max(init8_counts)
Full4x = np.zeros((len(initk_paths) * 2, max_init8 * 2), dtype=float) 
I_matrix = np.zeros((len(initk_paths), max_init8),  dtype=float)  
init8_paths_all = []

for i, initk_path in enumerate(initk_paths):
    init8_paths = sorted(glob.glob(os.path.join(initk_path, "init_*")))
    
    init8_paths_all.append(init8_paths)
    for j, init8_path in enumerate(init8_paths):
        target_file = os.path.join(init8_path, "splitop_tilted-EP1D.txt")
        # Remove all other files and folders except splitop_tilted-EP1D.txt
        for item in os.listdir(init8_path):
            item_path = os.path.join(init8_path, item)
            # Only remove if the folder contains splitop_tilted-EP1D.txt
            if os.path.exists(target_file):
                if item != "splitop_tilted-EP1D.txt":
                        os.system(f'rm -r {item_path}')
                        print(f"remove file {item_path}")
        if os.path.exists(target_file):
            t, ct = Ct(target_file)
            ω, I = Iw(t, ct)
            I_matrix[i, j] = I[0] if np.ndim(I) else I
# Plot that I matrix as a 2D plot

Full4x[len(initk_paths):, max_init8:] = I_matrix
Full4x[len(initk_paths):, :max_init8] = I_matrix[:,::-1]
Full4x[:len(initk_paths), :max_init8] = I_matrix[::-1,::-1]
Full4x[:len(initk_paths), max_init8:] = I_matrix[::-1,:]



maxkx = par.kx[int(initk_paths[-1].split("_")[-1]) - int(initk_paths[0].split("_")[-1])] * 1000
maxky = par.kx[int(init8_paths[-1].split("_")[-1]) - int(init8_paths[0].split("_")[-1])] * 1000

cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom',['#020024', '#547bff', '#00d4ff', '#a9f0ff', '#a9ffe2'])
plt.figure(figsize=(5, 5/ 1.034))
plt.imshow(Full4x, aspect='auto',cmap=cmap, interpolation='nearest', extent=[-maxkx, maxkx, -maxky, maxky])
#plt.colorbar(label='Intensity')
plt.xlabel('$k_x$ (1/a × $10^3$)')

plt.ylabel('$k_y$ (1/a × $10^3$)')
plt.clim(0, max(Full4x.flatten()) * 1.2)  # Set color limits to enhance visibility
#plt.title('2D Absorbtion Spectra For Tilted Multilayered System')
plt.tight_layout()
plt.savefig(os.path.join(base_path, "I_matrix_2d_plot.png"))
plt.close()
# Save the I_matrix to a text file
np.savetxt(os.path.join(base_path, "I_matrix.txt"), I_matrix, fmt='%.6f', delimiter='\t',
           header='I matrix for ω = {:.2f} eV'.format(ω[0] * 27.2114))
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
init_files = glob.glob(os.path.join(script_dir, '**/serial.py'), recursive=True)
# List all found __init__.py files
# font

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'font.family': 'sans-serif'})

for file in init_files:
     dir_path = os.path.dirname(file)
     os.chdir(dir_path)
     print(f'Running avg.py in {os.system("pwd")}')
     #os.system('python3 __init__.py')
     os.system('mv output/splitop_tilted-EP1D.txt ./splitop_tilted-EP1D.txt')

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
θ = round(par.θtilt, 2)
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
plt.savefig(os.path.join(base_path, "2D_Abs_Spec_θ=" + str(θ) + ".pdf"))
plt.close()
# Save the I_matrix to a text file
np.savetxt(os.path.join(base_path, "2D_Abs_Spec_θ=" + str(θ) + ".txt"), I_matrix, fmt='%.6f', delimiter='\t',
           header='I matrix for ω = {:.2f} eV'.format(ω[0] * 27.2114))
# Need to plot 1D absorbtion spectra for Ky = 0 and Kx = 0
#Ky = 0
##Find all the init8_* folders with Ky = 0
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
initk0 = glob.glob(os.path.join(base_path, "initk_0", "init_*"))
initk0.sort(key=extract_number)
print(initk0)
ω =  np.arange(0.0, 5.0, 0.01)/27.2114
I2Dky = np.zeros((len(ω), len(initk0)))
def Iw2(t, Ct):

    dt = t[1] - t[0]
    
    I = ω * 0

    for iω, ωj in enumerate(ω): 
        #print(iω, ωj)   
        expiωt = np.exp(1j * ωj * t)
        Ω = np.pi / (2 * t[-1])
        cosωt = np.cos( Ω * t)
        I[iω] = (np.sum(Ct * expiωt * cosωt * dt)).real
    return ω, I
for i, ifile in enumerate(initk0):
    file = os.path.join(ifile, "splitop_tilted-EP1D.txt")
    if not os.path.exists(file):
        print(f"File {file} does not exist, skipping.")
    else:
        t, ct = Ct(os.path.join(ifile, "splitop_tilted-EP1D.txt"))
        ω, I = Iw2(t, ct)
        I2Dky[:,i] = I
        print(f"Processing {ifile} with shape {I.shape}")
file_name = "1D_Abs_Spect_ky=0_θ=" + str(θ)  + ".txt"

np.savetxt(file_name, I2Dky)   
M2 = np.zeros((I2Dky[::-1,:].shape[0], I2Dky[::-1,:].shape[1]*2))
M2[:, :I2Dky[::-1,:].shape[1]] = I2Dky[::-1,::-1]
M2[:, I2Dky[::-1,:].shape[1]:] = I2Dky[::-1,:]
plt.figure(figsize=(5, 5/ 1.034))
plt.imshow(M2, aspect='auto' , extent=[ -maxky, maxky, 0.0, 5.0],cmap=cmap, interpolation='gaussian')
plt.clim(0, 250)
plt.tight_layout()
plt.xlabel('$k_x$ (1/a × $10^3$)')
plt.ylabel('Energy (eV)')
plt.ylim(2.0, 4.0)
plt.xticks([])
plt.yticks(np.linspace(2.0, 4.0, 5))
plt.savefig("1D_Abs_Spect_ky=0_θ=" + str(θ)  + ".pdf", bbox_inches='tight')
plt.close()
# For Kx = 0 
initkx0 = glob.glob(os.path.join(base_path, "initk_*"))
initkx0.sort(key=extract_number)
I2Dkx = np.zeros((len(ω), len(initkx0)))
print(initkx0)
for i, c in enumerate(initkx0):
    InitPht = par.Nsite*par.Nlayer 
    fold_kx = "init_" + str(InitPht)
    file2 = os.path.join(c, fold_kx, "splitop_tilted-EP1D.txt")
    if not os.path.exists(file2):
        print(f"File {file2} does not exist, skipping.")
    else:
        t, ct = Ct(file2)
        ω, I = Iw2(t, ct)
        I2Dkx[:,i] = I
        print(f"Processing {i} with shape {I.shape}")
file_name = '1D_Abs_Spect_kx=0_θ='+ str(θ)  + '.txt'
np.savetxt(file_name, I2Dkx)

M2 = np.zeros((I2Dkx[::-1,:].shape[0], I2Dkx[::-1,:].shape[1]*2))
M2[:, :I2Dkx[::-1,:].shape[1]] = I2Dkx[::-1,::-1]
M2[:, I2Dkx[::-1,:].shape[1]:] = I2Dkx[::-1,:]
plt.figure(figsize=(5, 5/ 1.034))
plt.imshow(M2, aspect='auto' , extent=[ -maxkx, maxkx, 0.0, 5.0],cmap=cmap, interpolation='gaussian')
plt.clim(0, 250)
plt.xlabel('$k_x$ (1/a × $10^3$)')
plt.ylim(2.0, 4.0)
plt.tight_layout()
plt.yticks(np.linspace(2.0, 4.0, 5))
plt.ylabel('Energy (eV)')
plt.savefig('1D_Abs_Spect_kx=0_θ=' + str(θ)  + '.pdf', bbox_inches='tight')   
plt.close()
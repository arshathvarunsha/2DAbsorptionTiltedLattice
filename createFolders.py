#!/sw/eb/sw/Anaconda3/2024.02-1/bin/python
##NECESSARY JOB SPECIFICATIONS
#SBATCH --partition=shared
#SBATCH --job-name=Run_NAMD         #Set the job name to "JobExample1"
#SBATCH --time=1:00:00          #Set the wall clock limit to 1 hour
#SBATCH --ntasks=1               #Request 1 task
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M              
#SBATCH --output=output/out.%j  
import os
import sys
import numpy as np
import glob
 
sys.path.append(os.popen("pwd").read().replace("\n","")+"/Model")
exec(f"import EP1D as model")
par = model.parameters()

frac = 0.01
npoints = 40
df = frac / npoints

NMAX = int(par.N_Pht * frac)
kx = np.arange(0, NMAX, int(NMAX/npoints))  # kx values for the x-axis

 
 
base_dir = os.getcwd()
model_dir = os.path.join(base_dir, "Model")
method_dir = os.path.join(base_dir, "Method")

outputname = "splitop_tilted-EP1D.txt"

#  = kx[int(sys.argv[1])] # Get the index from command line argument or default to 0
dirs = open("commands.in", "w")
for yk in kx:
    i = yk
    print(f"Creating initKy={i}")

    for k in kx:
        initState = int(par.Nsite * par.Nlayer + int(k))    
        fold = os.path.join(f"initk_{i}", f"init_{initState}")
        

        fout = glob.glob(fold+f"/{outputname}")

        if len(fout) > 0:
            print(f"Simulation for initState={initState}, initKy={i} already completed in folder {fold}. Skipping...")
        else:
            os.makedirs(fold, exist_ok=True)
            #os.system(f"cp -r {method_dir} {fold}")      
            os.system(f"cp -r {model_dir} {fold}")
            # os.system(f"cp {base_dir}/run.py {fold}")
            os.system(f"cp {base_dir}/serial.py {fold}")
            #os.system(f"cp {base_dir}/__init__.py {fold}")
            os.system(f"cp {base_dir}/input.txt {fold}")
            inputfile = open(os.path.join(model_dir, "EP1D.py"), "r").readlines()
            inputfile = [line if "    initState = " not in line else f"    initState  = {initState}  # Updated γ value\n" for line in inputfile]
            inputfile = [line if "    initKy = " not in line else f"    initKy = {i}  # Updated γ value\n" for line in inputfile]            
            with open(f"{fold}/Model/EP1D.py", "w") as fob:
                fob.writelines(inputfile)
            os.chdir(fold)
            # os.system("python run.py")
            os.chdir("../../")
            dirs.write(f"cd {fold}; python serial.py input.txt output; cd ../../\n")

dirs.close()
print(f"Directories for simulations written to commands.in")
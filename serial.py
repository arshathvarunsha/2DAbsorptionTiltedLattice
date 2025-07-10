#!/sw/eb/sw/Anaconda3/2024.02-1/bin/python

##NECESSARY JOB SPECIFICATIONS
#SBATCH --partition=shared
#SBATCH --job-name=NAMD         #Set the job name to "JobExample1"
#SBATCH --time=00:10:00          #Set the wall clock limit to 10 minutes
#SBATCH --ntasks=1               #Request 1 task
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M              
#SBATCH --output=output/out.%j            


import os
import sys

import numpy as np
mainfold = os.popen("pwd").read().replace("\n","")
# go back two directories to find Method and Model
sys.path.append(os.path.join(mainfold, "../../Method"))
sys.path.append(os.popen("pwd").read().replace("\n","")+"/Model")

#-------------------------
try:
    inputtxt = open(sys.argv[1], 'r').readlines()
    print(f"Reading {sys.argv[1]}")
except:
    print("Reading input.txt")
    inputtxt = open('input.txt', 'r').readlines()


def getInput(input,key):
    try:
        txt = [i for i in input if i.find(key)!=-1][0].split("=")[1].split("#", 1)[0].replace("\n","")
    except:
        txt = ""
    return txt.replace(" ","")

model_ =  getInput(inputtxt,"Model")
print(f"Model: {model_}")
method_ = getInput(inputtxt,"Method").split("-")
exec(f"import {model_} as model")
exec(f"import {method_[0]} as method")
try:
    stype = method_[1]
except:
    stype = "_"
#-------------------------
import time

import numpy as np

t0 = time.time()

try:
    fold = sys.argv[2]
except:
    fold = "./output"


os.system(f"mkdir -p {fold}")
ID = ''
try :
    ID = sys.argv[3]
    ID = "-" + ID
except:
    pass

t1 = time.time()

NTraj = model.parameters.NTraj
NStates = model.parameters.NStates

#------ Arguments------------------
par = model.parameters() 
print("par was read, γ=", par.γ)
par.ID     = np.random.randint(0,100)
par.SEED   = np.random.randint(0,100000000)
    
#---- methods in model ------
#par.dHel = model.dHel
par.dHel0 = model.dHel0
par.initR = model.initR
par.HelR   = model.HelR
par.stype = stype

if method_[0]=="nrpmd":
    par.initHel0 = model.initHel0
    

#---- overriden parameters ------

parameters = [i for i in inputtxt if i.split("#")[0].split("=")[0].find("$") !=- 1]
for p in parameters:
    exec(f"par.{p.split('=')[0].split('$')[1]} = {p.split('=')[1].split('#')[0]}")
    print(f"Overriding parameters: {p.split('=')[0].split('$')[1]} = {p.split('=')[1].split('#')[0]}")
#--------------------------------

#------------------- run --------------- 
rho_sum  = method.runTraj(par)
#--------------------------------------- 

try:
    PiiFilename =  f"{fold}/{method_[0]}-{method_[1]}-{model_}{ID}_{parameters.initState}.txt"
    psiFilename =  f"{fold}/psi-{method_[0]}-{method_[1]}-{model_}{ID}.npy"  
    #PiiFile = open(PiiFilename,"w+") 
except:
    PiiFilename = f"{fold}/{method_[0]}-{model_}{ID}.txt"
    psiFilename = f"{fold}/psi-{method_[0]}-{model_}{ID}.npy"
    #

NTraj = par.NTraj

#times = np.arange(rho_sum.shape[-1]) * par.nskip * par.dtN, par.nskip * par.dtN


rho_sum[:] = rho_sum[:]/ NTraj 
np.savetxt(PiiFilename, np.c_[rho_sum.real, rho_sum.imag])


# def MSD(NStates, a, rhodiag):
#     Matrix_r = np.arange(L)*a
#     return sum(rhodiag*(np.arange(L)**2)*(a**2) - (rhodiag*np.arange(L)*a)**2)


t2 = time.time()-t1
print(f"Total Time: {t2}")
print(f"Time per trajectory: {t2/NTraj}")
 
##  2D Absorption Spectra Workflow

This repository provides a framework for computing and plotting 2D absorption spectra using a parameterized model.  

---

## 🛠️ Step 1: Configure Simulation Parameters

Before running the job, edit the following file to set the model parameters:  
```Model/EP1D.py ```
,update the variables as needed to match your simulation setup.  

**Trial:**  
For testing, the parameters are set to a non-tilted system of 8000 excitons.  

---

## 🚀 Step 2: Run the Simulation

This code is written to run in a parallel manner using the **TAMULauncher** cluster system.  

- For detailed information, check the [HPRC TAMU TAMULauncher guide](https://hprc.tamu.edu/kb/Software/tamulauncher/).  
- For users outside TAMU, the process can be replaced with **SLURM array job submission**, which is equivalent. See [reference](https://stackoverflow.com/questions/71969482/running-parallel-jobs-in-slurm).  

### If you are using TAMULauncher
As an initial step before running jobs, in the CWD run:  
```bash tamulauncher --remove-logs commands.in ```
This will clear the cache memory associated with previous simulations. Otherwise, the new job will not run.  

To generate raw data for the absorption spectra, submit the main job script. If you are using SLURM:  

```bash sh start.sh ```
### (A) This script will:


(1) Create subfolders for each $k_y$ value.


(2) Generate the commands.in file, which describes jobs to run parellely.


(3) Launch simulations over different ($k_x$, $k_y$) values.


### (B) Submit supporting scripts:

(1) tamulauncher.sh → Submits all parallel runs.


(2) clean.sh → Cleans temporary files, prevents exceeding file limits, and plots data once simulations complete.
    Allows on-the-fly analysis while jobs are running.

    
The above steps will produce the raw absorption spectra data and intermediate output files.


## Note: How Modify Energy for 2D-cut


Open clean.sh and set your desired energy range for plotting by changing *E = *.
This can be adjusted at any time, even after the simulations are finished.

## 🧠 Additional Notes
You can modify the logic and plotting styles by editing the following files:

Method/


Model/


Plotting.py


Each script is well-commented to facilitate customization and understanding.


You will able to reproduce this result if you run this code with the test parameters given in the model, this is the 2D absorbtion spectra for a system with out any tilt at E= 3 eV .


![image](https://github.com/user-attachments/assets/10cee291-8220-4b73-852f-1f859a33f45e)



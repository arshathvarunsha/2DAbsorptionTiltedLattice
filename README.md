2D Absorption Spectra Workflow
This repository provides a framework for computing and plotting 2D absorption spectra using a parameterized model.

🛠️ Step 1: Configure Simulation Parameters
Before running the job, edit the following file to set the model parameters:
Model/EP1D.py
Update the variables as needed to match your simulation setup. 
Trial
For testing the parameters are set to a non-tilted system of 8000 excitons.

🚀 Step 2: Run the Simulation
To generate raw data for the absorption spectra, submit the main job script. If you're using SLURM:
sbatch run_abs_2D.py
This script will run the required simulations and output a number of intermediate files.

🧹 Step 3: Clean and Average the Output
After the simulation is complete, you can clean and average the results to prepare for plotting.
3.1 Modify Energy Range
Set the desired energy range for plotting within the clean.sh script. You can adjust this at any time.
3.2 Run the Cleanup Script
Submit the script using SLURM:
sbatch clean.sh
This will remove unnecessary files and produce an averaged dataset suitable for visualization.

🧠 Additional Notes
You can modify the logic and plotting styles by editing the following files:
Method/
Model/
Plotting.py
Each script is well-commented to facilitate customization and understanding.
![image](https://github.com/user-attachments/assets/10cee291-8220-4b73-852f-1f859a33f45e)

##📬 Questions or Contributions
For any issues, suggestions, or contributions, feel free to open an issue or submit a pull request.
Happy Simulating!

### LMME: Coherent Frequency Conversion

The light–matter moiré effect (LMME) can be used for **coherent frequency conversion**.  
You can initialize the system in a photonic state <span>A&#x0302;<sub>k'<sub>||</sub></sub> &nbsp;|&nbsp;<span>0&#x0305;</span>&#x27E9;</span></span> and obtain an output at
<span>A&#x0302;<sub>k'<sub>||</sub> + 2∆k</sub> &nbsp;|&nbsp;<span>0&#x0305;</span>&#x27E9;</span></span>.

We prepare an initial superposition of two single-photon modes with a relative phase ϕ as shown below:

<img width="347" height="56" alt="image" src="https://github.com/user-attachments/assets/0fb4db73-06bf-44a7-a7a9-9120cb4b1826" />

Here, ϕ is the relative phase between the two modes.

---

#### Selecting <span> A&#x0302;<sup>&dagger;</sup><sub>k&#x2225;&#x2032;</sub>  &nbsp;|&nbsp;<span>0&#x0305;</span>&#x27E9;</span>


Edit the file:

 ```Model/EP1D.py ``` 

 
Set the desired initial momentum state via `initKstate`, e.g.:

<img width="425" height="57" alt="image" src="https://github.com/user-attachments/assets/f6117632-aa32-4f75-85e7-00aad3fb5c61" />

---

#### Setting the phase difference

Adjust the phase parameter (the relative phase \(\phi\)) as shown:

<img width="853" height="47" alt="image" src="https://github.com/user-attachments/assets/04328af9-a048-425c-ac30-3481bc60e9bc" />

---
### Running the Simulation
In the current working directory, run this command in terminal:
 ```python serial.py input.txt``` or if you are using slurm  ```sbatch serial.py input.txt```
#### Outputs

Running the simulation produces **four outputs** and **two figures** corresponding to the **input** and **output** superpositions:

<img width="558" height="284" alt="image" src="https://github.com/user-attachments/assets/d1a09cb1-d09e-4906-877c-097836af19ad" />

<img width="557" height="284" alt="image" src="https://github.com/user-attachments/assets/5215ea97-9d24-497e-a2eb-8e2e13baa0ae" />

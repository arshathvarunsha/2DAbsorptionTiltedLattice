##  Tilted Material in an Optical Cavity: Light-Matter Moiré Effect and Coherent Frequency Conversion

Overview


This repository contains the simulation framework and analysis tools for studying the Light–Matter Moiré Effect (LMME) that emerges when a tilted 2D material is embedded inside a planar optical cavity.

Unlike conventional moiré materials formed by twisting layered crystals, here the moiré-like modulation arises purely from geometric tilt, leading to displaced polariton dispersions, emergent flat bands, and robust coherent frequency conversion.

![image](https://github.com/user-attachments/files/21908472/Fig0.pdf)

 *Comparison of twisted graphene and tilted material in an optical cavity} (a) Illustration of a graphene bilayer offset by a twist angle. (b) Illustration of a tilted material within a Fabry-Perot cavity coupling to cavity radiation. (c) Schematic band structure of a graphene bilayer without (left) and with (right) a twist of θ. (d) Schematic band structure of a tilted material without (left) and with (right) a tilt of  θ.*

## Theoretical Frame Work and General Work Flow
The system is modeled by a light–matter Hamiltonian that captures the coupling between cavity photons and excitons in a tilted 2D lattice. The tilt introduces a spatial modulation in the exciton–photon coupling, giving rise to the LMME.
The Hamiltonian for a fixed  momentum vector component along $\vec{y}, i.e, $$k_y$ is expressed as:
<img width="901" height="176" alt="image" src="https://github.com/user-attachments/assets/292f4e49-b9a6-4465-8c2d-8aafbae3f0a0" />

This decomposition allows us to perform independent 2D simulations for each 
$k_y$. The full 3D dynamics of the system can then be reconstructed by combining the results across all the $k_y$ values.This modular approach simplifies computations and allows efficient parallelization.
    
## Repository Structure and Usage
(1) 2D_ABS_SPECTRA/        --------------------------------------- # Scripts for absorption spectra calculations

(2) COHERENT_FREQ_CONV/    ---------------------------------------# Simulations of coherent frequency conversion

(3) FIGURES/               ---------------------------------------# Scripts & data for reproducing paper figures

### References
A. Manjalingal, S. R. Koshkaki, L. Blackham, A. Mandal,
Tilted Material in an Optical Cavity: Light–Matter Moiré Effect and Coherent Frequency Conversion (2025), [arXiv:2508.11237](https://arxiv.org/abs/2508.11237).
### 📬 Questions or Contributions
For any issues, suggestions, or contributions, feel free to open an issue or submit a pull request.
contact: 


✉️ arkajit@tamu.edu

✉️ avmanjalingal@tamu.edu

✉️ rahmanian@tamu.edu


<div style="display: flex; justify-content: space-between;">
  <img src="https://github.com/user-attachments/assets/407ab9e2-d7ab-427e-a645-9e1afe513b56" width="300" />
  <img src="https://github.com/user-attachments/assets/6d026d1f-5a9b-4890-8874-7485d2ead16a" width="300" />
</div>


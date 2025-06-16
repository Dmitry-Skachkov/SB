
# Preparation calculation

There is necessary to do several DFT (or any first-principles) calculations and prepare the following input files:

* input file *dos_bulk_.dat*. DFT calculation for the bulk of semiconductor in order to prepare density of states (DOS) file, energies should be shifted with respect to valence band maximum (VBM)
* to calculate projected density of states (PDOS) separated by k-points for interface system containing the metal surface and six layers of the semiconductor, and prepare the files *k_pdos_int_0.dat*, *k_pdos_int_3.dat*, and *k_pdos_int_4.dat* for surface layer and for third and fourth layers of the semiconductor. Energies should be shifted with respect to VBM of the bulk     
* input file *cbs.dat* containing the decaying rates of complex band structure (CBS) separated by light and heavy holes bands. k-point mesh for CBS calculation should corresponds to interface calculation
* input file *polarization.dat* containing the polarization of the bulk semiconductor with respect to electric field. It may contain only one number for static dielectric constant (see [Example_1](https://github.com/Dmitry-Skachkov/SB/tree/main/Examples/Example_1))
* input file *k_mesh.dat* containing k-mesh points from interface calculation in crystal representation
* input file *input.dat* containing the following information:
  - V_D0 in a.u. (output from QE),  
  - Lz (A)
  - Lz_int (A)
  - CNL (eV)
  - z3 (A)
  - z4 (A)
  - alat (au)
  - reciprocal vector b1
  - b2
  - b3
  - cz

# Running the SB program

In the folder with the input files run the program:

> shottky 300. 1.42 inf 0. > output

where 300. is the temperature, 1.42 is the Fermi level of the system, 'inf' corresponds to semiinfinite size of the semiconductor, and 0. corresponds to the surface charge on the gating electrode in cm-2

The *output* file contains the information about the Schottky barrier parameters: Schottky barrier height (SBH), surface filling level &Delta;E, charge density on the surface of the interface, the electric field close to the surface of the semiconductor.

There are two additional output files:
- *po_.dat* containing the charge densities for electrons, holes, MIGS and total (in cm-3)
- *elpot_.dat* containing the electrostatic potential (in V) and electric field (in V/A)



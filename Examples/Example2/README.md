# Example how to calculate Shottky Barrier without CBS effects


## Run example:
> schottky 300. 1.42 inf 0. > output.txt

here

300 - temperature (K)

1.42 - Fermi level (eV)

inf - length of the semiinfinite semiconductor (A)

0 - gating charge density (cm-2). If length of the semiconductor is infinite, then gating charge density will set to zero.

## Input files:

input.dat - input parameters of the system

dos_bulk_.dat - DOS of the bulk semiconductor

DOStot.dat - DOS of the surface

k_mesh.dat - mesh of k-points

kpdos_int_3.dat - PDOS separated by k-points for 3d layer of the interface

polarization.dat - contains dielectric constant for the bulk semiconductor

## Compare with the original output
> diff output.txt output_orig_noCBS.txt






# GaAs-Graphene

![GitHub Logo](https://github.com/Dmitry-Skachkov/SB/blob/main/Docs/SB_geometry_1.jpg)

## To run the example:
> schottky 300. 1.42 inf 0. > output

here

300 - temperature (K)

1.42 - Fermi level (eV)

inf - length of the semiconductor (inf or length in A, inf means semiinfinite)

0 - gating charge density (cm-2). If length of the semiconductor is infinite, then gating charge density will set to zero.

## Input files:

input.dat - input parameters of the system

dos_bulk_.dat - DOS of the bulk semiconductor

DOStot.dat - DOS of the surface containing one layer of metal and first layer of the semiconductor

k_mesh.dat - mesh of k-points

kpdos_int_3.dat - PDOS separated by k-points for third layer of the semiconductor in the interface

kpdos_int_4.dat - PDOS separated by k-points for fourth layer of the semiconductor in the interface

cbs.dat - file containing CBS with light and heavy complex bands. Folder [QE_CBS](QE_CBS) contains the input files for calculation of CBS with QE

polarization.dat - this file contains the set of calculated values of polarization with respect to electric field P(E). This file may contain only one number, the static dielectric constant, see [Example_1](https://github.com/Dmitry-Skachkov/SB/tree/main/Examples/Example_1).   

## Compare the results with the original:
> diff output output_orig


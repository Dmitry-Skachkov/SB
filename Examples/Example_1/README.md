
# Example how to calculate Shottky Barrier without CBS effects

![GitHub Logo](https://github.com/Dmitry-Skachkov/SB/blob/main/Docs/SB_geometry_1.jpg)

## Run example:
> schottky 300. 1.42 inf 0.0 > output.txt

here

300 - temperature (K)

1.42 - Fermi level (eV)

inf - length of the semiconductor (inf or length in A). inf means semiinfinite [0,∞)

0.0 - gating charge density (cm-2). If length of the semiconductor is infinite, then gating charge density will set to zero.

## Input files:

[input.dat](input.dat) - input parameters of the system

[dos_bulk_.dat](dos_bulk_.dat) - DOS of the bulk semiconductor

[DOStot.dat](DOStot.dat) - DOS of the surface containing one layer of metal and first layer of the semiconductor

[polarization.dat](polarization.dat) - contains dielectric constant for the bulk semiconductor or P(E). In this example, this file contains only one number for static dielectric constant. In [Example_2/polarization.dat](https://github.com/Dmitry-Skachkov/SB/tree/main/Examples/Example_2/polarization.dat), this file contains the set of calculated values of polarization with respect to electric field P(E).   

## Compare with the original output
> diff output.txt output_orig_noCBS.txt






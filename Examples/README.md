# Examples how to run the SB program:


[GaAs-Gr](GaAs-Gr) - example for GaAs-Graphene system

[Example2](Example2) - the same system but without including complex band structure (CBS) effects 



# How to prepare for the calculation:

1. Calculate DOS of the bulk
2. Calculate PDOS with k-point separated of the interface. The interface is the contact of few layers of metal with few layers of the semiconductor H-terminated on the surface exposed to vacuum.
3. If it is necessary to include complex band structure (CBS) effects, then to calculate CBS with QE and prepare cbs.dat file with two columns for light and heavy complex bands. See [cbs.dat](GaAs-Gr/cbs.dat) example. 



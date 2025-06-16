# Examples how to run the SB program:

[Example_1](Example_1) - GaAs-Graphene system, calculating Schottky barrier parameters without including complex band structure (CBS) effects 

[Example_2](Example_2) - GaAs-Graphene system, including CBS effects

[Example_3](Example_3) - calculate Schottky barrier with gate electrode



# How to prepare for the calculation:

1. Calculate DOS of the bulk
2. Calculate PDOS for each layer of the interface. The interface is the contact of few layers of metal with few layers of the semiconductor H-terminated on the surface exposed to vacuum.
3. If it is necessary to include CBS effects, then to calculate CBS with QE and prepare cbs.dat file with two columns for light and heavy complex bands. See [cbs.dat](Example_2/cbs.dat) example. 

For more details see [here](https://github.com/Dmitry-Skachkov/SB/blob/main/Docs/USAGE.md)

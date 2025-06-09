# Example how to calculate Shottky Barrier without CBS effects


## Run example:
> schottky 300. 1.42 1e4 0.0 > output.txt

here

300 - temperature (K)

1.42 - Fermi level (eV)

1e4 - length of the semiconductor (A), 1 micron

0.0 - gating charge density (cm-2)

Run several calculations for gating surface charges 10^6, 10^7 cm^-2
```
Sg (cm-2)     DLW (A)
0             1322
10^6          1435
10^7          4394
```




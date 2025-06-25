# Shottky Barrier for finite length semiconductor

![GitHub Logo](https://github.com/Dmitry-Skachkov/SB/blob/main/Docs/SB_geometry3.jpg)

## Run example:
> schottky 300. 1.42 1e4 1e7 > output.txt

here

300 - temperature (K)

1.42 - Fermi level (eV)

1e4 - length of the semiconductor (10<sup>4</sup> A = 1 micron)

1e7 - gating charge density (10<sup>7</sup> cm<sup>-2</sup>)

Run several calculations for gating surface charges 0, 10<sup>6</sup>, 10<sup>7</sup>, -10<sup>7</sup> cm<sup>-2</sup>
```
Sigma_g (cm-2)   DLW (A)
0                1367
10^6             1484
10^7             4394
-10^7             984
```




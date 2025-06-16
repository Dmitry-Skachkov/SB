# First-Principles Method for Schottky Barrier

* [Abstract](README.md#abstract)
* [Usage](README.md#usage)
* [Authors](README.md#authors)
* [License](README.md#license)
* [Contributing](README.md#contributing)
* [How to cite](README.md#how-to-cite)
* [Versions](VERSION.md)
* [Features](FEATURES.md)

# Abstract

We develop a first-principles method for calculating Schottky barrier parameters. The Poisson equation is solved self-consistently with the electrostatic charge density over the entire barrier using the density functional theory (DFT) electronic structure converged locally, allowing computation of a Schottky barrier entirely from DFT involving thousands of atomic layers in the semiconductor. The induced charge in the bulk consists of conduction and valence band charges from doping and band bending, as well as charge from the evanescent states in the gap of the semiconductor. The Schottky barrier height is determined when the induced charge density and the induced electrostatic potential reach self-consistency. The Schottky barrier height, width, along with depletion and inversion layers obtained self-consistently as functions of temperature and bulk doping.

![GitHub Logo](https://github.com/Dmitry-Skachkov/SchottkyBarrier/blob/main/Docs/logo.jpg)

# Usage

Please see [Docs](Docs) folder for instructions how to install the code and [Examples](Examples) folder for some sample calculations. 

# Authors

[Dmitry Skachkov](mailto:dmitry.skachkov@dsedu.org)

[Xiao-Guang Znang](mailto:xgz@ufl.edu)

[Hai-Ping Cheng](mailto:ha.cheng@northeastern.edu)

# License

[AFL 3.0 License](https://github.com/Dmitry-Skachkov/SB/blob/main/LICENSE.md) 

# Contributing

For contributing please contact Dr. Hai-Ping Cheng at m2qm.efrc[at]phys.ufl.edu

# How to cite

If you use the SB package in your research, please cite the following publications:

* Dmitry Skachkov, Shuang-Long Liu, Yan Wang, Xiao-Guang Zhang, Hai-Ping Cheng  
"First-Principles Theory for Schottky Barrier Physics" Phys. Rev. B, 104, 045429 (2021) [DOI 10.1103/PhysRevB.104.045429](https://doi.org/10.1103/PhysRevB.104.045429)  
[![arXiv](https://img.shields.io/badge/arXiv-2001.00710-b31b1b.svg?style=plastic)](https://arxiv.org/abs/2001.00710)
* D. Skachkov, X.-G. Zhang, H.-P. Cheng (2021) SB: GitHub 

See [BibTex](BibTex.md) entry for the Github repository and the publications
  
## For information about versions and updates see [VERSION.md](VERSION.md) 

## For features of the SB program see [FEATURES.md](FEATURES.md) 






# First-Principles Method for Schottky Barrier

We develop a first-principles theory for Schottky barrier physics. The Poisson equation is solved self-consistently with the electrostatic charge density over the entire barrier using the density functional theory (DFT) electronic structure converged locally, allowing computation of a Schottky barrier entirely from DFT involving thousands of atomic layers in the semiconductor. The induced charge in the bulk consists of conduction and valence band charges from doping and band bending, as well as charge from the evanescent states in the gap of the semiconductor. The Schottky barrier height is determined when the induced charge density and the induced electrostatic potential reach self-consistency. Tests on the GaAs – graphene and Si/Al heterostructures yield Schottky barrier height, width, along with depletion and inversion layers obtained self-consistently as functions of temperature and bulk doping.

![GitHub Logo](https://github.com/Dmitry-Skachkov/SchottkyBarrier/blob/main/Docs/logo.jpg)


# Installation

Copy sb_v1.0.tar file with the archive into the selected folder. Untar the archive:

> tar -xvf sb_v1.0.tar

Load the appropriate modules for compilers. Compile the program:

> make

# Usage

Please see the file ![USAGE.md](https://github.com/Dmitry-Skachkov/SchottkyBarrier/blob/main/Docs/USAGE.md) in the docs folder for instructions and the examples folder for some sample calculations.

# License

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# Contributing

Please contact at m2qm.efrc@phys.ufl.edu

# How to cite

If you use the FPSchottkyBarrier package in your research, please cite the following papers:

"First-Principles Theory for Schottky Barrier Physics"
Dmitry Skachkov, Shuang-Long Liu, Yan Wang, Xiao-Guang Zhang, Hai-Ping Cheng
https://arxiv.org/abs/2001.00710


<   @misc{VASPsol-Software,
     title        = {VASPsol: Implicit solvation and electrolyte model for density-functional theory},
     author       = {K. Mathew and V. S. Chaitanya Kolluru and R. G. Hennig},
     year         = 2018,
     publisher    = {GitHub},
     journal      = {GitHub repository},
     howpublished = {\url{https://github.com/henniggroup/VASPsol}},
     url          = {https://github.com/henniggroup/VASPsol},
     doi          = {10.5281/zenodo.2555053}
   }
   
   @article{VASPsol2014-Dielectric,
     title        = {Implicit solvation model for density-functional study of nanocrystal surfaces
                     and reaction pathways.},
     author       = {K. Mathew and R. Sundararaman and K. Letchworth-Weaver and T. A. Arias and
                     R. G. Hennig},
     year         = 2014,
     journal      = {J. Chem. Phys.},
     volume       = 140,
     pages        = {084106},
     doi          = {10.1063/1.4865107}
   }
   
   @article{VASPsol2019-Electrolyte,
     title        = {Implicit self-consistent electrolyte model in plane-wave density-functional theory.},
     author       = {K. Mathew and V. S. C. Kolluru and S. Mula and S. N. Steinmann and R. G. Hennig},
     year         = 2019,
     journal      = {J. Chem. Phys.},
     volume       = 151,
     pages        = {234101},
     doi          = {10.1063/1.5132354}
   }
>





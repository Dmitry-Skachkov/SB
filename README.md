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

# Authors

Dmitry Skachkov

Xiao-Guang Znang

Hai-Ping Cheng

# License

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# Contributing

Please contact at m2qm.efrc[at]phys.ufl.edu

# How to cite

If you use the SB package in your research, please cite the following papers:

Dmitry Skachkov, Shuang-Long Liu, Yan Wang, Xiao-Guang Zhang, Hai-Ping Cheng  
"First-Principles Theory for Schottky Barrier Physics"  
https://arxiv.org/abs/2001.00710

BibTex entry for the Github repository and the publications:

    @misc{SB-Software,
     title        = {SB: First-Principles Method for Schottky Barrier},
     author       = {D. Skachkov and Xiao-Guang Zhang and Hai-Ping Cheng},
     year         = 2021,
     publisher    = {GitHub},
     journal      = {GitHub repository},
     howpublished = {\url{https://github.com/m2qm-efrc/SB}},
     url          = {https://github.com/m2qm-efrc/SB},
     doi          = {}
    }
   
    @article{SB-paper,
     title        = {First-Principles Theory for Schottky Barrier Physics},
     author       = {D. Skachkov and S.-L. Liu and Y. Wang and X.-G. Zhang and H.-P. Cheng},
     year         = 2021,
     journal      = {ArXiv},
     volume       = 2001,
     pages        = {00710},
     doi          = {}
    }
   







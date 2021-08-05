
# Installation

Copy <arxive>.tar file with the archive into the selected folder. Untar the archive:

> tar -xvf <arxive>.tar

Load the appropriate modules for Fortran compilers and compile the program:

> make
  
You can modify makefile in order to use the selected compiler. Currently GNU/4.8.5, Intel/19.0.5, and PGI NVIDIA/20.1 compilers are tested.

# Usage

Please see file [USAGE.md](USAGE.md) for instructions and [Examples](../Examples) folder for some sample calculations.

# Further development

Please see file [Developer_documentation.pdf](https://github.com/Dmitry-Skachkov/SB/blob/main/Docs/Devepoler_documentation.pdf) for the detailed scheme of the SB program. For any questions related to the code please contact Dr. Dmitry Skachkov at d.skachkov[_at_]ufl.edu. 

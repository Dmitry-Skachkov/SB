
# Installation

Copy sb.tar file with the archive into the selected folder. Untar the archive:

> tar -xvf sb.tar

Load the appropriate module for Fortran compiler. You can modify makefile in order to use the selected compiler. Currently GNU/4.8.5, Intel/19.0.5, and PGI NVIDIA/20.1 compilers are tested. Compile the program:

> make
  

# Usage

Please see file [USAGE.md](USAGE.md) for instructions and [Examples](../Examples) folder for some sample calculations.

# Further development

Please see file [Developer_documentation.pdf](https://github.com/Dmitry-Skachkov/SB/blob/main/Docs/Devepoler_documentation.pdf) for the detailed scheme of the SB program. Please also see [COMMENTS](COMMENTS.md) for description of the variables in the program. For questions related to the code please contact [Dr. Dmitry Skachkov](<mailto:dmitry.skachkov@DSedu.org?subject=SB code on GitHub>). 

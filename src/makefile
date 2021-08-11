default: schottky

include gfortran.make                # use GNU compiler
#include pgi.make                    # use PGI NVIDIA compiler
#include intel.make                  # use Intel compiler

OBJS = obj/SBMathLibrary.o \
       obj/SBParameters.o  \
       obj/SBSpline_functions.o \
       obj/SBCalc_i.o \
       obj/SBCalc_z.o \
       obj/SBPlot.o \
       obj/SB.o

schottky: $(OBJS)
	$(FC) $(FFLAGS) -o ~/bin/schottky $(OBJS)

clean:
	rm -rf ~/bin/schottky obj/*.o obj/*.mod

SBSpline_functions.o: SBParameters.o  SBMathLibrary.o
SBCalc_i.o:           SBParameters.o  SBMathLibrary.o   SBSpline_functions.o
SBCalc_z.o:           SBParameters.o  SBCalc_i.o        SBSpline_functions.o
SBPlot.o:             SBParameters.o  SBCalc_i.o        SBSpline_functions.o 
SB.o:                 SBParameters.o  SBCalc_i.o        SBSpline_functions.o  SBPlot.o  SBCalc_z.o




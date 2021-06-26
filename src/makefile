default: schottky

obj/SBMathLibrary.o: SBMathLibrary.f90
	@mkdir -p obj
	mpif90 -c SBMathLibrary.f90      -o obj/SBMathLibrary.o         -module obj
obj/SBParameters.o: SBParameters.f90
	mpif90 -c SBParameters.f90       -o obj/SBParameters.o          -module obj
obj/SBSpline_functions.o: SBSpline_functions.f90
	mpif90 -c SBSpline_functions.f90 -o obj/SBSpline_functions.o    -module obj
obj/SBCalc_i.o: SBCalc_i.f90
	mpif90 -c SBCalc_i.f90           -o obj/SBCalc_i.o              -module obj
obj/SBCalc_z.o: SBCalc_z.f90
	mpif90 -c SBCalc_z.f90           -o obj/SBCalc_z.o              -module obj
obj/SBPlot.o: SBPlot.f90
	mpif90 -c SBPlot.f90             -o obj/SBPlot.o                -module obj
obj/SB.o: SB.f90
	mpif90 -c SB.f90                 -o obj/SB.o                    -module obj

OBJS = obj/SBMathLibrary.o \
       obj/SBParameters.o  \
       obj/SBSpline_functions.o \
       obj/SBCalc_i.o \
       obj/SBCalc_z.o \
       obj/SBPlot.o \
       obj/SB.o

schottky: $(OBJS)
	mpif90 -o ~/owens/bin/schottky_v1 $(OBJS)

clean:
	rm -rf ~/owens/bin/schottky_v1 obj/*.o obj/*.mod

SBSpline_functions.o: SBParameters.o SBMathLibrary.o
SBCalc_i.o:           SBParameters.o SBMathLibrary.o   SBSpline_functions.o
SBCalc_z.o:           SBParameters.o SBCalc_i.o        SBSpline_functions.o
SBPlot.o:             SBParameters.o SBCalc_i.o        SBSpline_functions.o 
SB.o:                 SBParameters.o SBCalc_i.o        SBSpline_functions.o SBPlot.o SBCalc_z.o



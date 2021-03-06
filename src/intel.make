FC = ifort
FFLAGS = 

obj/SBMathLibrary.o: SBMathLibrary.f90
	@mkdir -p obj
	$(FC) $(FFLAGS) -c SBMathLibrary.f90      -o obj/SBMathLibrary.o         -module obj
obj/SBParameters.o: SBParameters.f90
	$(FC) $(FFLAGS) -c SBParameters.f90       -o obj/SBParameters.o          -module obj
obj/SBSpline_functions.o: SBSpline_functions.f90
	$(FC) $(FFLAGS) -c SBSpline_functions.f90 -o obj/SBSpline_functions.o    -module obj
obj/SBCalc_i.o: SBCalc_i.f90
	$(FC) $(FFLAGS) -c SBCalc_i.f90           -o obj/SBCalc_i.o              -module obj
obj/SBCalc_z.o: SBCalc_z.f90
	$(FC) $(FFLAGS) -c SBCalc_z.f90           -o obj/SBCalc_z.o              -module obj
obj/SBPlot.o: SBPlot.f90
	$(FC) $(FFLAGS) -c SBPlot.f90             -o obj/SBPlot.o                -module obj
obj/SB.o: SB.f90
	$(FC) $(FFLAGS) -c SB.f90                 -o obj/SB.o                    -module obj


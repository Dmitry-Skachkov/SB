FC = gfortran
FFLAGS = -ffree-line-length-512

obj/SBMathLibrary.o: SBMathLibrary.f90
	@mkdir -p obj
	$(FC) $(FFLAGS) -c SBMathLibrary.f90      -o obj/SBMathLibrary.o         -J./obj
obj/SBParameters.o: SBParameters.f90
	$(FC) $(FFLAGS) -c SBParameters.f90       -o obj/SBParameters.o          -J./obj
obj/SBSpline_functions.o: SBSpline_functions.f90
	$(FC) $(FFLAGS) -c SBSpline_functions.f90 -o obj/SBSpline_functions.o    -J./obj
obj/SBCalc_i.o: SBCalc_i.f90
	$(FC) $(FFLAGS) -c SBCalc_i.f90           -o obj/SBCalc_i.o              -J./obj
obj/SBCalc_z.o: SBCalc_z.f90
	$(FC) $(FFLAGS) -c SBCalc_z.f90           -o obj/SBCalc_z.o              -J./obj
obj/SBPlot.o: SBPlot.f90
	$(FC) $(FFLAGS) -c SBPlot.f90             -o obj/SBPlot.o                -J./obj
obj/SB.o: SB.f90
	$(FC) $(FFLAGS) -c SB.f90                 -o obj/SB.o                    -J./obj





       Program SB                                                   ! calculate Schottky barrier parameters
        use SBParameters
        use SBCalc_i                                              
        use SBCalc_z
        use SBSpline_functions
        use SBPlot
        call read_data                                              ! read data from the command line and from input files
        call spline_start_1                                         ! calculate spline coefficients for ImKs,DOS_Ms,DOS_SCs,pos
        call calc_zero_EF                                           ! calculate EFermi corresponding to undopped GaAs (stored in Efermi_00)
        call calc_po00                                              ! check po_h and po_e on infinity
        call calc_DOS_int                                           ! calculate integral of DOS to check number of valence electrons
        call calc_CBS_limits                                        ! calculate energy limits for each CBS band (using spline functions)
        call calc_scf                                               ! find SCF solution for the system M-SC
        call write_results
        call print_results
       end program SB














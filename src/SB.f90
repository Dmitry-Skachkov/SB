
!
!      Program Schottky Barrier to calculate Schottky contact parameters from first-principles
!
! (C)  Copyright 2021, Skachkov, Zhang, Cheng, Center for Molecular Magnetic Quantum Materials (M2QM),
!                      Quantum Theory Project (QTP),
!                      Department of Physics, University of Florida, Gainesville, FL 32611, USA 
!                      https://efrc.ufl.edu/
!                      Description of the method:      https://arxiv.org/abs/2001.00710
!                      Documentation and license:      https://github.com/m2qm-efrc/SB
!                      Requests for contributing:      Dr. Hai-Ping Cheng m2qm.efrc@phys.ufl.edu
!                      Questions regarding the code:   Dr. Dmitry Skachkov d.g.skachkov@gmail.com
!

       Program SB                                                   ! calculate Schottky barrier parameters
        use SBParameters                                            ! module with parameters of the system
        use SBCalc_i                                                ! module for integration
        use SBCalc_z                                                ! module for calculation on z-mesh and scf cycle
        use SBSpline_functions                                      ! module with all spline functions
        use SBPlot                                                  ! module for printing and writing results into the output files
        call print_logo                                             ! print common information about the method
        call read_data                                              ! read data from the command line and from input files
        call spline_start_1                                         ! calculate spline coefficients for ImKs,DOS_Ms,DOS_SCs,pos
        call calc_zero_EF                                           ! calculate EFermi corresponding to undopped GaAs (stored in Efermi_00)
        call calc_po00                                              ! check po_h and po_e on infinity
        call calc_DOS_int                                           ! calculate integral of DOS to check number of valence electrons
        call calc_CBS_limits                                        ! calculate energy limits for each CBS band (using spline functions)
        call calc_scf                                               ! find SCF solution for the system M-SC
        call calc_check_scf                                         ! check accuracy of scf cycle
        call write_results                                          ! write density po and potential V into gthe files
        call print_results                                          ! print all SB parameters
       end program SB














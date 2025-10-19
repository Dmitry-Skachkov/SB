
!
!      Program Schottky Barrier to calculate Schottky contact parameters from first-principles
!
! Copyright (C) 2021   Center for Molecular Magnetic Quantum Materials (M2QM),
!                      Quantum Theory Project (QTP),
!                      Department of Physics, University of Florida, Gainesville, FL 32611, USA 
!                      https://efrc.ufl.edu/
!                      Description of the method:      https://arxiv.org/abs/2001.00710
!                      Documentation and license:      https://github.com/m2qm-efrc/SB
!                      Requests for contributing:      Dr. Hai-Ping Cheng <m2qm.efrc@phys.ufl.edu>
!                      Questions regarding the code:   Dr. Dmitry Skachkov <dmitry.skachkov@dsedu.org>
!

       Program SB                                                   ! calculate Schottky barrier parameters
        use SBmpi
        use SBParameters                                            ! module with parameters of the system
        use SBCalc_i                                                ! module for integration
        use SBCalc_z                                                ! module for calculation on z-mesh and scf cycle
        use SBSpline_functions                                      ! module with all spline functions
        use SBPlot                                                  ! module for printing and writing results into the output files
        call P_start                                                !
        call P_calc_group(Nz)                                       !
        call print_logo                                             ! print common information about the method
        call read_data                                              ! read data from the command line and from input files
! call P_wait
!        write(iPrint,*) 'Spline_start'
        call spline_start                                           ! calculate spline coefficients for ImKs,DOS_Ms,DOS_SCs,pos
        if(CNL==0.d0) call calc_CNL                                 ! calculate CNL for interface
! call P_wait
!        write(iPrint,*) 'calc_zero_EF'
        call calc_zero_EF                                           ! calculate EFermi corresponding to undopped GaAs
! call P_wait
!        write(iPrint,*) 'calc_po00'
        call calc_po00                                              ! check po_h and po_e on infinity
        if(L_debug) call calc_DOS_int                               ! calculate integral of DOS to check number of valence electrons
!        write(iPrint,*) ' calc_scf'
! call P_wait
        call calc_scf                                               ! find SCF solution for the system M-SC
!        write(iPrint,*) ' calc_check_scf'
! call P_wait
        call calc_check_scf                                         ! check accuracy of scf cycle
!        write(iPrint,*) ' calc_LW'
! call P_wait
        call calc_LW                                                ! calculate DLW and ILW
!        write(iPrint,*) ' write_results'
! call P_wait
        call write_results                                          ! write density po and potential V into the files
!        write(iPrint,*) ' print_results'
! call P_wait
        call print_results                                          ! print all SB parameters
!        write(iPrint,*) ' P_stop'
! call P_wait
        call P_stop
       end program SB
















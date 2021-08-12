

!      Program Schottky Barrier to calculate Schottky contact parameters from first-principles
!
! Copyright (C) 2021   Center for Molecular Magnetic Quantum Materials (M2QM),
!                      Quantum Theory Project (QTP),
!                      Department of Physics, University of Florida, Gainesville, FL 32611, USA 
!                      https://efrc.ufl.edu/

     Module SBParameters
      implicit none
      real(8), parameter     :: pi = 3.141592653589793238462643383279502884197169d0
      real(8), parameter     :: ee = 2.7182818284590452353602874713527d0
      integer, parameter     :: Nitscf  = 1                    ! number of cycles for scf deltaE 
      integer, parameter     :: Nitscf0 = 4                    ! number of attempts to find approximate solution
      integer, parameter     :: Nitscf1 = 3                    ! number of attempts to find accurate solution
      integer, parameter     :: Nitscf2 = 1                    ! number of cycles for post scf
      integer, parameter     :: Nitscf3 = 3                    ! number of cycles for experimental scf
      integer, parameter     :: Nitmax  = 600                  ! max number of iterations
      logical, parameter     :: L_debug = .false.              ! for detailed printing
      logical, parameter     :: L_super_debug = .false.        ! for super detailed printing
      real(8), parameter     :: Lsc_inf = 1.d8                 ! numerical length of semi-infinite semiconductor (in A)
      real*8,  parameter     :: kb=8.6173430060d-5             ! Boltsman const in eV/K
      real(8), parameter     :: e0_SI = 8.8541878128d-12       ! in C/(V*m)
      real(8), parameter     :: C1 = 6.24150975d18             ! 1 Coulomb in e
      real(8), parameter     :: e0 = e0_SI*C1*1.d-10           ! e0 in e/(V*Angstrom)
      real(8), parameter     :: BohrA = 0.52917721067121d0     ! Bohr to A
      real(8), parameter     :: Eau = 36.3609d0                ! 1 a.u. for electric field in V/Ang
      integer, parameter     :: Nz = 704                       ! number of points in z (10-4, 10-3, 10-2, 0.1....)
      integer, parameter     :: Nptm = 10000                   ! max number of energy points
      real(8), allocatable   :: PDOS3(:,:)                     ! for separation D_L and D_H
      real(8), allocatable   :: PDOS4(:,:)                     !
      real(8), allocatable   :: DOS_M(:,:)                     ! DOS for interface (E,kx,ky)
      real(8), allocatable   :: dpol(:)                        ! polarization (a.u.)
      real(8), allocatable   :: Epol(:)                        ! electric field points for polarization (a.u.)
      real(8), allocatable   :: kp(:,:)                        ! k-points for integration
      real(8), allocatable   :: wk(:)                          ! weights for k-space integration 
      real(8), allocatable   :: ImKL1(:,:)                     ! light holes
      real(8), allocatable   :: ImKH1(:,:)                     ! heavy holes
      real(8), allocatable   :: Ef1(:)                         ! energy points for CBS
      integer                :: Nk                             ! number of k-points                  
      integer                :: Npt1                           ! number of points for CBS
      integer                :: N_DOS_M                        ! number of points for DOS_M interface
      integer                :: N_DOS_SC                       ! number of points for DOS_SC semiconductor
      integer                :: Npol                           ! number points for polarization
      integer                :: iter                           ! iteration number
      integer                :: k9                             ! nearest k-point to Gamma
      real(8)                :: DOS_Mtot(Nptm)                 ! DOS for interface (E) integrated over kx,ky
      real(8)                :: DOS_Sc(Nptm)                   ! DOS for bulk semiconductor
      real(8)                :: Ef_DOS_M(Nptm)                 ! E points for DOS_M
      real(8)                :: Ef_DOS_SC(Nptm)                ! E points for DOS_SC
      real(8)                :: Efi3(Nptm)                     ! for separation D_L and D_H
      real(8)                :: Efi4(Nptm)                     !
      real(8)                :: DOS0(Nptm),Efi0(Nptm)          ! PDOS of the surface
      real(8)                :: Sig                            ! charge surface density on the interface (in e/Angstrom**2), see Eqn.(7)
      real(8)                :: SigS                           ! initial surface charge on the interface
      real(8)                :: SigMIGS                        ! equivalent charge surface density of the MIGS (in e/Angstrom**2)
      real(8)                :: Sig_e                          ! equivalent charge surface density of the electrons (in e/Angstrom**2)
      real(8)                :: Sig_h                          ! equivalent charge surface density of the holes (in e/Angstrom**2)
      real(8)                :: zp                             ! z as a parameter for DMIGS (in Angstrom)
      real(8)                :: ECBM                           ! CBM (in eV)
      real(8)                :: EVBM                           ! VBM (in eV)
      real(8)                :: Temp                           ! temperature in K
      real(8)                :: kbT                            ! = kb*Temp (in eV)
      real(8)                :: er                             ! dielectric constant of SC from DFT
      real(8)                :: kappa                          ! kappa for P = kappa*E from DFT calculation  er = 1 + kappa/e0
      real(8)                :: cz                             ! size of the cell in z direction for CBS (for QE units)
      real(8)                :: ckA                            ! coefficient for ImK to convert to 1/A
      real(8)                :: EFermi1                        ! E Fermi in SC
      real(8)                :: EFermi2                        ! E Fermi in M
      real(8)                :: Emin1,Emax1                    ! CBS limits
      real(8)                :: V_el1(Nz)                      ! electrostatic potential 
      real(8)                :: V_el0(Nz)                      ! electrostatic potential (from previous step) 
      real(8)                :: V_eln(Nz)                      ! electrostatic potential on next iteration (mixed with alfa) 
      real(8)                :: V0                             ! parameters of initial distribution V(z)   
      real(8)                :: alfa                           ! mixing parameter for po
      real(8)                :: alfa_V                         ! mixing parameter for V
      real(8)                :: alfa_Sig                       ! mixing parameter for Sigma
      real(8)                :: alfa_MIGS                      ! mixing parameter for MIGS
      real(8)                :: alfa_dipole                    ! mixing parameter for dipole correction
      real(8)                :: po_0(Nz)                       ! charge density initial
      real(8)                :: po1(Nz)                        ! charge density calculated by formulas from DOS
      real(8)                :: po_new(Nz)                     ! charge density after mixing
      real(8)                :: po_MIGS(Nz)                    ! po_MIGS to calculate Sigma_MIGS
      real(8)                :: po_h(Nz)                       ! po_h holes 
      real(8)                :: po_e(Nz)                       ! po_e electrons 
      real(8)                :: Zz(Nz)                         ! z mesh
      real(8)                :: po0                            ! parameters of initial distribution, see Eqn.(20)     
      real(8)                :: V_D0                           ! volume (in A^3) of the interface
      real(8)                :: V_DSC                          ! volume (in A^3) of the semiconductor 
      real(8)                :: po00                           ! - doping concentration 
      real(8)                :: po00_h                         ! density of holes at +infinity
      real(8)                :: po00_e                         ! density of electrons at +infinity
      real(8)                :: po00_h0                        ! density of holes at +infinity
      real(8)                :: po00_e0                        ! density of electrons at +infinity
      real(8)                :: po00_hS                        ! density of holes at interface
      real(8)                :: po00_eS                        ! density of electrons at interface
      real(8)                :: Nee1                           ! norm for DOS_SC
      real(8)                :: Nee2                           ! norm for DOS_M
      real(8)                :: poh_max                        ! max of poh(z)
      real(8)                :: poe_max                        ! max of poe(z)
      real(8)                :: poMe_max                       ! max of poSIGS_e(z)
      real(8)                :: delta_po                       ! delta (po_1 - po0)
      real(8)                :: delta_V                        ! delta (V_el1 - V_el0)
      real(8)                :: El_f(Nz)                       ! gradients of potential (electric field)
      real(8)                :: El_f2(Nz)                      ! gradients of potential (electric field)   for TEST (from dipole contribution)
      real(8)                :: Volpr                          ! volume of the primitive cell for GaAs
      real(8)                :: DS0                            ! parameter for Sigma on the surface
      real(8)                :: CNL                            ! charge neutrality level
      real(8)                :: SigmaS                         ! initial charge on the interface
      real(8)                :: SBH                            ! Schottky barrier height (SBH)
      real(8)                :: sumk                           ! sum of wk(k)
      real(8)                :: a2p                            ! 2pi/a for CBS
      real(8)                :: alat                           ! lattice parameter (A)
      real(8)                :: b1(3),b2(3),b3(3)              ! reciprocal vectors
      real(8)                :: zconnect                       ! connection numerical and analytical solution
      real(8)                :: z3,z4                          ! 3d and 4th layers
      real(8)                :: za                             ! parameter z0 in Eqn.(20)
      real(8)                :: EFermi_00                      ! charge neutrality level in the bulk semiconductor
      real(8)                :: Efermi_input                   ! Fermi level in input to set up in the system
      real(8)                :: Exxh1(4),Exxh2(4)              ! limits for DOS of bulk for holes (4 intervals for intagration)
      real(8)                :: Exxe1(4),Exxe2(4)              ! limits for DOS of bulk for electrons
      real(8)                :: Exxm1(4),Exxm2(4)              ! limits for PDOS for MIGS
      real(8)                :: dEf                            ! filling level of the surface delta_E, see Eqn.(18) 
      real(8)                :: Surface                        ! surface size of the cell
      real(8)                :: Lz                             ! length of the cell in the interface calculation (in A) 
      real(8)                :: Lz_int                         ! width of the interfacial layer (in A) 
      real(8)                :: Sig_gate                       ! surface charge on the gate (in cm-2)  
      real(8)                :: Lsc                            ! the length of the SC (A)
      real(8)                :: E_00                           ! electric field close to SC surface
      real(8)                :: P_00                           ! polarization close to SC surface
      real(8)                :: gap                            ! band gap of semiconductor (eV)
      real(8)                :: DLW                            ! depletion layer width (DLW)
      real(8)                :: ILW                            ! inversion layer width (ILW)
      character(10)          :: LscA                           ! the length of the SC
      logical                :: L_restart                      ! restart from previous solution
      logical                :: L_inf                          ! infinite or finite SC
      logical                :: L_p_type                       ! p-type semiconductor
      logical                :: L_n_type                       ! n-type semiconductor
      logical                :: L_conv = .false.               ! reach convergency
      logical                :: L_scf  = .false.               ! divergency
      logical                :: L_exp_scf  = .false.           ! experimental scf
     contains







       subroutine read_data
        call getargR(1,Temp)                                        ! temperature                                          
        call getargR(2,EFermi_input)                                ! Fermi level                               
        call getarg(3,LscA)                                         ! length of SC ('inf' or the value in A)
        if(L_debug) print *,' LscA=',LscA
        if(trim(adjustl(LscA)) == 'inf') then
         L_inf = .true.
         Lsc = Lsc_inf
        else
         L_inf = .false.
         read(LscA,*) Lsc                                           ! convert character(10) LscA to real(8) Lsc
        endif  
        if(L_debug) print *,' Lsc=',Lsc
        call getargR(4,Sig_gate)                                    ! gating charge density (in cm-2)
        if(L_inf) Sig_gate = 0.d0
        if(L_debug) print *,' Sig_gate=',Sig_gate
        Sig_gate = Sig_gate*1.D-16                                  ! convert to A^-2
        if(Temp == 0.d0) then
         print *,'Please use nonzero temperature'
         stop
        endif
        kbT = kb*Temp 
        if(L_check_file('restart.dat')) then
         L_restart = .true.                                         ! continue from previous calculation
         if(L_super_debug) print *,'continue calculations from previous step'
         call read_restart_dat
        else
         L_restart = .false.                                        ! start from scratch
         if(L_super_debug) print *,'start calculation from scratch'
        endif
        call read_input_dat  
        SigS    =  SigS**1.D-16                                     ! convert to A^-2
        EVBM    =  0.00d0                                           ! VBM
        ECBM    =  gap                                              ! CBM
        alat    =  alat*BohrA                                       ! convert to A
        V_D0    =  V_D0*BohrA**3                                    ! convert to A^3
        Surface =  V_D0/Lz                                          ! surface of the cell
        V_D0    =  V_D0/Lz*Lz_int                                   ! volume of the interfacial layer for PDOS (for MIGS)
        V_DSC   =  V_DSC*BohrA**3                                   ! volume of bulk semiconductor cell (in A^3) 
        ckA     =  2.d0*pi/cz                                       ! coefficient for ImK to convert to 1/A
        call calc_reciprocal_param
        call read_k_mesh    
        call read_pol                                               ! read polarization data                                    
        call print_input_parameters
        call read_CBS_data    
        call calc_CBS_limits                                        ! calculate energy limits for each CBS band (using spline functions)
        call set_limits_DOS                                         ! set limits for integration
        call read_PDOS                                              ! read DOS of SC and PDOS of interfacial layer   
        EFermi2 = EFermi1 
        alfa   = 1.0d0                                              ! mixing po
        alfa_MIGS = 1.d0
        alfa_dipole = 1.d0
   !     alfa_V   = 0.0d0                                           ! mixing V
   !     alfa_Sig =-0.01d0                                          ! mixing Sigma
        alfa_V = 1.d0
        alfa_Sig = 1.d0  
       end subroutine read_data




       subroutine read_input_dat
        if(L_super_debug) print *,'read input.dat'
        open(unit=1,file='input.dat')    
         read(1,*) V_D0                                            ! V_D0 in a.u. (from QE output)
         read(1,*) Lz                                              ! Lz (A)
         read(1,*) Lz_int                                          ! Lz_int (A)
         read(1,*) CNL                                             ! CNL (eV)
         read(1,*) z3                                              ! z3 (A)
         read(1,*) z4                                              ! z4 (A)
         read(1,*) alat                                            ! alat (au from QE output) 
         read(1,*) b1(1:3)                                         ! reciprocal vector b1 (in crystal representation from QE output)
         read(1,*) b2(1:3)                                         !
         read(1,*) b3(1:3)                                         !
         read(1,*) cz                                              ! cz (A)
         read(1,*) V_DSC                                           ! V_DSC volume of semiconductor in a.u. (from QE output)
         read(1,*) gap                                             ! band gap of the bulk (eV)
         read(1,*) SigS                                            ! initial surface charge on interface (in cm-2)
        close(unit=1)
       end subroutine read_input_dat




       subroutine read_restart_dat
        real(8)            :: Temp0,EFermi_input0,Lsc0,Sig_gate0
        if(L_super_debug) print *,'read restart.dat'
        open(unit=1,file='restart.dat')    
         read(1,*) za
         read(1,*) Sig
         read(1,*) Temp0
         read(1,*) EFermi_input0
         read(1,*) Lsc0
         read(1,*) Sig_gate0
        close(unit=1)
        if(dabs((Temp0-Temp)/Temp) > 0.2d0) L_restart = .false.                  ! if Temp change more then 20% start from scratch
        if(Sig_gate /= 0.d0 .and. Sig_gate0 /= 0.d0) then                       
         if(dabs((Sig_gate0-Sig_gate)/Sig_gate) > 0.5d0) L_restart = .false.     ! if Sig_gate change more then 50% start from scratch
        else
         L_restart = .true.
        endif
        if(dabs(EFermi_input0-EFermi_input) > 0.25d0) L_restart = .false.        ! if EFermi change more then 0.25 eV start from scratch
        if(trim(adjustl(LscA)) /= 'inf') then
         if(Lsc <= 1.d5) then
          if(dabs((Lsc0-Lsc)/Lsc) > 10.d0) L_restart = .false.                   ! if Lsc change more then x 10 start from scratch
         else
          L_restart = .true.
         endif
        endif
        if(.not.L_restart) print 1
 1      format(/' Parameters of the systems are changed significantly from previous calculation'/' starting from scratch'/)
       end subroutine read_restart_dat



       subroutine write_restart_dat
        if(L_super_debug) print *,'read restart.dat'
        open(unit=1,file='restart.dat')    
         write(1,*) za
         write(1,*) Sig
         write(1,*) Temp
         write(1,*) EFermi_input
         write(1,*) Lsc
         write(1,*) Sig_gate
        close(unit=1)
       end subroutine write_restart_dat





       subroutine print_input_parameters
         if(trim(adjustl(LscA)) == 'inf') then
          print 20
          if(L_debug) print 22,Lsc_inf 
         else
          print 21,Lsc
         endif
         print 2,Lz
         print 3,Lz_int
         print 5,V_D0
         print 19,Surface
         print 7,z3
         print 8,z4
         print 9,alat
         print 11,cz
         print 12,V_DSC
         print 15,EVBM
         print 16,ECBM
         print 17,EFermi_input
         print 18,Temp
         print 23,er
         print 10,b1,b2,b3
         if(L_debug) then
          print 13,e0
          print 14,ckA
         endif
 2       format(' length of the cell in the interface calculation (Lz) = ',   F12.4,' A')
 3       format(' width of the interfacial layer              (Lz_int) = ',   F12.4,' A') 
 5       format(' volume of the interfacial layer               (V_D0) = ',   F12.4,' A^3')
 7       format(' position of 3d layer of the interface           (z3) = ',   F12.4,' A')
 8       format(' position of 4th layer of the interface          (z4) = ',   F12.4,' A')
 9       format(' crystal parameter                             (alat) = ',   F12.4,' A')
10       format(' reciprocal vectors (in 2*pi/alat)'/3(3F15.3/))
11       format(' size of the cell in z direction for CBS         (cz) = ',   F12.4,' A')
12       format(' volume of the cell of the bulk semiconductor (V_DSC) = ',   F12.4,' A^3')
13       format(' constant e0                                          = ',1p,E12.4,' e/(V*A)')
14       format(' ckA                                                  = ',1p,E12.4,' 1/A')
15       format(' Valence band maximum                          (EVBM) = ',   F12.4,' eV')
16       format(' Conduction band minimum                       (ECBM) = ',   F12.4,' eV')
17       format(' Fermi level for the system                  (EFermi) = ',   F12.4,' eV')
18       format(' Temperature                                   (Temp) = ',   F12.4,' K')
19       format(' Surface of interface cell                  (Surface) = ',   F12.4,' A^2')
20       format(' length of the semiconductor                    (Lsc) =      infinity')
21       format(' length of the semiconductor                    (Lsc) = ',1p,E12.1,' A')
22       format(' numerical length of "infinite" semiconductor   (Lsc) = ',1p,E12.1,' A')
23       format(' dielectric constant of the bulk                 (er) = ',   F12.4)
       end subroutine print_input_parameters





     subroutine read_CBS_data    
      integer                :: rNk,k,k1,i
      if(L_super_debug) then
       print *
       print *,'READ CBS'
       print 4
      endif
      open(unit=1,file='cbs.dat') 
       read(1,1) rNk,Npt1
       if(L_super_debug) then 
        print *,'Npt1=',Npt1
        print *,'allocate arrays'
       endif
       allocate(ImKL1(Npt1,Nk))
       allocate(ImKH1(Npt1,Nk))
       allocate(Ef1(Npt1))
       if(rNk/=Nk) then
        print *,'ERROR with Nk'
        stop
       endif
       do k=1,Nk
        read(1,2) k1
        if(k1/=k) then
         print *,'ERROR with k1'
         stop
        endif
        do i=1,Npt1
         read(1,3) Ef1(i),ImKL1(i,k),ImKH1(i,k)
         if(L_super_debug) print 3,Ef1(i),ImKL1(i,k),ImKH1(i,k)
        enddo
       enddo
      close(unit=1)
 1    format(2I5)
 2    format(I5)
 3    format(3F14.6)
 4    format('   Ef1      ImKL1      ImKH1')
     end subroutine read_CBS_data







       subroutine set_limits_DOS                 ! set limits for DOS integration
        Exxh1(1) =  EVBM - 6.d0                  ! limits for integration for h
        Exxh2(1) =  EVBM - 4.d0 
        Exxh1(2) =  EVBM - 4.d0
        Exxh2(2) =  EVBM - 2.d0 
        Exxh1(3) =  EVBM - 2.d0
        Exxh2(3) =  EVBM - 1.d0 
        Exxh1(4) =  EVBM - 1.d0 
        Exxh2(4) =  EVBM
        Exxe1(1) =  ECBM                         ! limits for integration for e
        Exxe2(1) =  ECBM + 1.d0  
        Exxe1(2) =  ECBM + 1.d0  
        Exxe2(2) =  ECBM + 2.d0  
        Exxe1(3) =  ECBM + 2.d0  
        Exxe2(3) =  ECBM + 4.d0  
        Exxe1(4) =  ECBM + 4.d0  
        Exxe2(4) =  ECBM + 6.d0  
        Exxm1(1) =  EVBM                         ! set limits for MIGS integration (4 intervals for accurate integration)
        Exxm2(1) =  0.10d0*(EVBM+ECBM)
        Exxm1(2) =  0.10d0*(EVBM+ECBM)  
        Exxm2(2) =  0.25d0*(EVBM+ECBM)
        Exxm1(3) =  0.25d0*(EVBM+ECBM)
        Exxm2(3) =  0.99d0*(EVBM+ECBM)
        Exxm1(4) =  0.99d0*(EVBM+ECBM)
        Exxm2(4) =  ECBM
       end subroutine set_limits_DOS







     subroutine calc_reciprocal_param
      a2p = 2.d0*pi/alat 
      b1 = b1*a2p
      b2 = b2*a2p
      b3 = b3*a2p
     end subroutine calc_reciprocal_param







       subroutine plot_DOS_Mtot_test
        integer     :: j
        open(unit=3,file='DOS_Mtot.dat')
        write(3,3) 'DOS_Mt',N_DOS_M
        do j=1,N_DOS_M
          write(3,2) Ef_DOS_M(j),DOS_Mtot(j)
        enddo
        close(unit=3)
 2      format(5x,F8.3,E12.4)
 3      format(/'VARIABLES = "E", "PDOS"'/'ZONE T="',A7,'" I=',I4,' F=POINT')
       end subroutine plot_DOS_Mtot_test









       subroutine calc_DOS_Mtot(PDOS3)
        real(8)     :: DOSx
        real(8)     :: PDOS3(Nptm,Nk)
        integer     :: j
        integer     :: k
        if(L_super_debug) print *,'calc_DOS_Mtot:'
        do j=1,N_DOS_M
         DOSx = 0.d0
         do k=1,Nk
          DOSx = DOSx + PDOS3(j,k)*wk(k)
         enddo
         DOS_Mtot(j) = DOSx/sumk
        enddo
       end subroutine calc_DOS_Mtot







       subroutine read_k_mesh   
        integer      :: k
        real(8)      :: kr
        real(8)      :: kr9
        if(L_debug) print *,'open file k_mesh.dat'
        open(unit=2,file='k_mesh.dat')
         read(2,*) Nk
         if(L_debug) print *,'Nk=',Nk
         allocate(kp(3,Nk))
         allocate(wk(Nk))
         kp(1:3,1:Nk) = 0.d0
         do k=1,Nk
          read(2,*) kp(1:2,k),wk(k)
         enddo
        close(unit=2)
        if(L_debug) print *,'read k-mesh with weights'
        sumk = 0.d0
        do k=1,Nk
         sumk = sumk + wk(k)
        enddo
        if(L_debug) print *,'TEST wk:', sumk
        if(L_debug) print *,'k kx ky kr'
        kr9 = 100.d0
        do k=1,Nk
         kr = dsqrt(kp(1,k)**2+kp(2,k)**2)
         if(k > 1 .and. kr < kr9) then
          k9 = k 
          kr9 = kr
         endif
         if(L_debug) print 4,k,kp(1,k),kp(2,k),kr
        enddo
 1      format(20x,3F12.7,7x,F12.7)
 2      format(3F16.10)
 3      format(I5)
 4      format(I4,3F15.5)
       end subroutine read_k_mesh








   subroutine getargR(i,R)                ! read real from argument line
    integer       :: i
    real(8)       :: R
    character(10)  :: Rc
    call getarg(i,Rc)
    read(Rc,*) R
   end subroutine getargR




   subroutine getargI(i,N)                ! read integer from argument line
    integer       :: i
    integer       :: N
    character(7)  :: Nc
    call getarg(i,Nc)
    read(Nc,*) N
   end subroutine getargI






     subroutine read_PDOS   
      integer       :: j                       ! energy 
      integer       :: k                       ! k-points
      allocate(PDOS3(Nptm,Nk))
      allocate(PDOS4(Nptm,Nk))
      allocate(DOS_M(Nptm,Nk))
      call read_pdos_1(Efi3,PDOS3,3)           ! PDOS of 3d layer of interface
      call read_pdos_1(Efi4,PDOS4,4)           ! PDOS of 4th layer of interface
      call read_pdos_0                         ! PDOS of the surface
      if(L_debug) print *,'N_DOS_M=',N_DOS_M
      do k=1,Nk
       do j=1,N_DOS_M
        Ef_DOS_M(j) = Efi3(j)                  ! copy PDOS3
        DOS_M(j,k) = PDOS3(j,k)
       enddo
      enddo
      call calc_DOS_Mtot(PDOS3)                ! calculate integrated DOS_M of interfacial layer
      call open_file(2,'dos_bulk_.dat')        ! DOS of bulk SC
      read(2,*) N_DOS_SC
      if(L_debug) print *,'N_DOS_SC=',N_DOS_SC
      do j=1,N_DOS_SC
       read(2,*) Ef_DOS_SC(j),DOS_SC(j)
       if(DOS_SC(j) < 0.d0) then
        DOS_SC(j) = 0.d0
       endif
      enddo
      close(unit=2)
      if(L_super_debug) print *,'read_PDOS: read ',N_DOS_SC,' points of bulk'
     end subroutine read_PDOS






     subroutine read_pdos_1(Efi,PDOS,ilayer)                                  ! read PDOS for layer 3(4)
      real(8)           :: PDOS(Nptm,Nk),Efi(Nptm)
      integer           :: ilayer
      integer           :: k,j,kr
      if(L_debug) print *,'read_PDOS:'
      call open_file(2,'kpdos_int_'//trim(adjustl(fstr(ilayer)))//'.dat')     ! PDOS of interface (E,k)
      read(2,*) kr,N_DOS_M 
      if(L_debug) print *,'read N_DOS_M =',N_DOS_M
      do k=1,Nk                                                               ! read Nk points
       read(2,*)
       do j=1,N_DOS_M
        read(2,*) Efi(j),PDOS(j,k)
       enddo
      enddo
      if(L_debug) print *,'read_PDOS: read ',N_DOS_M,' points of interface and ',Nk,' k-points'
      close(unit=2)
     end subroutine read_pdos_1




     subroutine read_pdos_0                    ! read PDOS for surface (1st layer + metal surface)
      integer           :: j
      if(L_debug) print *,'read_pdos_0:'
      call open_file(2,'DOStot.dat')    
       read(2,*) 
       do j=1,N_DOS_M
        read(2,*) Efi0(j),DOS0(j)
       enddo
       if(L_super_debug) then
        print *,'read_pdos_0: read ',N_DOS_M,' points of surface layer'
        do j=1,N_DOS_M
         print 11,Efi0(j),DOS0(j)
        enddo
       endif
      close(unit=2)
 1    format(15x,I4)
11    format(F14.5,E17.5)
     end subroutine read_pdos_0






     character(len=7) function fstr(k)                             !   Convert an integer to character*7
      integer, intent(in) :: k
      write (fstr,'(I7)') k
     end function fstr






      subroutine open_file(un,name)
       integer      :: un
       character(*) :: name
       if(L_check_file(trim(adjustl(name)))) then
        if(L_debug) print *,' Open file ',name
        open(unit=un,file=trim(adjustl(name)))
!        L_file = .true.
       else
        print *,'ERROR:  File ',trim(adjustl(name)),' does not exist'
        stop
       endif
      end subroutine open_file





    logical function L_check_file(name)
     character(*)   name
     inquire(file=name,EXIST=L_check_file)
    end function L_check_file








    subroutine read_pol                                             ! read polarization data    
     integer                :: i
     if(L_super_debug) then
      print *
      print *,' open polarization.dat'
      print *,'    E     polarization'
     endif
     open(unit=11,file='polarization.dat')
      read(11,*) Npol
      if(Npol == 1) then
       read(11,*) er
       kappa = e0*(er-1.d0)
      elseif(Npol > 1) then
       if(L_debug) print *,'Npol=',Npol
       allocate(Epol(Npol))
       allocate(dpol(Npol))
       do i=1,Npol
        read(11,*) Epol(i),dpol(i)
        if(L_super_debug) print 1,Epol(i),dpol(i)
       enddo
       kappa = 0.d0
       do i=2,Npol
        kappa = kappa + (dpol(i)/Epol(i))
       enddo
       kappa = 1.d0/dfloat(Npol)*kappa
       if(L_debug) print *,'read_pol:   kappa=',kappa
       if(L_debug) print *,'read_pol:   kappa/e0=',kappa/e0
       er = kappa/e0 + 1.d0
       deallocate(Epol)
       deallocate(dpol)
      else
       print *,'read_pol: *ERROR* Npol = 0'
       print *,'read_pol: no data in polarization.dat file'
       stop
      endif
     close(unit=11)
 1   format(F11.6,F12.7)
    end subroutine read_pol 





      subroutine calc_CBS_limits
       Emin1 = Ef1(1)
       Emax1 = Ef1(Npt1)
       if(L_debug) then
       print *
       print *,'calc_CBS_limits'
       print 2
       print 1,Emin1,Emax1
       endif
 1     format(2E12.4)
 2     format('     Emin1           Emax1 ')
      end subroutine calc_CBS_limits
 





  end module SBParameters












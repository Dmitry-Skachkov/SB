

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
      integer, parameter     :: Nptm = 10000                   ! max number of energy points
      integer                :: Nk                             ! number of k-points                    ! 27 for GaAs    121 for Si-Al
      real(8), allocatable   :: DOS_M(:,:)                     ! DOS for interface (E,kx,ky)
      real(8)                :: DOS_Mtot(Nptm)                 ! DOS for interface (E) integrated over kx,ky
      real(8)                :: DOS_Sc(Nptm)                   ! DOS for bulk semiconductor
      real(8)                :: Ef_DOS_M(Nptm)                 ! E points for DOS_M
      real(8)                :: Ef_DOS_SC(Nptm)                ! E points for DOS_SC
      real(8), allocatable   :: PDOS3(:,:)                     ! for separation D_L and D_H
      real(8)                :: Efi3(Nptm)                     ! for separation D_L and D_H
      real(8), allocatable   :: PDOS4(:,:)
      real(8)                :: Efi4(Nptm)
      real(8)                :: DOS0(Nptm),Efi0(Nptm)          ! PDOS of the surface
      integer                :: N_DOS_M                        ! number of points for DOS_M interface
      integer                :: N_DOS_SC                       ! number of points for DOS_SC semiconductor
      real(8)                :: Sig                            ! charge surface density on the interface (in e/Angstrom**2)
      real(8)                :: SigMIGS                        ! equivalent charge surface density of the MIGS (in e/Angstrom**2)
      real(8)                :: Sig_e                          ! equivalent charge surface density of the electrons (in e/Angstrom**2)
      real(8)                :: Sig_h                          ! equivalent charge surface density of the holes (in e/Angstrom**2)
      real(8)                :: zp                             ! z as a parameter for DMIGS (in Angstrom)
      real(8)                :: ECBM                           ! CBM (in eV)
      real(8)                :: EVBM                           ! VBM (in eV)
      real*8,  parameter     :: kb=8.6173430060d-5             ! Boltsman const in eV/K
      real(8)                :: Temp                           ! temperature in K
      real(8)                :: kbT                            != kb*Temp                  ! in eV
      real(8), parameter     :: e0_SI = 8.8541878128d-12       ! in C/(V*m)
      real(8), parameter     :: C1 = 6.24150975d18             ! 1 Coulomb in e
      real(8), parameter     :: e0 = e0_SI*C1*1.d-10           ! e0 in e/(V*Angstrom)
      real(8)                :: er                             ! dielectric constant of GaAs (static) from DFT
      real(8)                :: kappa                          ! kappa for P = kappa*E from DFT calculation  er = 1 + kappa/e0
      real(8), parameter     :: BohrA = 0.52917721067121d0     ! Bohr to A
      real(8)                :: cz                             != 9.7917162d0               ! A 15.d0*BohrA*1.896d0       ! size of the cell in z direction for CBS (for QE units)
      real(8)                :: ckA                            ! = 2.d0*pi/cz               ! coefficient for ImK to convert to 1/A
      real(8)                :: EFermi1                        ! E Fermi in GaAs
      real(8)                :: EFermi2                        ! E Fermi in Graphene
      logical                :: L_conv = .false.               ! reach convergency
      logical                :: L_scf  = .false.               ! divergency
      real(8)                :: Emin1,Emax1                    ! CBS limits
      integer, parameter     :: Nz = 704                       ! number of points in z (10-6, 10-5, 10-4, 10-3, 10-2, 0.1....)
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
      real(8)                :: po_h(Nz)                       ! po_h holes (to write to the file only)
      real(8)                :: po_e(Nz)                       ! po_e electrons (to write to the file only)
      real(8)                :: Zz(Nz)                         ! z mesh
      integer                :: filen                          ! number of file for results
      real(8)                :: po0                            ! parameters of initial distribution po(z) = po0 * exp(- ImKz(E*,0)z)      !z/z0)
      real(8)                :: V_D0                           ! volume (in A^3) of the interface
      real(8)                :: V_DSC                          ! volume (in A^3) of the semiconductor 
      real(8)                :: po00                           ! density of natural doped electrons of SC (on infinity)
      real(8)                :: po00_h                         ! density of holes at +infinity
      real(8)                :: po00_e                         ! density of electrons at +infinity
      real(8)                :: po00_h0                        ! density of holes at +infinity
      real(8)                :: po00_e0                        ! density of electrons at +infinity
      real(8)                :: Nee1                           ! norm for DOS_SC
      real(8)                :: Nee2                           ! norm for DOS_M
      real(8)                :: poh_max                        ! max of poh(z)
      real(8)                :: poe_max                        ! max of poe(z)
      real(8)                :: poMe_max                       ! max of poSIGS_e(z)
      real(8)                :: delta_po                       ! delta (po_1 - po0)
      real(8)                :: delta_V                        ! delta (V_el1 - V_el0)
      integer                :: Nitscf                         ! number of iterations for scf
      integer                :: Nitscf0                        ! number of iterations for initial search
      integer                :: Nitscf2                        ! number of iterations second loop
      character(1)           :: Calc                           ! s - start; c - continue calculation
      integer                :: Npol                           ! number points for polarization
      real(8), allocatable   :: dpol(:)                        ! polarization (a.u.)
      real(8), allocatable   :: Epol(:)                        ! electric field points for polarization (a.u.)
      real(8)                :: El_f(Nz)                       ! gradients of potential (electric field)
      real(8)                :: El_f2(Nz)                      ! gradients of potential (electric field)   for TEST (from dipole contribution)
      real(8), parameter     :: Eau = 36.3609d0                ! 1 a.u. for electric field in V/Ang
      real(8)                :: Volpr                          ! volume of the primitive cell for GaAs
      real(8)                :: DS0                            ! parameter for Sigma on the surface
      real(8)                :: CNL                            ! charge neutrality level
      real(8)                :: SigmaS                         ! charge on the interface
      real(8)                :: SBH                            ! Schottky barrier height
      real(8), allocatable   :: kp(:,:)                        ! k-points for integration
      real(8), allocatable   :: wk(:)                          ! weights for k-space integration 
      real(8)                :: sumk                           ! sum of wk(k)
      real(8)                :: a2p                            ! 2pi/a for CBS
      real(8)                :: alat                           ! lattice parameter (A)
      real(8)                :: b1(3),b2(3),b3(3)              ! reciprocal vectors
      real(8)                :: zconnect                       ! connection numerical and analytical solution
      integer                 :: Npt1                          ! number of points for CBS
      real(8), allocatable   :: ImKL1(:,:)                     ! light holes
      real(8), allocatable   :: ImKH1(:,:)                     ! heavy holes
      real(8), allocatable   :: Ef1(:)                      ! energy points for CBS
      real(8)                :: z3,z4                          ! 3d and 4th layers
      real(8)                :: za                             ! initial approximation (width)
      real(8)                :: EFermi_00                      ! charge neutrality level in the bulk semiconductor
      real(8)                :: Efermi_input                   ! Fermi level in input to set up in the system
      real(8)                :: Exxh1(7),Exxh2(7)              ! limits for DOS of bulk for holes
      real(8)                :: Exxe1(7),Exxe2(7)              ! limits for DOS of bulk for electrons
      real(8)                :: Exxm1(7),Exxm2(7)              ! limits for PDOS for MIGS
      logical                :: L_p_type                       ! p-type semiconductor
      logical                :: L_n_type                       ! n-type semiconductor
      real(8)                :: dEf                            ! filling level of the surface 
      real(8)                :: Surface                        ! surface size of the cell
      real(8)                :: Lz                             ! length of the cell in the interface calculation (in A) 
      real(8)                :: Lz_int                         ! width of the interfacial layer (in A) 
      real(8)                :: Sig_gate                       ! surface charge on the gate (in cm-2)  
      real(8)                :: Lsc                            ! the length of the SC
      character(10)          :: LscA                           ! the length of the SC
      logical                :: L_inf                          ! infinite or finite SC
      real(8)                :: E_00                           ! electric field close to SC surface
      real(8)                :: P_00                           ! polarization close to SC surface
      logical                :: L_debug                        ! for detailed printing
      logical                :: L_super_debug                  ! for super detailed printing
     contains







       subroutine read_data
        call getargR(1,Temp)                    ! temperature                                          
        call getargR(2,EFermi_input)            ! Fermi level                               
        call getarg(3,LscA)                     ! length of SC ('inf' or the value in A)
        if(L_debug) print *,' LscA=',LscA
        if(trim(adjustl(LscA)) == 'inf') then
         L_inf = .true.
         Lsc = 1.d6
        else
         L_inf = .false.
         read(LscA,*) Lsc                       ! convert character(10) LscA to real(8) Lsc
        endif  
        if(L_debug) print *,' Lsc=',Lsc
        call getargR(4,Sig_gate)               ! gating charge density
        if(L_inf) Sig_gate = 0.d0
        if(L_debug) print *,' Sig_gate=',Sig_gate
        Sig_gate = Sig_gate*1.D-16             ! convert to A^-2
        kbT = kb*Temp
        Calc = 's'
        if(Calc=='s') print *,'start calculations'
        if(Calc=='c') print *,'continue calculations from previous step'
        call read_k_mesh      
        print *,'open input.dat'
        open(unit=1,file='input.dat')    
         print *,'open input.dat'
         read(1,*) V_D0
         print *,'V_D0=',V_D0
         read(1,*) Lz
         print *,'Lz=',Lz
         read(1,*) Lz_int
         print *,'Lz_int=',Lz_int
         read(1,*) CNL
         print *,'CNL=',CNL
         read(1,*) z3
         print *,'z3=',z3
         read(1,*) z4
         print *,'z4=',z4
         read(1,*) alat
         print *,'alat=',alat
         read(1,*) b1(1:3)
         print *,'b1=',b1
         read(1,*) b2(1:3)
         print *,'b2=',b2
         read(1,*) b3(1:3)
         print *,'b3=',b3
         read(1,*) cz
         print *,'cz=',cz
         read(1,*) V_DSC
         print *,'V_DSC=',V_DSC
        close(unit=2)
        call calc_reciprocal_param
        call read_CBS_data    
        call set_limits_GaAs                                       ! set limits for integration
        call read_PDOS                                             ! read DOS of SC and PDOS of interfacial layer   
        call read_pol                                              ! read polarization data                                    
        EFermi2 = EFermi1
        print *,'Temp=',Temp,' K'
        print *,'EFermi1=',EFermi1,' eV'
        print *,'EFermi2=',EFermi2,' eV'
        alfa   = 1.0d0                                             ! mixing po
        alfa_MIGS = 1.d0
        alfa_dipole = 1.d0
   !     alfa_V   = 0.0d0                                          ! mixing V
   !     alfa_Sig =-0.01d0                                         ! mixing Sigma
        alfa_V = 1.d0
        alfa_Sig = 1.d0  
        print *
        print *
        print *
        print *,'e0=',e0,' e/(V*A)'
        print *,'cz=',cz,' A'
        print *,'ckA=',ckA,' 1/A'
        print *,'EVBM=',EVBM,' eV'
        print *,'ECBM=',ECBM,' eV'
        print *,'cz=',cz,' A'
        print *,'V0=',V0,' eV'
        print *,'EFermi_input=',EFermi_input,' eV'
        Nitscf  = 1                                      ! cycle over deltaE
        Nitscf0 = 3                                      ! pre scf cycle
        Nitscf2 = 1                                      ! post scf cycle
        L_debug = .false.                                 ! print all intermediate results
        L_super_debug = .false.                          ! extreme printing
       end subroutine read_data






     subroutine read_CBS_data    
      integer                :: rNk,k,k1,i
      if(L_debug) then
       print *
       print *,'READ CBS'
       print 4
      endif
      open(unit=1,file='cbs.dat') 
       read(1,1) rNk,Npt1
       if(L_debug) then 
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
         if(L_debug) print 3,Ef1(i),ImKL1(i,k),ImKH1(i,k)
        enddo
       enddo
      close(unit=1)
      ckA = 2.d0*pi/cz                                    ! coefficient for ImK to convert to 1/A
 1    format(2I5)
 2    format(I5)
 3    format(3F14.6)
 4    format('   Ef1      ImKL1      ImKH1')
     end subroutine read_CBS_data







       subroutine set_limits_GaAs                               ! set limits for DOS
      Exxh1(1) = -6.900d0 
      Exxh2(1) = -6.561d0 
      Exxh1(2) = -6.561d0 
      Exxh2(2) = -4.061d0 
      Exxh1(3) = -4.061d0 
      Exxh2(3) = -3.811d0 
      Exxh1(4) = -3.811d0 
      Exxh2(4) = -2.761d0 
      Exxh1(5) = -2.761d0 
      Exxh2(5) =  0.00d0                         ! VBM
      Exxe1(1) =   1.420d0                       ! CBM
      Exxe2(1) =   3.340d0  
      Exxe1(2) =   3.340d0  
      Exxe2(2) =   5.370d0  
      Exxe1(3) =   5.370d0  
      Exxe2(3) =   7.8946d0  
      Exxm1(1) = EVBM    
      Exxm2(1) = 0.06762d0
      Exxm1(2) = 0.06762d0
      Exxm2(2) = 0.14200d0
      Exxm1(3) = 0.14200d0
      Exxm2(3) = 0.71000d0
      Exxm1(4) = 0.71000d0
      Exxm2(4) = 1.35238d0
      Exxm1(5) = 1.35238d0
      Exxm2(5) = ECBM     
       end subroutine set_limits_GaAs







     subroutine calc_reciprocal_param
      alat = alat*BohrA                         ! convert to A
      print *,'alat=',alat
      a2p = 2.d0*pi/alat 
      b1 = b1*a2p
      b2 = b2*a2p
      b3 = b3*a2p
      print *,'reciprocal vectors'
      print 1,b1
      print 1,b2
      print 1,b3
 1    format(3F15.4)
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
        print *,'calc_DOS_Mtot:'
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
        if(L_debug) then
         print *,'k kx ky kr'
         do k=1,Nk
          kr = dsqrt(kp(1,k)**2+kp(2,k)**2)
          print 4,k,kp(1,k),kp(2,k),kr
         enddo
        endif
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
      call read_pdos_1(Efi3,PDOS3,3)       ! PDOS of 3d layer of interface
      call read_pdos_1(Efi4,PDOS4,4)       ! PDOS of 4th layer of interface
      call read_pdos_0                     ! PDOS of the surface
      if(L_debug) print *,'N_DOS_M=',N_DOS_M
      do k=1,Nk
       do j=1,N_DOS_M
        Ef_DOS_M(j) = Efi3(j)                  ! copy PDOS3
        DOS_M(j,k) = PDOS3(j,k)
       enddo
      enddo
      call calc_DOS_Mtot(PDOS3)                                     ! calculate integrated DOS_M of interfacial layer
      V_D0  = V_D0*BohrA**3          ! convert to A^3
      Surface = V_D0/Lz             ! surface of the cell
      V_D0 = V_D0/Lz*Lz_int         ! volume of the interfacial layer for PDOS (for MIGS)
      print *,'Surface=',Surface
      print *,'Volume of the cell for interface calculation'
      print *,'V_D0=',V_D0,' A'
      call open_file(2,'dos_bulk_.dat')     ! DOS of bulk SC
      read(2,*)
      read(2,*)
      read(2,2) N_DOS_SC
      do j=1,N_DOS_SC
       read(2,*) Ef_DOS_SC(j),DOS_SC(j)
       if(DOS_SC(j) < 0.d0) then
        DOS_SC(j) = 0.d0
       endif
      enddo
      close(unit=2)
      print *,'read_PDOS: read ',N_DOS_SC,' points of bulk'
      EVBM    =  0.00d0                                                     ! VBM
      ECBM    =  1.42d0                                                     ! CBM
!      V_DSC = 304.8114d0*BohrA**3                                           ! volume of semiconductor in A^3 primitive cell
      V_DSC = V_DSC*BohrA**3                                           ! volume of semiconductor in A^3 primitive cell
      print *,'Volume of the cell for primitive cell of GaAs'
      print *,'V_DSC=',V_DSC,' A'
 1    format(15x,I4)
 2    format(19x,I4)
     end subroutine read_PDOS






     subroutine read_pdos_1(Efi,PDOS,ilayer)                 ! read PDOS for layer 3(4)
      real(8)           :: PDOS(Nptm,Nk),Efi(Nptm)
      integer           :: ilayer
      integer           :: k,j
      print *,'read_PDOS:'
      call open_file(2,'kpdos_int_'//trim(adjustl(fstr(ilayer)))//'.dat')     ! PDOS of interface (E,k)
      do k=1,Nk                             ! read Nk points
       read(2,*)
       read(2,*)
       read(2,*)
       if(k==1) then
        read(2,1) N_DOS_M
        if(L_debug) print *,'read N_DOS_M =',N_DOS_M
       else
        read(2,*)
       endif    
       do j=1,N_DOS_M
        read(2,*) Efi(j),PDOS(j,k)
       enddo
      enddo
      print *,'read_PDOS: read ',N_DOS_M,' points of interface and ',Nk,' k-points'
      close(unit=2)
 1    format(15x,I4)
     end subroutine read_pdos_1




     subroutine read_pdos_0                    ! read PDOS for surface (1st layer + metal surface)
      integer           :: j
      print *,'read_pdos_0:'
      call open_file(2,'DOStot.dat')    
       read(2,*)
       read(2,*)
       read(2,*)
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
        print *,' Open file ',name
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
     print *
     print *,'    E     polarization'
     open(unit=11,file='polarization.dat')
     read(11,*) Npol
     if(L_debug) print *,'Npol=',Npol
     allocate(Epol(Npol))
     allocate(dpol(Npol))
     do i=1,Npol
      read(11,1) Epol(i),dpol(i)
      if(L_super_debug) print 1,Epol(i),dpol(i)
     enddo
     close(unit=11)
     kappa = 0.d0
     do i=2,Npol
      kappa = kappa + (dpol(i)/Epol(i))
     enddo
     kappa = 1.d0/dfloat(Npol)*kappa
     print *,'read_pol:   kappa=',kappa
     print *,'read_pol:   kappa/e0=',kappa/e0
 1   format(F11.6,F12.7)
     er = kappa/e0 + 1.d0
    end subroutine read_pol 




  end module SBParameters













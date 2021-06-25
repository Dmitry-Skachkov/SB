

!      Program Schottky Barrier to calculate Schottky contact parameters from first-principles
!
! (C)  Copyright 2021  Center for Molecular Magnetic Quantum Materials (M2QM),
!                      Quantum Theory Project (QTP),
!                      Department of Physics, University of Florida, Gainesville, FL, USA 32611
!                      https://efrc.ufl.edu/


  Module SBSpline_functions
   use SBParameters
   use SBMathLibrary
   implicit none
   real*8, allocatable, dimension(:,:) :: bspl2,cspl2,dspl2                        ! spline coefficients for DOS_M(E,k)
   real*8, allocatable, dimension(:)   :: bspl2t,cspl2t,dspl2t                     ! spline coefficients for DOS_Mtot(E)
   real*8, allocatable, dimension(:)   :: bspl20,cspl20,dspl20                     ! spline coefficients for DOS0(E)
   real*8, allocatable, dimension(:)   :: bspl3,cspl3,dspl3                        ! spline coefficients for DOS_SC
   real*8, allocatable, dimension(:)   :: bspl13,cspl13,dspl13                     ! PDOS3
   real*8, allocatable, dimension(:)   :: bspl14,cspl14,dspl14                     ! PDOS4
   real(8), dimension(Nz)              :: bspl4,cspl4,dspl4                        ! spline coefficients for Vel (electrostatic potential)
   real(8), dimension(Nz)              :: bspl5,cspl5,dspl5                        ! spline coefficients for po (charge density)
   real(8), dimension(Nz)              :: bspl7,cspl7,dspl7                        ! spline coefficients for El_f2
   real(8), dimension(Nz)              :: bspl15,cspl15,dspl15                     ! spline coefficients for po_MIGS
   real(8), dimension(Nz)              :: bspl16,cspl16,dspl16                     ! spline coefficients for po_e
   real(8), dimension(Nz)              :: bspl17,cspl17,dspl17                     ! spline coefficients for po_h
   real*8, allocatable, dimension(:,:) :: bspl21,cspl21,dspl21                     ! spline coefficients for ImKL1
   real*8, allocatable, dimension(:,:) :: bspl22,cspl22,dspl22                     ! spline coefficients for ImKH1
  contains






   subroutine set_initial_po_V
    integer            :: i
    real(8)            :: r
    real(8)            :: er11                      ! dielectric constant
    real(8)            :: a_init
    er11 = er
    if(Calc=='s') then                              ! start new calculation
     a_init = 1.d0/za
     V0 = -(EFermi1 - CNL) 
     po0 = -V0*a_init**2*e0*er
     do i=1,Nz
      po_new(i) = po0*dexp(-a_init*Zz(i))
      V_eln(i) = V0*dexp(-a_init*Zz(i))
     enddo
     Sig = -po0/a_init
    elseif(Calc=='c') then                           ! continue calculations
     call open_file(1,'filen.dat')
      read(1,*) filen                                ! read number of last iteration
     close(unit=1)
     call open_file(2,'po_'//trim(adjustl(fstr(filen)))//'.dat')
     call open_file(3,'elpot_'//trim(adjustl(fstr(filen)))//'.dat')
     do i=1,Nz
      read(2,*) Zz(i),r,r,r,po_new(i)
      read(3,*) Zz(i),V_eln(i)
     enddo
     do i=1,Nz
      V_eln(i) = -V_eln(i)
     enddo
     close(unit=2)
     close(unit=3)
    endif
    do i=1,Nz
     V_el0(i) = V_eln(i)
     po_0(i) = po_new(i)
    enddo
 2  format(F17.5,5E21.12e3)
   end subroutine set_initial_po_V








    subroutine spline_start_1
      integer           :: k     
      print *
      print *
      print *
      print *,'spline_start'
      print *,'allocate'
      print *,'Nk=',Nk
      print *,'N_DOS_M=',N_DOS_M
      print *,'N_DOS_SC=',N_DOS_SC
      print *,'filen=',filen
      allocate(bspl2(N_DOS_M,Nk),cspl2(N_DOS_M,Nk),dspl2(N_DOS_M,Nk))            ! DOS_M(E,k)
      allocate(bspl2t(N_DOS_M),cspl2t(N_DOS_M),dspl2t(N_DOS_M))                  ! DOS_Mtot(E)
      allocate(bspl20(N_DOS_M),cspl20(N_DOS_M),dspl20(N_DOS_M))                  ! DOS0(E)
      allocate(bspl3(N_DOS_SC),cspl3(N_DOS_SC),dspl3(N_DOS_SC))
      allocate(bspl21(Npt1,Nk),cspl21(Npt1,Nk),dspl21(Npt1,Nk))                     ! CBS
      allocate(bspl22(Npt1,Nk),cspl22(Npt1,Nk),dspl22(Npt1,Nk))                     ! CBS
      allocate(bspl13(N_DOS_M),cspl13(N_DOS_M),dspl13(N_DOS_M))                 ! CBS    k=1 Gamma for Heavy
      allocate(bspl14(N_DOS_M),cspl14(N_DOS_M),dspl14(N_DOS_M))                 ! CBS    k=1 Gamma for Heavy
      print *,'1 N_DOS_M=',N_DOS_M
      call spline(Efi3(1:N_DOS_M), PDOS3(1:N_DOS_M,1),bspl13(1:N_DOS_M),cspl13(1:N_DOS_M),dspl13(1:N_DOS_M),N_DOS_M)       ! calculate spline coefficients for PDOS3
      print *,'2 N_DOS_M=',N_DOS_M
      call spline(Efi4(1:N_DOS_M), PDOS4(1:N_DOS_M,1),bspl14(1:N_DOS_M),cspl14(1:N_DOS_M),dspl14(1:N_DOS_M),N_DOS_M)       ! calculate spline coefficients for PDOS4
      print *,'3 N_DOS_M=',N_DOS_M
      do k=1,Nk
       print *,'k=',k 
       call spline(Ef1(1:Npt1),ImKL1(1:Npt1,k),bspl21(1:Npt1,k),cspl21(1:Npt1,k),dspl21(1:Npt1,k),Npt1)          ! calculate spline coefficients 
       call spline(Ef1(1:Npt1),ImKH1(1:Npt1,k),bspl22(1:Npt1,k),cspl22(1:Npt1,k),dspl22(1:Npt1,k),Npt1)          ! calculate spline coefficients 
      enddo
      print *,'DOS_M'
      do k=1,Nk
       call spline(Ef_DOS_M,DOS_M(1:N_DOS_M,k),bspl2(1:N_DOS_M,k),cspl2(1:N_DOS_M,k),dspl2(1:N_DOS_M,k),N_DOS_M)    ! calculate spline coefficients for DOS_M
      enddo
      call spline(Ef_DOS_M, DOS_Mtot,bspl2t,cspl2t,dspl2t,N_DOS_M)                ! calculate spline coefficients 
      call spline(Efi0, DOS0,bspl20,cspl20,dspl20,N_DOS_M)                        ! calculate spline coefficients 
      call spline(Ef_DOS_SC,DOS_SC, bspl3, cspl3, dspl3, N_DOS_SC)                ! calculate spline coefficients 
      call calc_z_mesh
    end subroutine spline_start_1






    subroutine spline_start_2
     call spline(Zz,V_eln,bspl4,cspl4,dspl4,Nz)                                  ! calculate spline coefficients for Vel
     call spline(Zz,po_new,bspl5,cspl5,dspl5,Nz)                                 ! calculate spline coefficients for po
    end subroutine spline_start_2







     subroutine calc_z_mesh                                     ! logarithmic z-mesh points
      real(8)          :: dz
      integer          :: i
      real(8)          :: pz,dp
      dp = dfloat(Nz-1)*0.01d0/dfloat(Nz-1)
      pz = -1.d0
      print *
      print *,'Z mesh:'
      print *,'dp=',dp
      Zz(1) = 0.d0
      do i=4,Nz
       Zz(i) = 10.d0**pz
       pz = pz + dp
      enddo
      Zz(1) = 1.D-4
      Zz(2) = 1.D-3
      Zz(3) = 1.D-2
      do i=1,Nz
       print 2,i,Zz(i)
      enddo
1     format(2F17.5)
2     format(I5,F12.3)
     end subroutine calc_z_mesh





   real(8) function Vels(z)
    real(8)             :: z
    Vels = ispline(z,Zz,V_eln,bspl4,cspl4,dspl4,Nz)
   end function Vels





   real(8) function pos(z)
    real(8)             :: z
    pos = ispline(z,Zz,po_new,bspl5,cspl5,dspl5,Nz)
   end function pos



   real(8) function poMIGSs(z)
    real(8)             :: z
    poMIGSs = ispline(z,Zz,po_MIGS,bspl15,cspl15,dspl15,Nz)
   end function poMIGSs



   real(8) function po_es(z)
    real(8)             :: z
    po_es = ispline(z,Zz,po_e,bspl16,cspl16,dspl16,Nz)
    po_es = po_es - po00_e0
   end function po_es



   real(8) function po_hs(z)
    real(8)             :: z
    po_hs = ispline(z,Zz,po_h,bspl17,cspl17,dspl17,Nz)
    po_hs = po_hs - po00_h0
   end function po_hs




   real(8) function dpols(z)                               ! polarization
    real(8)             :: z
    real(8)             :: Ex                              ! electric field at point z
    integer             :: i
     Ex = El_f2s(z)
     dpols = kappa*Ex
   end function dpols



   real(8) function El_f2s(z)
    real(8)             :: z
    El_f2s = ispline(z,Zz,El_f2,bspl7,cspl7,dspl7,Nz)
   end function El_f2s






   real(8) function ImKs2(E,k)
    real(8)               :: E                                                     ! Energy vs EFermi
    integer               :: k
    if(E >= Emin1 .and. E <= Emax1) then
     ImKs2 = ispline(E,Ef1(1:Npt1),ImKL1(1:Npt1,k),bspl21(1:Npt1,k),cspl21(1:Npt1,k),dspl21(1:Npt1,k),Npt1) 
    else
     ImKs2 = 1.d5
    endif
   end function ImKs2




   real(8) function ImKHs2(E,k)                           ! heavy holes for Gamma only GaAs
    real(8)               :: E                                                     ! Energy vs EFermi
    integer               :: k
    if(E >= Emin1 .and. E <= Emax1) then
     ImKHs2 = ispline(E,Ef1(1:Npt1),ImKH1(1:Npt1,k),bspl22(1:Npt1,k),cspl22(1:Npt1,k),dspl22(1:Npt1,k),Npt1) 
    else
     ImKHs2 = 1.d5
    endif
   end function ImKHs2



   real(8) function DOS_Ms(E,k)                               ! DOS_M(E,k) interfacial DOS separated by kx,ky
    real(8)               :: E                                ! Energy vs EFermi
    integer               :: k                                ! k-point
    if(E < Ef_DOS_M(1) .or. E > Ef_DOS_M(N_DOS_M)) then
     DOS_Ms = 0.d0
    else
     DOS_Ms = ispline(E,Ef_DOS_M,DOS_M(1:N_DOS_M,k),bspl2(1:N_DOS_M,k),cspl2(1:N_DOS_M,k),dspl2(1:N_DOS_M,k),N_DOS_M) 
    endif
    if(DOS_Ms < 0.d0) then
     DOS_Ms = 0.d0
    endif
   end function DOS_Ms




   real(8) function DOS_Mtots(E)                ! DOS_Mtot     integrated over kx,ky DOS of interfacial layer  
    real(8)               :: E                                                     ! Energy vs VBM
    if(E < Ef_DOS_M(1) .or. E > Ef_DOS_M(N_DOS_M)) then
     DOS_Mtots = 0.d0
    else
     DOS_Mtots = ispline(E,Ef_DOS_M,DOS_Mtot,bspl2t,cspl2t,dspl2t,N_DOS_M) 
    endif
    if(DOS_Mtots < 0.d0) then
     DOS_Mtots = 0.d0
    endif
   end function DOS_Mtots



   real(8) function DOS0s(E)                ! DOS_Mtot     integrated over kx,ky DOS of interfacial layer  
    real(8)               :: E                                                     ! Energy vs VBM
    if(E < Efi0(1) .or. E > Efi0(N_DOS_M)) then
     DOS0s = 0.d0
    else
     DOS0s = ispline(E,Efi0,DOS0,bspl20,cspl20,dspl20,N_DOS_M) 
    endif
    if(DOS0s < 0.d0) then
     DOS0s = 0.d0
    endif
   end function DOS0s




   real(8) function DOS_SCs(E)
    real(8)               :: E                                                     ! Energy vs VBM
    if(E < Ef_DOS_SC(1) .or. E > Ef_DOS_SC(N_DOS_SC)) then
     DOS_SCs = 0.d0
    else
     DOS_SCs = ispline(E,Ef_DOS_SC,DOS_SC,bspl3,cspl3,dspl3,N_DOS_SC) 
    endif
    if(DOS_SCs < 0.d0) then
     DOS_SCs = 0.d0
    endif
   end function DOS_SCs












      subroutine calc_CBS_limits
   !    real(8)                  :: X1,X2,X3
   !    integer                  :: L
   !    integer                  :: k
   !    do k=1,Nk
         Emin1 = Ef1(1)
         Emax1 = Ef1(Npt1)
      ! enddo
       if(L_debug) then
       print *
       print *,'calc_CBS_limits'
       print 2
      ! do k=1,Nk
        print 1,Emin1,Emax1
      ! enddo
       endif
 1     format(2E12.4)
 2     format('     Emin1           Emax1 ')
      end subroutine calc_CBS_limits
 








     subroutine calc_DL_DH(pds3,pds4,ImKz_Hx,ImKz_Lx,DL,DH,DLDH)
      real(8)     :: pds3,pds4
      real(8)     :: ImKz_Hx,ImKz_Lx
      real(8)     :: DL,DH,DLDH
      real(8)     :: a,b,c,d
      a = dexp(-2.d0*ImKz_Lx*z3)
      b = dexp(-2.d0*ImKz_Hx*z3)
      c = dexp(-2.d0*ImKz_Lx*z4)
      d = dexp(-2.d0*ImKz_Hx*z4)
       if(dabs(pds4*a - pds3*c) > 1.d-6) then
        DLDH = (pds3*d - pds4*b)/(pds4*a - pds3*c)
       else
        DLDH = 1.d0
       endif
       if(DLDH > 1.d0) DLDH = 1.d0
       if(DLDH < 0.d0) then
        DLDH = 0.d0
       endif
     end subroutine calc_DL_DH




     subroutine separate(Ex,ImKz_Hx,ImKz_Lx,DLDH)
      real(8)       :: Ex
      real(8)       :: pds3,pds4
      real(8)       :: ImKz_Hx,ImKz_Lx
      real(8)       :: DL,DH,DLDH
      pds3 = ispline(Ex,Efi3(1:N_DOS_M),PDOS3(1:N_DOS_M,1),bspl13(1:N_DOS_M),cspl13(1:N_DOS_M),dspl13(1:N_DOS_M),N_DOS_M)
      pds4 = ispline(Ex,Efi4(1:N_DOS_M),PDOS4(1:N_DOS_M,1),bspl14(1:N_DOS_M),cspl14(1:N_DOS_M),dspl14(1:N_DOS_M),N_DOS_M)
      call calc_DL_DH(pds3,pds4,ImKz_Hx,ImKz_Lx,DL,DH,DLDH)          
     end subroutine separate





   end module SBSpline_functions





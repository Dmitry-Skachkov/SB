


!      Program Schottky Barrier to calculate Schottky contact parameters from first-principles
!
! Copyright (C) 2021   Center for Molecular Magnetic Quantum Materials (M2QM),
!                      Quantum Theory Project (QTP),
!                      Department of Physics, University of Florida, Gainesville, FL 32611, USA 
!                      https://efrc.ufl.edu/




     Module SBCalc_z                                           ! calculate density and V on z-mesh points and do scf cycle
      use SBParameters
      use SBCalc_i
      use SBSpline_functions
      implicit none
      integer          :: Nitscf20
      logical          :: L_pre
     contains





     subroutine calc_elpot                                     ! calculate Schottky potential
      integer          :: i,j
      real(8)          :: delta_Vi
      real(8)          :: V1(Nz-1),V2(Nz-1)
      real(8)          :: Vel11
      real(8)          :: Vtest(3)
      do i=1,Nz                                   
       V_el1(i) = Vel(Zz(i))
      enddo
      do i=2,Nz-1                                              ! calculate derivative
       El_f(i) = -(V_el1(i+1)-V_el1(i-1))/(Zz(i+1)-Zz(i-1))    ! V/Ang             
      enddo
      El_f(Nz) = 0.d0
      El_f(1) = El_f(2) - (El_f(3) - El_f(2))*(Zz(2)-Zz(1))/(Zz(3)-Zz(2))    
      delta_V = 0.d0
      do i=1,Nz
       delta_Vi = dabs(V_el1(i)-V_el0(i))
       delta_V = delta_V + delta_Vi**2
      enddo
      delta_V = dsqrt(delta_V)/dfloat(Nz)
     end subroutine calc_elpot








     subroutine calc_po1                                  ! calculate charge density on z-mesh points
      integer          :: i
      real(8)          :: po2,po3,po4     
      poh_max = 0.d0
      poe_max = 0.d0
      poMe_max = 0.d0
      delta_po = 0.d0
       do i=1,Nz        
        if(L_super_debug) print *,'calc_po1: i=',i 
        po2 = poh(Zz(i))
        po_h(i) = po2
        poh_max = dmax1(poh_max,po2)
        po3 = poe(Zz(i))
        po_e(i) = po3
        poe_max = dmin1(poe_max,po3)
        po4 = poMIGS(Zz(i))                  
        po_MIGS(i) = po4
        poMe_max = dmin1(poMe_max,po4)
        po1(i) = po2 + po3 + po4 - po00                    ! h + e + MIGS + doping
        if(po_new(i)> 1.d-15) then
         delta_po = delta_po + (po1(i)-po_0(i))**2
        endif
       enddo
       delta_po = delta_po/dfloat(Nz)
       po_MIGS(1) = po_MIGS(2)
     end subroutine calc_po1








     subroutine mixing_po
      integer         :: i
      real(8)         :: po2,po3,po4
      if(L_super_debug) print *,'mixing_po:'
      do i=1,Nz                                                      ! mixing new and old densities
       po_new(i) = (1.d0-alfa)*po_0(i) + alfa*po1(i)
      enddo
      do i=1,Nz
       po_0(i) = po_new(i)
      enddo
      call spline (Zz,po_new,bspl5,cspl5,dspl5,Nz)                    ! calculate spline coefficients for po_new
      call spline (Zz,po_MIGS,bspl15,cspl15,dspl15,Nz)                ! calculate spline coefficients for po_new
      call spline (Zz,po_e,bspl16,cspl16,dspl16,Nz)                   ! calculate spline coefficients for po_e
      call spline (Zz,po_h,bspl17,cspl17,dspl17,Nz)                   ! calculate spline coefficients for po_h
     end subroutine mixing_po





     subroutine mixing_V
      integer                :: i
      real(8)                :: dV(Nz),dV1(Nz),dV0,ddV0,ddV01,dVn0,dF0
      real(8), dimension(10) :: bsplx,csplx,dsplx
      real(8)                :: z01
      real(8)                :: Fx,F1,F2
      real(8)                :: V_els(Nz)    ! add for test
      alfa_V = 1.d0
      alfa_Sig = 1.d0
      do i=1,Nz
       V_els(i) = V_el1(i)           ! just copy
       V_el1(i) = V_els(i)    ! + alfa_dipole*Vdip(i)      ! sum 
      enddo
      do i=3,Nz-1                                    ! calculate derivative
       El_f2(i) = -(V_el1(i+1)-V_el1(i-1))/(Zz(i+1)-Zz(i-1))            ! V/Ang             
      enddo
      El_f2(Nz) = 0.d0
      El_f2(2) = El_f2(3) - (El_f2(4) - El_f2(3))*(Zz(3)-Zz(3))/(Zz(4)-Zz(3))
      El_f2(1) = El_f2(2) - (El_f2(3) - El_f2(2))*(Zz(2)-Zz(1))/(Zz(3)-Zz(2))
      z01 = 100.d0
      dV0 = dV(1)
      call spline(Zz,dV,bsplx,csplx,dsplx,10)                        ! make cubic spline for 3 points in order to calculate derivative at z=0
      ddV0 = bsplx(1)                                                ! V(z) = V(i) + b(i)*(z-z(i)) + c(i)*(z-z(i))**2 + d(i)*(z-z(i))**3
                                                                     ! dV/dz(z=0) = b(1)
  !    F1 = dV0*(1.d0-alfa_V)
  !    F2 = 3.d0*ddV0*(1.d0-alfa_Sig)/10.d0     ! need to be removed
       F1 = 0.d0
       F2 = 0.d0
      dV1(1:Nz) = 0.d0
      do i=1,Nz                                                      ! mixing new and old V with correction
!       Fx = F1*dexp(-((Zz(i))/z01)**2) + F2*Bessel_j1(10.d0*Zz(i))
!       if(dabs(Fx) > 10.d0) then
!        print *
!        print *,'Problems with the correction scheme: '
!        print *,'Fx =',Fx
!        stop
!       endif
!       dV1(i) = dV(i) - Fx 
       Fx = 0.d0                
       V_eln(i) = V_el1(i) - Fx
      enddo
      call spline(Zz,V_eln,bsplx,csplx,dsplx,3)                        ! make cubic spline for 3 points in order to calculate derivative at z=0
      dVn0 = bsplx(1)*(er*e0)                                               ! V(z) = V(i) + b(i)*(z-z(i)) + c(i)*(z-z(i))**2 + d(i)*(z-z(i))**3
      do i=1,Nz
       Fx = F1*dexp(-((Zz(i))/z01)**2) + F2*Bessel_j1(Zz(i))
      enddo
      call spline (Zz,V_eln,bspl4,cspl4,dspl4,Nz)                      ! calculate spline coefficients for Vel
   !   call spline (Zz,El_f2,bspl7,cspl7,dspl7,Nz)                      ! calculate spline coefficients for El_f2
      do i=1,Nz
       V_el0(i) = V_eln(i)
      enddo
 2    format(F17.5,F18.10,E17.5e3)
 3    format('        z             V_el1(i)          V_el0(i)           dV(i)             ddV                F               dV1(i)           ddV1             V_eln         dV_eln/dz       dV/dz*(er*e0)')
 1    format(F12.4,8F18.8,2E15.5)
 7    format(F12.4,3F18.8,E18.5e3,F18.8,2E18.5e3,E15.5e3)
12    format(F17.5,F18.10,F18.10,F18.10,E17.5e3,E17.5e3)
33    format(I5)
     end subroutine mixing_V






     real(8) function Bessel_j1(z)
      real(8)     :: z
      if(dabs(z) < 1.d-20) then
       Bessel_j1 = 0.d0
      else
       Bessel_j1 = dsin(z)/z**2 - dcos(z)/z
      endif
     end function Bessel_j1







       subroutine calc_scf                                               ! find SCF solution for the system M-SC
        real(8)       :: zz1,zz2
        integer       :: is,is0,is1
        logical       :: Lsuc1,Lsuc2
        real(8)       :: diffV
        filen = 0
        print 5
        print 4
        call set_limits_zz1(zz1,zz2)                                     ! set initial limits for searching solution
        call set_initial_deltaE                                          ! set first approximation to deltaE
        do is = 1,Nitscf                                                 ! scf cycle over deltaE
         Nitscf20 = 1
         L_pre = .true.
         do is0 = 1,Nitscf0                                              ! pre cycle in order to find approximate solution
          call calc_look_po_z(is0,zz1,zz2,Lsuc1,Lsuc2)                  
          if(Lsuc1 .and. Lsuc2) then                                     ! if success
           print 1,is0,filen
           exit
          elseif((.not.Lsuc1 .or. .not.Lsuc2) .and. (is0==Nitscf0)) then ! not success
           call print_scf_error(is0,zz1,zz2,Lsuc1,Lsuc2)
           stop
          endif
          call reset_zz(zz1,zz2,Lsuc1,Lsuc2)
         enddo
         call set_new_limits(zz1,zz2)                                    ! narrow limits for next accurate scf
         Nitscf20 = Nitscf2
         L_pre = .false.
         do is1 = 1,Nitscf1                                              ! post scf cycle in order to find accurate solution
          call calc_look_po_z(is1,zz1,zz2,Lsuc1,Lsuc2)                  
          if(Lsuc1 .and. Lsuc2) then
           print 3,is1,filen
           exit
          elseif((.not.Lsuc1 .or. .not.Lsuc2) .and. (is1==Nitscf1)) then ! not success
           call print_scf_error(is1,zz1,zz2,Lsuc1,Lsuc2)
      !     print 2
           stop
          endif
          call reset_zz2(zz1,zz2,Lsuc1,Lsuc2)
         enddo
         call calc_charges                                               ! calculate charge on the interface
         call calc_deltaE                                                ! calculate filling level of the surface
        enddo
 1      format(' We have found the approximate solution at ',I3,' cycle of pre scf cycle using ',I4,' iterations'/)
 2      format(/' Please make limits in set_new_limits wider') 
 3      format(' We have found the accurate solution at ',I3,' cycle of post scf cycle using ',I4,' iterations'/)
 4      format(/'  it    -eV(0)      Sig(po)       dV/dz(er*e0)     poh_m          poe_m        poMe_m          po(0)      delta_po     delta_V')
 5      format(//' Start SCF cycle')
       end subroutine calc_scf





       subroutine print_scf_error(Natmpt,za1,za2,Lsuc1,Lsuc2)
        real(8)        :: za1,za2   
        integer        :: Natmpt
        logical        :: Lsuc1,Lsuc2
        print 3,Natmpt
        if(.not.Lsuc1) print 10,za1
        if(.not.Lsuc2) print 11,za2
        print 9
        print 6,EFermi1
        if(L_p_type) print 7,dabs(EFermi1-EFermi_00)
        if(L_n_type) print 8,dabs(EFermi1-EFermi_00)
        print 5,po00*1.d24
 3      format(/' We made ',I4,' attempts to find the solution, unfortunately it is outside of the region')
 5      format(' Doping concentration is ',1p,E12.3,' cm-3')
 6      format(' Fermi level =',F12.4,' eV')
 7      format(' p-type doped semiconductor with Fermi level of ',F12.4,' eV below the charge neutrality level for bulk')
 8      format(' n-type doped semiconductor with Fermi level of ',F12.4,' eV higher the charge neutrality level for bulk')
 9      format(' The parameters of the system are:')
10      format(' za < ',F20.4,' A')
11      format(' za > ',F20.4,' A')
       end subroutine print_scf_error






       subroutine reset_zz(zz1,zz2,Lsuc1,Lsuc2)                  ! increase searching range
        real(8)       :: zz1,zz2
        logical       :: Lsuc1,Lsuc2
        if(.not.Lsuc1) then
         zz1 = 0.1d0*zz1
        elseif(.not.Lsuc2) then
         zz2 = 10.d0*zz2
        endif
       end subroutine reset_zz




       subroutine reset_zz2(zz1,zz2,Lsuc1,Lsuc2)                 ! increase searching range
        real(8)       :: zz1,zz2
        logical       :: Lsuc1,Lsuc2
        if(.not.Lsuc1) then
         zz1 = 0.9d0*zz1
        elseif(.not.Lsuc2) then
         zz2 = 1.1d0*zz2
        endif
       end subroutine reset_zz2





       subroutine set_limits_zz1(zz1,zz2)
        real(8)     :: zz1,zz2
        real(8)     :: x,zz
        x = (EFermi1 - EFermi_00)/(ECBM - EFermi_00)
        zz = 1.d0/(0.00000001d0*x + 0.000001d0*x**2 + 0.00001d0*x**3 + 0.0010d0*x**4)
        zz1 = 0.1d0*zz
        zz2 = 10.d0*zz
       end subroutine set_limits_zz1




       subroutine set_new_limits(zz1,zz2)                                   ! update limits for next iterations
        real(8)     :: zz1,zz2
        zz1 = za*0.9d0
        zz2 = za*1.1d0
       end subroutine set_new_limits




       subroutine set_initial_deltaE                                        ! set first approximation to deltaE
        real(8)          :: x
        if(L_n_type) then
         x = (EFermi1 - EFermi_00)/(ECBM - EFermi_00)
         dEf =  0.005d0*x**2
        elseif(L_p_type) then
         x = (EFermi1 - EFermi_00)/(EVBM - EFermi_00)
         dEf = -0.005d0*x**2
        else
         dEf = 0.d0
        endif
       end subroutine set_initial_deltaE




       subroutine calc_look_po_z(Natmpt,za1,za2,Lsuc1,Lsuc2)                ! find za on the interval [za1,za2]
        real(8)        :: za1,za2   
        real(8)        :: eps
        real(8)        :: a,b
        real(8)        :: diffV
        integer        :: Natmpt
        logical        :: Lsuc1,Lsuc2
        eps = 0.0001d0                                                       ! 0.1 meV accuracy to find the Fermi level
        a = za1
        b = za2
        Lsuc1 = .true.
        Lsuc2 = .true.
        do while (dabs(a-b) > eps)
         filen = filen + 1
         za = (a+b)/2.d0
         if(L_debug) print *,'za=',za
         call calc_diff_eV0(diffV)
         if(L_debug) print *,'diffV=',diffV
        if(EFermi1 > CNL+dEf) then
         if(diffV > 0.d0) then
          b = za
         else
          a = za
         endif
        else
         if(diffV < 0.d0) then
          b = za
         else
          a = za
         endif
        endif
        enddo
        if( (dabs(za-za1).le.0.01d0) ) then
         if(.not.Natmpt==Nitscf0) print 1
         Lsuc1 = .false.
        elseif( (dabs(za-za2).le.0.01d0) ) then
         if(.not.Natmpt==Nitscf0) print 2
         Lsuc2 = .false.
        endif
        if(L_debug) then
         print *,'filen=',filen
         print *,'found za=',za
         print *,'-eV(0)=',-V_eln(1)
         print *,'EFermi1=',EFermi1
         print *,'CNL=',CNL
         print *,'EFermi1-(CNL+dEf)=',EFermi1-(CNL+dEf)
        endif
 1      format(/'***** ERROR *****'/' za is too close to za1'/' Let us decrease the value za1 and try one more time'/)
 2      format(/'***** ERROR *****'/' za is too close to za2'/' Let us increase the value za2 and try one more time'/)
       end subroutine calc_look_po_z




                                          
        subroutine calc_diff_eV0(diffV)
         real(8)         :: diffV
         integer         :: i
         call set_initial_po_V                                                ! calculate V(z)
         call spline_start_2                                                  ! spline coeff. for V_eln and po_new
         if(.not.L_pre) call calc_deltaE                                      ! calculate filling level of the surface
         do i=1,Nitscf20
          call calc_po1                                                       ! calculate charge density po1 (po_h, po_e, and po_MIGS)
          call mixing_po                                                      ! mixing new and old density (po_new = (1-a)po0+a*po1) + spline coeff. 
          call calc_elpot                                                     ! calculate electrostatic potential using po_new
          call mixing_V                                                       ! spline coeff. for V_eln   (teper' i ne nuzhno!)
         enddo
         diffV = -V_eln(1) - (EFermi1 - (CNL+dEf))                            ! -eV(0) - (EFermi1 - (CNL+dEf))
         if(L_super_debug) then
          print *,'-V_eln(1)=',-V_eln(1)
          print *,'filling level=',dEf
          print *,'(EFermi1 - (CNL+dEf))=',(EFermi1 - (CNL+dEf))
         endif
         print 2,filen,dabs(-V_eln(1))-dabs(dEf),Sig,bspl4(1)*er*e0,poh_max*1.d24,poe_max*1.d24,poMe_max*1.d24,po_new(1)*1.d24,delta_po*1.d24,delta_V
 2       format(I4,F11.5,2E15.5,6E14.4,F12.4)
        end subroutine calc_diff_eV0




        subroutine calc_check_scf    
         real(8)      :: acc_scf
         real(8)      :: ddEf
         call calc_po1                                                        ! calculate charge density po1 (po_h, po_e, and po_MIGS)
         acc_scf = delta_po/po_new(1)
         if(acc_scf < 1.d-5) then
          L_conv = .true.
         endif
         ddEf = dabs(dEf-dEf1)
         print 1,acc_scf
         print 2,ddEf
 1       format(' accuracy of scf cycle ',1p,E12.3) 
 2       format(' accuracy reached for deltaE ',1p,E12.3)
        end subroutine calc_check_scf




        subroutine calc_charges
         call calc_totq                                                       ! calculate total charge on interface
   !      call calc_MIGS_q                                                     ! calculate equivalent charge due to MIGS
   !      call calc_h_q                                                        ! calculate equivalent charge due to holes
   !      call calc_e_q                                                        ! calculate equivalent charge due to electrons
          if(L_debug) print *,'Sig_total=',Sig
   !      print *,'SigMIGS=',SigMIGS
   !      print *,'Sig_e=',Sig_e
   !      print *,'Sig_h=',Sig_h
        end subroutine calc_charges






   end module SBCalc_z





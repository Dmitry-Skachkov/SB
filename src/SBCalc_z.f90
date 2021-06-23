


!      Program Schottky Barrier to calculate Schottky contact parameters from first-principles
!
! (C)  Copyrigth 2021, Skachkov, Zhang, Cheng, Center for Molecular Magnetic Quantum Materials (M2QM),
!                      Quantum Theory Project (QTP),
!                      Department of Physics, University of Florida, Gainesville, FL, USA 32611
!                      https://efrc.ufl.edu/




     Module SBCalc_z                                           ! calculate density and V on z-mesh points and do scf cycle
      use SBParameters
      use SBCalc_i
      use SBSpline_functions
      implicit none
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
      delta_Vm = 0.d0
      do i=1,Nz
       delta_Vi = dabs(V_el1(i)-V_el0(i))
       delta_V = delta_V + delta_Vi**2
       delta_Vm = dmax1(delta_Vm,delta_Vi)
      enddo
      delta_V = dsqrt(delta_V)
     end subroutine calc_elpot








     subroutine calc_po1
      integer          :: i
      real(8)          :: po2,po3,po4,po1Nz
      poh_max = 0.d0
      poe_max = 0.d0
      poMe_max = 0.d0
      delta_po = 0.d0
       do i=Nz,1,-1     ! going to z(1)
        po2 = poh(Zz(i))
        po_h(i) = po2
        poh_max = dmax1(poh_max,po2)
        po3 = poe(Zz(i))
        po_e(i) = po3
        poe_max = dmin1(poe_max,po3)
        po4 = poMIGS(Zz(i))                  
        po_MIGS(i) = po4
        poMe_max = dmin1(poMe_max,po4)
        po1(i) = po2 + po3 + po4 - po00 
        if(i==Nz) po1Nz = po1(i)                  ! make po1(Nz) = 0, just to correct inaccuracy of calculations
        if(i==Nz) then
         if(dabs(po1(Nz)) > 1.d-24) then
          print *,'po1(Nz)=',po1(Nz)
         endif
        endif
        po1(i) = po1(i) - po1Nz                   ! make po1(Nz) = 0
        if(po_new(i)> 1.d-15) then
         delta_po = delta_po + (po1(i)-po_0(i))**2
        endif
       enddo
       po_MIGS(1) = po_MIGS(2)
     end subroutine calc_po1












     subroutine mixing_po
      integer         :: i
      real(8)         :: po2,po3,po4
       do i=1,Nz                                                      ! mixing new and old densities
        po_new(i) = (1.d0-alfa)*po_0(i) + alfa*po1(i)
!        po_new(i) =  po1(i)                           ! do not mix
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
      call spline (Zz,El_f2,bspl7,cspl7,dspl7,Nz)                      ! calculate spline coefficients for El_f2
      do i=1,Nz
       V_el0(i) = V_eln(i)
      enddo
2     format(F17.5,F18.10,E17.5e3)
 12    format(F17.5,F18.10,F18.10,F18.10,E17.5e3,E17.5e3)
33    format(I5)
 3    format('        z             V_el1(i)          V_el0(i)           dV(i)             ddV                F               dV1(i)           ddV1             V_eln         dV_eln/dz       dV/dz*(er*e0)')
 1    format(F12.4,8F18.8,2E15.5)
 7    format(F12.4,3F18.8,E18.5e3,F18.8,2E18.5e3,E15.5e3)
     end subroutine mixing_V






     real(8) function Bessel_j1(z)
      real(8)     :: z
      if(dabs(z) < 1.d-20) then
       Bessel_j1 = 0.d0
      else
       Bessel_j1 = dsin(z)/z**2 - dcos(z)/z
      endif
     end function Bessel_j1







       subroutine calc_scf                                             ! find SCF solution for the system M-SC
        real(8)       :: zz1,zz2
        integer       :: is
        call set_limits_zz1(zz1,zz2)
        do is = 1,Nitscf
        call calc_look_po_z(zz1,zz2)                                   ! find za corresonding to -eV0
        call calc_charges                                              ! calculate charge on the interface
        call calc_deltaE                                               ! calculate filling level of the surface
        print *,'CBM=',ECBM
        print *,'EFermi1=',EFermi1
        if(L_n_type) then
         SBH = -V_eln(1) + (ECBM - EFermi1)
        elseif(L_p_type) then
         SBH =  V_eln(1) + (EFermi1 - EVBM)
        endif
        print *,'is=',is,' SBH=',SBH
        enddo
       end subroutine calc_scf






       subroutine set_limits_zz1(zz1,zz2)
        real(8)     zz1,zz2
        if(1.40d0 < EFermi1 .and. EFermi1 <= 1.42d0) then
         zz1 = 100.d0 
         zz2 = 300.d0    !5000.d0 
        elseif(1.35d0 < EFermi1 .and. EFermi1 <= 1.40d0) then
         zz1 = 100.d0 
         zz2 = 5000.d0 
        elseif(1.25d0 < EFermi1 .and. EFermi1 <= 1.35d0) then
         zz1 = 1000.d0 
         zz2 = 30000.d0 
        elseif(1.15d0 < EFermi1 .and. EFermi1 <= 1.25d0) then
         zz1 = 5000.d0 
         zz2 = 100000.d0
        elseif(0.99d0 < EFermi1 .and. EFermi1 <= 1.15d0) then
         zz1 = 5000.d0 
         zz2 = 100000.d0
        elseif(0.00d0 <= EFermi1 .and. EFermi1 < 0.20d0) then
         zz1 = 3.d0 
         zz2 = 300.d0
        elseif(0.20d0 <= EFermi1 .and. EFermi1 < 0.24d0) then
         zz1 = 10.d0 
         zz2 = 5000.d0
        elseif(0.24d0 <= EFermi1 .and. EFermi1 < 0.27d0) then
         zz1 = 2000.d0 
         zz2 = 2500.d0
        elseif(0.27d0 <= EFermi1 .and. EFermi1 < 0.33d0) then
         zz1 = 1000.d0 
         zz2 = 30000.d0
        else
         print *,'define limits for zz1,zz2'
         stop
        endif 
       end subroutine set_limits_zz1








       subroutine calc_look_po_z(za1,za2)                                     ! find za corresonding to -eV0 on the interval [za1,za2]
        real(8)        :: za1,za2   
        real(8)        :: eps
        real(8)        :: a,b
        real(8)        :: diffV
        eps = 0.0001d0                                                        ! 0.1 meV accuracy to find the Fermi level
        a = za1
        b = za2
         filen = 0
        do while (dabs(a-b) > eps)
         filen = filen + 1
         za = (a+b)/2.d0
         print *,'za=',za
         call calc_diff_eV0(diffV)
         print *,'diffV=',diffV
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
        if( (dabs(za-za1).le.0.01d0) .or. &
            (dabs(za-za2).le.0.01d0) ) then
         print 1
         print *,'za=',za
         print *,'za1=',za1
         print *,'za2=',za2
         stop
        endif
        print *,'filen=',filen
        print *,'found za=',za
        print *,'-eV(0)=',-V_eln(1)
        print *,'EFermi1=',EFermi1
        print *,'CNL=',CNL
        print *,'EFermi1-(CNL+dEf)=',EFermi1-(CNL+dEf)
 1      format(' za is too close to za1 or za2')
       end subroutine calc_look_po_z




                                          
        subroutine calc_diff_eV0(diffV)
         real(8)         :: diffV
         call set_initial_po_V                                                ! calculate V(z)
         call spline_start_2                                                  ! spline coeff. for V_eln and po_new
         filen = filen + 1
         call calc_po1                                                        ! calculate charge density po1 (po_h, po_e, and po_MIGS)
         call mixing_po                                                       ! mixing new and old density (po_new = (1-a)po0+a*po1) + spline coeff. 
         call calc_elpot                                                      ! calculate electrostatic potential using po_new
         call mixing_V                                                        ! spline coeff. for V_eln   (teper' i ne nuzhno!)
         diffV = -V_eln(1) - (EFermi1 - (CNL+dEf))                            ! -eV(0) - (EFermi1 - (CNL+dEf))
         print *,'-V_eln(1)=',-V_eln(1)
         print *,'filling level=',dEf
         print *,'(EFermi1 - (CNL+dEf))=',(EFermi1 - (CNL+dEf))
        end subroutine calc_diff_eV0



        subroutine calc_charges
         call calc_totq                                                       ! calculate total charge on interface
   !      call calc_MIGS_q                                                     ! calculate equivalent charge due to MIGS
   !      call calc_h_q                                                        ! calculate equivalent charge due to holes
   !      call calc_e_q                                                        ! calculate equivalent charge due to electrons
         print *,'Sig_total=',Sig
   !      print *,'SigMIGS=',SigMIGS
   !      print *,'Sig_e=',Sig_e
   !      print *,'Sig_h=',Sig_h
        end subroutine calc_charges






   end module SBCalc_z





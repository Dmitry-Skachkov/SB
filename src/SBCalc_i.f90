


!      Program Schottky Barrier to calculate Schottky contact parameters from first-principles
!
! (C)  Copyright 2021, Skachkov, Zhang, Cheng, Center for Molecular Magnetic Quantum Materials (M2QM),
!                      Quantum Theory Project (QTP),
!                      Department of Physics, University of Florida, Gainesville, FL, USA 32611
!                      https://efrc.ufl.edu/



     Module SBCalc_i                                           ! integration, calculate charge density and potential
      use SBParameters
      use SBSpline_functions
      use SBMathLibrary
      implicit none
      real(8)                :: eVz
      real(8)                :: eV0
      real(8)                :: EFermi111
      real(8)                :: dEf1
     contains




     real(8) function Vel(z)                                   ! electrostatic potential in V
      real(8)            :: z
      real(8)            :: eps
      real(8)            :: intpo2,intpo21,intpo22,intpo23,intpo24,intpo25,intpo26
      eps = 1.d-10
      zp = z
      if(z < 10.d0) then                                       
       intpo21 = 0.d0
       intpo22 = 0.d0
       intpo23 = 0.d0
       intpo24 = 0.d0
       intpo25 = 0.d0
       intpo26 = 0.d0
       if(Lsc > 10.d0) then
        call QSL3D(intpo21,z,10.d0,F6,eps)
        if(Lsc > 100.d0) then
         call QSL3D(intpo22,10.d0,100.d0,F6,eps)
         if(Lsc > 1000.d0) then
          call QSL3D(intpo23,100.d0,1000.d0,F6,eps)
          if(Lsc > 10000.d0) then
           call QSL3D(intpo24,1000.d0,10000.d0,F6,eps)
           if(Lsc > 100000.d0) then
            call QSL3D(intpo25,10000.d0,100000.d0,F6,eps)
            call QSL3D(intpo26,100000.d0,Zz(Nz),F6,eps)
           else
            call QSL3D(intpo25,10000.d0,Lsc,F6,eps)
           endif
          else
           call QSL3D(intpo24,1000.d0,Lsc,F6,eps)
          endif
         else
          call QSL3D(intpo23,100.d0,Lsc,F6,eps)
         endif
        else
         call QSL3D(intpo22,10.d0,Lsc,F6,eps)
        endif
       else
        call QSL3D(intpo21,z,Lsc,F6,eps)
       endif
       intpo2 = intpo21+intpo22+intpo23+intpo24+intpo25+intpo26
      elseif(z < 100.d0) then 
       intpo22 = 0.d0
       intpo23 = 0.d0
       intpo24 = 0.d0
       intpo25 = 0.d0
       intpo26 = 0.d0
       if(Lsc > 100.d0) then
        call QSL3D(intpo22,z,100.d0,F6,eps)
        if(Lsc > 1000.d0) then
         call QSL3D(intpo23,100.d0,1000.d0,F6,eps)
         if(Lsc > 10000.d0) then
          call QSL3D(intpo24,1000.d0,10000.d0,F6,eps)
          if(Lsc > 100000.d0) then
           call QSL3D(intpo25,10000.d0,100000.d0,F6,eps)
           call QSL3D(intpo26,100000.d0,Zz(Nz),F6,eps)
          else
           call QSL3D(intpo25,10000.d0,Lsc,F6,eps)
          endif
         else
          call QSL3D(intpo24,1000.d0,Lsc,F6,eps)
         endif
        else
         call QSL3D(intpo23,100.d0,Lsc,F6,eps)
        endif
       else
        call QSL3D(intpo22,z,Lsc,F6,eps)
       endif
       intpo2 = intpo22+intpo23+intpo24+intpo25+intpo26
      elseif(z < 1000.d0) then  
       intpo23 = 0.d0
       intpo24 = 0.d0
       intpo25 = 0.d0
       intpo26 = 0.d0
       if(Lsc > 1000.d0) then
        call QSL3D(intpo23,z,1000.d0,F6,eps)
        if(Lsc > 10000.d0) then
         call QSL3D(intpo24,1000.d0,10000.d0,F6,eps)
         if(Lsc > 100000.d0) then
          call QSL3D(intpo25,10000.d0,100000.d0,F6,eps)
          call QSL3D(intpo26,100000.d0,Zz(Nz),F6,eps)
         else
          call QSL3D(intpo25,10000.d0,Lsc,F6,eps)
         endif
        else
         call QSL3D(intpo24,1000.d0,Lsc,F6,eps)
        endif
       else
        call QSL3D(intpo23,z,Lsc,F6,eps)
       endif
       intpo2 = intpo23+intpo24+intpo25+intpo26
      elseif(z < 10000.d0) then  
       intpo24 = 0.d0
       intpo25 = 0.d0
       intpo26 = 0.d0
       if(Lsc > 10000.d0) then
        call QSL3D(intpo24,z,10000.d0,F6,eps)
        if(Lsc > 100000.d0) then
         call QSL3D(intpo25,10000.d0,100000.d0,F6,eps)
         call QSL3D(intpo26,100000.d0,Zz(Nz),F6,eps)
        else
         call QSL3D(intpo25,10000.d0,Lsc,F6,eps)
        endif
       else
        call QSL3D(intpo24,z,Lsc,F6,eps)
       endif
       intpo2 = intpo24+intpo25+intpo26
      elseif(z < 100000.d0) then                      
       intpo25 = 0.d0
       intpo26 = 0.d0
       if(Lsc > 100000.d0) then
        call QSL3D(intpo25,z,100000.d0,F6,eps)
        call QSL3D(intpo26,100000.d0,Zz(Nz),F6,eps)
       else
        call QSL3D(intpo25,z,Lsc,F6,eps)
       endif
       intpo2 = intpo25+intpo26
      elseif(z > 100000.d0) then                      
       call QSL3D(intpo26,z,Zz(Nz),F6,eps)
       intpo2 = intpo26
      endif
      Vel = (1.d0/(er*e0)*intpo2)  + (1.d0/(er*e0)*Sig_gate*(z-Lsc))
     end function Vel





     real(8) function F6(zs)
      real(8)          :: zs
      F6 = -pos(zs)*(zs-zp)
     end function F6





     subroutine calc_po00                                    ! charge density of doped SC
      if(L_debug) print *,'calc_po00:'
      po00_h = poh0(0.d0)
      po00_e = poe0(0.d0)
      po00 =   po00_h + po00_e                           
      if(L_n_type) then                       
       po00_e0 = po00
       po00_h0 = 0.d0
      elseif(L_p_type) then                   
       po00_h0 = po00
       po00_e0 = 0.d0
      endif
      if(dabs(EFermi1-EFermi_00) .lt. 0.1d0) print 9
      if(L_super_debug) then
       print *,'   EFermi1     po_n'
       print 88,EFermi1,po00*1.d24
      endif
88    format(F12.3,E20.4)
 9    format(/'*** WARNING ***'/' This is too small doping concentration. Are you sure about the input parameters?'/)
     end subroutine calc_po00





     subroutine set_eps(eps,name)                                        ! set accuracy of integartion for poe and poh
      character(*)    :: name
      real(8)         :: eps
      real(8)         :: x,pp
      if(EFermi1 > EFermi_00) then
       x = (EFermi1-EFermi_00)/(ECBM-EFermi_00)
       pp = 20*(1-x)+16*x
       eps = 10**(-pp)
      elseif(EFermi1 < EFermi_00) then
       x = (EFermi_00-EFermi1)/(EFermi_00-EVBM)
       pp = 20*(1-x)+16*x
       eps = 10.D0**(-pp)
      else
       eps = 1.D-20
      endif
     end subroutine set_eps






     real(8) function poe0(z) 
      real(8)       :: eps
      real(8)       :: z
      real(8)       :: R13,R14,R15,R16,R17,R18
      call set_eps(eps,'poe0')
      call QSL3D(R13,Exxe1(1),Exxe2(1),F2,eps)        
      call QSL3D(R14,Exxe1(2),Exxe2(2),F2,eps)        
      call QSL3D(R15,Exxe1(3),Exxe2(3),F2,eps)        
      call QSL3D(R16,Exxe1(4),Exxe2(4),F2,eps)        
      call QSL3D(R17,Exxe1(5),Exxe2(5),F2,eps)        
      call QSL3D(R18,Exxe1(6),Exxe2(6),F2,eps)        
      poe0 = -(R13+R14+R15+R16+R17+R18)/V_DSC
     end function poe0




     real(8) function F2(E)
      real(8)              :: E
      F2 = DOS_SCs(E)/(1.d0+dexp((E-EFermi1)/kbT))
     end function F2




     real(8) function poh0(z) 
      real(8)       :: eps
      real(8)       :: z
      real(8)       :: R11,R12,R13,R14,R15,R16,R17
      call set_eps(eps,'poh0')
      call QSL3D(R11,Exxh1(1),Exxh2(1),F1,eps)        
      call QSL3D(R12,Exxh1(2),Exxh2(2),F1,eps)        
      call QSL3D(R13,Exxh1(3),Exxh2(3),F1,eps)        
      call QSL3D(R14,Exxh1(4),Exxh2(4),F1,eps)        
      call QSL3D(R15,Exxh1(5),Exxh2(5),F1,eps)        
      call QSL3D(R16,Exxh1(6),Exxh2(6),F1,eps)        
      call QSL3D(R17,Exxh1(7),Exxh2(7),F1,eps)        
      poh0 =  (R11+R12+R13+R14+R15+R16+R17)/V_DSC
     end function poh0




     real(8) function F1(E)
      real(8)              :: E
      real(8)              :: dexpx
      dexpx = dexp((E-EFermi1)/kbT)
      F1 = (dexpx/(1.d0+dexpx))*DOS_SCs(E)
     end function F1





     real(8) function poh(z)                                    ! charge density of holes
      real(8)        :: z
      real(8)        :: R11,R12,R13,R14,R15,R16,R17
      real(8)        :: eps
      if(L_p_type) then
       call set_eps(eps,'poh')
      else
       eps = 1.d-22
      endif
      eVz = -Vels(z)
      call QSL3D(R11,Exxh1(1)+eVz,Exxh2(1)+eVz,poh_2,eps)        
      call QSL3D(R12,Exxh1(2)+eVz,Exxh2(2)+eVz,poh_2,eps)        
      call QSL3D(R13,Exxh1(3)+eVz,Exxh2(3)+eVz,poh_2,eps)        
      call QSL3D(R14,Exxh1(4)+eVz,Exxh2(4)+eVz,poh_2,eps)        
      call QSL3D(R15,Exxh1(5)+eVz,Exxh2(5)+eVz,poh_2,eps)        
      call QSL3D(R16,Exxh1(6)+eVz,Exxh2(6)+eVz,poh_2,eps)        
      call QSL3D(R17,Exxh1(7)+eVz,Exxh2(7)+eVz,poh_2,eps)        
      poh =  (R11+R12+R13+R14+R15+R16+R17)/V_DSC                ! devide by volume of DOS for semiconductor
     end function poh



     real(8) function poh_2(E)                                   ! charge density of holes
      real(8)        :: E
      real(8)        :: dexpx
      dexpx = dexp((E-EFermi1)/kbT)
      poh_2 =  DOS_SCs(E-eVz)*(dexpx/(1.d0+dexpx)) 
     end function poh_2






     real(8) function poe(z)                                     ! charge density of electrons
      real(8)        :: z
      real(8)        :: R11,R12,R13,R14,R15,R16
      real(8)        :: eps
      if(L_n_type) then
       call set_eps(eps,'poe')
      else
       eps = 1.d-22
      endif
      eVz = -Vels(z)
      call QSL3D(R11,Exxe1(1)+eVz,Exxe2(1)+eVz,poe_2,eps)        
      call QSL3D(R12,Exxe1(2)+eVz,Exxe2(2)+eVz,poe_2,eps)        
      call QSL3D(R13,Exxe1(3)+eVz,Exxe2(3)+eVz,poe_2,eps)        
      call QSL3D(R14,Exxe1(4)+eVz,Exxe2(4)+eVz,poe_2,eps)        
      call QSL3D(R15,Exxe1(5)+eVz,Exxe2(5)+eVz,poe_2,eps)        
      call QSL3D(R16,Exxe1(6)+eVz,Exxe2(6)+eVz,poe_2,eps)        
      poe =   -(R11+R12+R13+R14+R15+R16)/V_DSC                   ! devide by volume of DOS for semiconductor      
     end function poe




     real(8) function poe_2(E)                                   ! charge density of electrons
      real(8)        :: E
      poe_2 = DOS_SCs(E-eVz)/(1.d0+dexp((E-EFermi1)/kbT))
     end function poe_2






     real(8) function poMIGS(z)                                            ! charge density of MIGS
      real(8)        :: z
      real(8)        :: R3,R31,R32,R33
      real(8)        :: eps
      zp = z
      if(z < 1000.d0) then                          
       call set_epsMIGS(z,eps)
       eVz = -Vels(z)
       call QSL3D(R31,Exxm1(1)+eVz,Exxm2(1)+eVz,poMIGS_2,eps)             
       call QSL3D(R32,Exxm1(2)+eVz,Exxm2(2)+eVz,poMIGS_2,eps)         
       call QSL3D(R33,Exxm1(3)+eVz,Exxm2(3)+eVz,poMIGS_2,eps)         
       R3 = (R31+R32+R33)
      else
       R3 = 0.d0                                                           
      endif 
      poMIGS =   -R3/V_D0           
     end function poMIGS




     real(8) function poMIGS_2(E)                                          ! charge density of MIGS
      real(8)        :: E
      poMIGS_2 =   DMIGS(E-eVz)/(1.d0+dexp((E-EFermi2)/kbT))              
     end function poMIGS_2





     subroutine set_epsMIGS(z,eps)                                         ! set accuracy of integration for poMIGS calculation
      real(8)         :: eps
      real(8)         :: z,x,pp
      if(z < 3000.d0) then
       x = z/3000.d0
       pp = 30*x+12*(1-x)
       eps = 10**(-pp)
      else
       eps = 1.d-30
      endif
     end subroutine set_epsMIGS








    real(8) function DMIGS(E)                                  ! density of MIGS in states/(eV*Vcell)
      real(8)           :: E                                   ! in eV
      integer           :: k
      real(8)           :: Dx
      real(8)           :: DMIGS1,DMIGS_G1
      real(8)           :: DLDH,cx,Imks2x
      if(zp < z3) then
       DMIGS = 0.d0
      else
!      call calc_connection(E)                                 ! calculate z0 for connection analytical and numerical results
       zconnect = 1.5d0
       if(E > 0.d0 .and. E < ECBM) then
        call separate(E,ImKHs2(E,k),ImKs2(E,k),DLDH)
       elseif(E <= 0.d0) then
        DLDH = 1.0d0 
       elseif(E >= ECBM) then
        DLDH = 0.0d0 
       endif
       if(zp < zconnect) then                                 ! int dk is calculated only for small z
        Dx = 0.d0
        do k=1,Nk                                             ! integrate over kx,ky
         if(k==1) then
          cx = DLDH
         else
          cx = 1.d0
         endif
         Dx = Dx + wk(k)*DOS_Ms(E,k)*exp(-2.d0*Imks2(E,k)*ckA*abs(zp))*cx
        enddo
        DMIGS1 = Dx/sumk
       else
        DMIGS1 = 1.d20
       endif
       if(abs(zp) > 0.d0) then
        Imks2x = Imks2(E,1)
        DMIGS_G1 = DOS_Ms(E,1)*exp(-2.d0*Imks2x*ckA*abs(zp))*(Imks2x*ckA)/abs(zp)/(a2p**2)*(1.d0-exp(-((a2p**2)/(Imks2x*ckA))*abs(zp)))
       else
        DMIGS_G1 = DOS_Ms(E,1)
       endif
       DMIGS_G1 = DMIGS_G1*DLDH
       if(zp < zconnect+0.5d0) then     
        if(DMIGS1 > DMIGS_G1) then                   ! to connect smoothly, Gaussian for small z may be larger
         DMIGS = DMIGS_G1                            ! zconnect may be calculated inaccurate, and in order to avoid bigger DMIGS_G then DMIGS 
        else
         DMIGS = DMIGS1
        endif 
       else
        DMIGS = DMIGS_G1                             ! Gaussian for large z
       endif 
      endif  
    end function DMIGS





    real(8) function DMIGS_G(E)                      ! MIGS from Gussian distribution around (kx,ky) = (0,0)
     real(8)           :: E                          ! in eV
     if(abs(zp) > 0.d0) then
      DMIGS_G = DOS_Ms(E,1)*exp(-2.d0*Imks2(E,1)*ckA*abs(zp))*(Imks2(E,1)*ckA)/abs(zp)/(a2p**2)*(1.d0-exp(-((a2p**2)/(Imks2(E,1)*ckA))*abs(zp)))
     else
      DMIGS_G = DOS_Ms(E,1)
     endif
    end function DMIGS_G




    real(8) function DiG(E,k)                        ! DOS interface with Gaussian approximation (for large distances)
     real(8)      :: E
     integer      :: k
     real(8)      :: kr2,d2
     real(8)      :: kvec(3)
     d2 = Imks2(E,1)*ckA/abs(zp)
     kvec(1:3) = kp(1,k)*b1(1:3) + kp(2,k)*b2(1:3)
     kr2 = (kvec(1)**2 + kvec(2)**2)              
     DiG = DOS_Ms(E,1)*exp(-2.d0*Imks2(E,1)*ckA*abs(zp))*dexp(-kr2/d2)
    end function DiG



    real(8) function Di(E,k)                         ! DOS interface (E,k,z)
     real(8)      :: E
     integer      :: k
     Di = DOS_Ms(E,k)*exp(-2.d0*Imks2(E,k)*ckA*abs(zp))
    end function Di

 


    real(8) function Gaus(E,k)
     real(8)    :: E
     integer    :: k
     real(8)    :: kr2,d2
     real(8)    :: kvec(3)
     d2 = Imks2(E,1)*ckA/abs(zp)
     kvec(1:3) = kp(1,k)*b1(1:3) + kp(2,k)*b2(1:3)
     kr2 = (kvec(1)**2 + kvec(2)**2)              
     Gaus = dexp(-kr2/d2)    
    end function Gaus





     subroutine calc_gaus(E)                            ! check integral
      real(8)     :: E
      integer     :: k
      real(8)     :: int
      int = 0.d0
      do k=1,Nk                                         ! integrate over kx,ky
       int = int + wk(k)*Gaus(E,k)
      enddo
      int = int/sumk*(1.d0/(pi*a2p**2))
      print *,'calc_gaus: int=',int
      print *,'calc_gaus: anal=',1.d0/(a2p**2)*(Imks2(E,1)*ckA)/abs(zp)*(1.d0 - exp(-a2p**2/(Imks2(E,1)*ckA)*abs(zp)))
      print *,'ckA=',ckA
     end subroutine calc_gaus





     subroutine calc_totq
      real(8)            :: Sig2
      real(8)            :: eps
      eps = 1.d-14
      call QSL3D(Sig2,0.d0,Zz(Nz),pos,eps)
      Sig = -Sig2
     end subroutine calc_totq






     subroutine calc_DSC_int
      real(8)            :: Ne21,Ne22,Ne23,Ne24,Ne25,Ne26,Ne27
      real(8)            :: Ne31,Ne32,Ne33,Ne34,Ne35,Ne36
      real(8)            :: eps
      real(8)            :: Exx1(5),Exx2(5)
      eps = 1.d-12
      print *
      print *,'Bulk'
      call QSL3D(Ne21,Exxh1(1),Exxh2(1),DOS_SCs,eps)
      call QSL3D(Ne22,Exxh1(2),Exxh2(2),DOS_SCs,eps)
      call QSL3D(Ne23,Exxh1(3),Exxh2(3),DOS_SCs,eps)
      call QSL3D(Ne24,Exxh1(4),Exxh2(4),DOS_SCs,eps)
      call QSL3D(Ne25,Exxh1(5),Exxh2(5),DOS_SCs,eps)
      call QSL3D(Ne26,Exxh1(6),Exxh2(6),DOS_SCs,eps)
      call QSL3D(Ne27,Exxh1(7),Exxh2(7),DOS_SCs,eps)

      call QSL3D(Ne31,Exxe1(1),Exxe2(1),DOS_SCs,eps)
      call QSL3D(Ne32,Exxe1(2),Exxe2(2),DOS_SCs,eps)
      call QSL3D(Ne33,Exxe1(3),Exxe2(3),DOS_SCs,eps)
      call QSL3D(Ne34,Exxe1(4),Exxe2(4),DOS_SCs,eps)
      call QSL3D(Ne35,Exxe1(5),Exxe2(5),DOS_SCs,eps)
      call QSL3D(Ne36,Exxe1(6),Exxe2(6),DOS_SCs,eps)

      Nee1 = Ne21+Ne22+Ne23+Ne24+Ne25+Ne26+Ne27+ &
             Ne31+Ne32+Ne33+Ne34+Ne35+Ne36
      print *,'Nee1=',Nee1
     end subroutine calc_DSC_int





     subroutine calc_DM_int                                  ! calculate integrals of DOS 
      real(8)            :: eps
      real(8)            :: NeC
      real(8)            :: Emin,Emax
      print *
      print *,'Interfacial layer'
      eps = 1.d-12
      Emin = -12.55d0
      Emax =  2.446d0 
      call QSL3D(Nee2 ,Emin, Emax, DOS_Mtots,eps)
      print *,'Int over -inf to +inf    Nee2=',Nee2
      call QSL3D(NeC,EVBM,ECBM,DOS_Mtots,eps)
      print *,'Int over EVBM to ECBM    NeC=',NeC
      call QSL3D(NeC,EVBM,ECBM,DOS0s,eps)
      print *,'Int over EVBM to ECBM for DOS0  NeC=',NeC
      DS0 = NeC/(ECBM-EVBM)/V_D0                             ! density of states 
      print *,'DS0=',DS0
     end subroutine calc_DM_int





     subroutine calc_deltaE                                  ! calculate filling level of the surface
      real(8)            :: eps
      real(8)            :: E1,E2
      dEf1 = dEf                                             ! store previous value
      eps = 0.000001d0
      if(L_n_type) then
       E1 =  0.d0
       E2 =  0.2d0
      elseif(L_p_type) then
       E1 = -0.2d0
       E2 =  0.0d0
      else
       E1 =  0.0d0
       E2 =  0.0d0
      endif
      call zero12(E1,E2,dEf,eps)
      if(dabs(dEf-0.2d0) .le. 0.001d0) then
       print *,'*** ERROR: calc_deltaE: increase searching range'
       stop
      endif
     end subroutine calc_deltaE    





     subroutine zero12(a,b,za,eps)                           ! search zero by 1/2 method
      real(8)    :: a,b,za,eps
        do while (dabs(a-b) > eps)
         za = (a+b)/2.d0
        if(L_n_type) then
         if(Fsigma(za) > 0.d0) then
          b = za
         else
          a = za
         endif
        elseif(L_p_type) then
         if(Fsigma(za) > 0.d0) then
          a = za
         else
          b = za
         endif
        endif
        enddo
     end subroutine zero12




     real(8) function Fsigma(dE)
      real(8)            :: dE
      real(8)            :: Nsigma,SInt
      real(8)            :: eps
      real(8)            :: E1,E2
      real(8)            :: SigS                                            ! initial surface charge on the interface
      eps = 1.d-12
      EFermi111 = CNL + dE
      if(dE > 0.d0) then
       E1 = CNL
       E2 = CNL+dE
      elseif(dE < 0.d0) then
       E1 = CNL+dE
       E2 = CNL
      else
       E1 = CNL
       E2 = CNL
      endif
      call QSL3D(SInt,E1,E2,F99,eps)
      SigS = 0.d-3                      ! 1E-3 in e/A^2     = 1E13 cm-2   (x1E16)
      Nsigma = (dabs(Sig)+dabs(SigS))*Surface
      Fsigma = Sint - Nsigma 
     end function Fsigma




     real(8) function F99(E)                                  
      real(8)        :: E
      real(8)        :: dexpx
      if(L_n_type) then
       F99 = DOS0s(E)/( 1.d0+dexp((E-EFermi111)/kbT) )
      elseif(L_p_type) then
       dexpx = dexp((E-EFermi111)/kbT)
       F99 = DOS0s(E)*( dexpx/( 1.d0+dexpx ) )
      endif
     end function F99






       subroutine calc_zero_EF                                             ! calculate EFermi for intrinsic SC
        integer :: i
        if(L_debug) then
         print *,'calc_zero_EF:'
         print *,'looking for EF between EVBM=',EVBM,' and ECBM=',ECBM
        endif
        call calc_EFermi(1.d-15)
        if(L_debug) then
         print *
         print *,'EFermi1=',EFermi1
         print *,'po00=',po00                                               ! this should be close to zero
         print *,'po00_h=',po00_h                                           ! h and e compensate each other
         print *,'po00_e=',po00_e
        endif
        EFermi_00 = EFermi1
        print *,'Temperature =',Temp
        print *,'Fermi level for intrinsic SC =',EFermi1
        print *
        if(EFermi_input > EFermi_00) then                                  ! n-type
         L_n_type = .true.
         L_p_type = .false.
        elseif(EFermi_input < EFermi_00) then                              ! p-type
         L_n_type = .false.
         L_p_type = .true.
        else                                                               ! intrinsic
         L_n_type = .false.
         L_p_type = .false.
        endif
        EFermi1 = EFermi_input                                             ! set Fermi level for the system
       end subroutine calc_zero_EF





       subroutine calc_EFermi(eps)
        real(8)    :: eps
        real(8)    :: a,b
        if(L_debug) print *,'calc_EFermi'
        a = 0.20d0                    ! EVBM
        b = 1.40d0     ! ECBM
        do while (dabs(a-b) > eps)
         EFermi1 = (a+b)/2.d0
         call calc_po00
         if(dabs(po00_h) > dabs(po00_e)) then
          a = (a+b)/2.d0
         else
          b = (a+b)/2.d0
         endif
        enddo
       end subroutine calc_EFermi









      subroutine calc_connection(Ex)                               ! calculate connection between numerical and analytical solutions
       real(8)      :: Ex
       real(8)      :: d2                                          ! numerical width
       real(8)      :: G1,G2
       real(8)      :: kvec(3)
       real(8)      :: kr2
       integer      :: k9
       real(8)      :: logG2G1
       real(8)      :: Kz0,Kz9
       k9 = 9                                                      ! k=9 and 16 point in k-mesh for 15x15 
       G1 = DOS_Ms(Ex,1)                                           ! value of Di(E,k) at Gamma
       G2 = DOS_Ms(Ex,k9)                                          ! value of Di(E,k) at nearest point at Gamma-K pathway (no. 9 in 15x15 k-mesh) 
       kvec(1:3) = kp(1,k9)*b1(1:3) + kp(2,k9)*b2(1:3)             ! 1/A
       kr2 = (kvec(1)**2 + kvec(2)**2)                             ! kr**2     kvec(3) = 0
       logG2G1 = dlog(G2/G1)
       d2 = -kr2/dlog(G2/G1)
       Kz9 = ImKs2(Ex,k9)*ckA
       Kz0 = ImKs2(Ex,1)*ckA
       zconnect = logG2G1/(2.d0*(Kz9-Kz0)-kr2/Kz0)
       if(zconnect < 0.d0 .or. zconnect > 10.d0) then
        print *
        print *,'calc_connection:'
        print *,'Ex=',Ex
        print *,'kp=',kp(1:3,k9)
        print *,'kvec=',kvec(1:3)
        print *,'kr2=',kr2,' 1/A^2'
        print *,'logG2G1=',logG2G1
        print *,'d2 numeric=',d2,'  1(A^2)'
        print *,'Kz0=',Kz0,' 1/A'
        print *,'Kz9=',Kz9,' 1/A'
        print *,'**** WARNING:'
        print *,'zconnect=',zconnect
       endif
      end subroutine calc_connection




       subroutine calc_DOS_int
        call calc_DM_int
        call calc_DSC_int
       end subroutine calc_DOS_int




     end Module SBCalc_i











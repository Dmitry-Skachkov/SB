



     Module SBPlot                                           ! plot results
      use SBParameters
      use SBSpline_functions
      implicit none
     contains




     subroutine write_results
      integer       :: i
      open(unit=2,file='elpot_.dat' )
       do i=1,Nz
        write(2,2) Zz(i),-V_eln(i),El_f2(i)                 
       enddo
      close(unit=2) 
      open(unit=2,file='po_.dat' )
       write(2,3) Nz
       do i=1,Nz
        write(2,1) Zz(i),po_h(i)*1.d24,po_e(i)*1.d24,po_MIGS(i)*1.d24,po_new(i)*1.d24   
       enddo
      close(unit=2)
 1    format(F17.5,4E21.12e3)
 2    format(F17.5,F18.10,E17.5e3)
 3    format('VARIABLES = "z", "poh", "poe", "poMIGS", "po"'/'ZONE  I=',I4,' F=POINT')
     end subroutine write_results




      subroutine print_results
        E_00 = -(V_eln(4)-V_eln(3))/(Zz(4)-Zz(3))
        print 5,E_00
        P_00 = kappa*E_00
        print 6,P_00
        print 1
        print 2,filen,dabs(-V_eln(1))-dabs(dEf),Sig,bspl4(1)*er*e0,poh_max*1.d24,poe_max*1.d24,poMe_max*1.d24,poMh_max*1.d24,po_new(1)*1.d24,delta_po*1.d24,delta_V,delta_Vm
        if(L_n_type) then
         SBH = dabs(-V_eln(1))-dEf + (ECBM-EFermi1) 
        elseif(L_p_type) then
         SBH = dabs(-V_eln(1))-dabs(dEf) + EFermi1
        endif
        print *,'SBH=',SBH 
        print 7,Sig*1.d16 
        print 8,dEf
        if(L_scf) print 3,filen,V_eln(1)
        if(L_conv) print 4
1       format('  it    -eV(0)      Sig(po)     dV/dz(er*e0)         poh_m          poe_m        poMe_m        poMh_m         po(0)        delta_po2     delta_V     delta_Vm')
2       format(I4,F9.5,2E15.5,7E14.4,F12.4)
3       format(I4,F9.5,'     divergency due to V'/'STOP')
4       format('We have reached convergency with respect to po')
5       format(' E(0)=',E15.5,'  V/A')
6       format(' P(0)=',E15.5,'  e/A^2')
7       format('  Sig=',1p,E15.5,' cm-2')
8       format(' delta E = ',1p,E15.4,' eV')
      end subroutine print_results




   end module SBPlot





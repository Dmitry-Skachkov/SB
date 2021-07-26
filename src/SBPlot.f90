



     Module SBPlot                                           ! plot and print results
      use SBParameters
      use SBSpline_functions
      implicit none
     contains



     subroutine print_logo
      print 1
  1   format(///                                                                                           &
             '   ************************************************************************************' /   &
             '   *                                                                                  *' /   &
             '   *               Program SB: First-Principles Method for Schottky Barrier           *' /   &
             '   *               This program is the open-source software, please cite              *' /   &
             '   *               Phys. Rev. B, 2021                                                 * '/   &
             '   *                                                                                  *' /   &
             '   *               Center for Molecular Magnetic Quantum Materials (M2QM)             *' /   &
             '   *               Quantum Theory Project (QTP)                                       *' /   &
             '   *               Department of Physics                                              *' /   &
             '   *               University of Florida, Gainesville, FL 32611, USA                  *' /   &
             '   *               https://efrc.ufl.edu/                                              *' /   &
             '   *               For contributing please contact <m2qm.efrc@phys.ufl.edu>           *' /   &
             '   *                                                                                  *' /   &
             '   ************************************************************************************' /// )
     end subroutine print_logo




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
      call write_restart_dat
 1    format(F17.5,4E21.12e3)
 2    format(F17.5,F18.10,E17.5e3)
 3    format(I5,'     !  "z", "poh", "poe", "poMIGS", "po"')
     end subroutine write_results




      subroutine print_results
        if(L_scf) print 3,iter,V_eln(1)
        if(L_conv) print 4
        E_00 = -(V_eln(4)-V_eln(3))/(Zz(4)-Zz(3))
        P_00 = kappa*E_00
        if(L_p_type) print 12
        if(L_n_type) print 13
        print 14,Temp
        print 10,po00*1.d24
        print 11,EFermi1
        print 21,CNL
        print 16,EFermi_00
        print 9,SBH 
        print 18,-V_eln(1)
        print 8,dEf
        print 15,Sig_gate*1.d16
        print 7,Sig*1.d16 
        print 5,E_00
        print 6,P_00
        print 19,DLW
        print 20,ILW
        print 17
3       format(I4,F11.5,'     divergency due to V'/'STOP')
4       format(' Convergence with respect to charge density is reached')
5       format(' Electric field close to the surface of semiconductor E(0) =',1p,E14.5,' V/A')
6       format(' Polarization close to the surface of semiconductor   P(0) =',1p,E14.5,' e/A^2')
7       format(' Charge on the interface                               Sig =',1p,E14.5,' cm-2')
8       format(' Energy filling level on the surface               delta E =',1p,E14.5,' eV')
9       format(' Schottky barrier height                               SBH =',   F14.5,' eV')
10      format(' Doping concentration                                 po00 =',1p,E14.4,' cm-3')
11      format(' Fermi level                                        EFermi =',   F14.5,' eV')
12      format(///' p-type doping')
13      format(///' n-type doping')
14      format(' Temperature                                               =',   F14.3,' K')
15      format(' Surface charge density of the gate electrode     Sig_gate =',1p,E14.4,' cm-2')
16      format(' Fermi level for intrinsic semiconductor         EFermi_00 =',   F14.5,' eV')
17      format(///' NORMAL TERMINATION'/)
18      format(' Amplitude of Schottky potential energy               -eV0 =',   F14.5,' eV')
19      format(' Depletion layer width                                 DLW =',   F14.2,' A')
20      format(' Inversion layer width                                 ILW =',   F14.2,' A')
21      format(' Charge neutrality level for interface vs VBM          CNL =',   F14.5,' eV')
      end subroutine print_results




   end module SBPlot





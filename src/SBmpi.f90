
 Module SBmpi                                                                  ! MPI parallel calculation
      use mpi
      !   include 'mpif.h'
      integer, parameter   :: sumtype = 17
      integer, parameter   :: sizetype = 10
      integer, parameter   :: masternode = 0
      integer              :: Process             ! 0 .. Numnodes-1
      integer              :: NumNodes
      integer              :: MPIerror
      integer              :: MPIstatus(3)
      real(8)              :: Pt0,Ptf
      integer              :: P_MPI_group,P_MPI_start
      real(8), allocatable :: A1(:),A2(:),A3(:),A4(:) ! working arrays
      real(8)              :: R1,R2,R3,R4
      integer              :: iPrint
 contains


   subroutine P_start                                                          ! initialization of parallel calculation
        call MPI_INIT( MPIerror )                                        
          if(MPIerror/=0) print *,' P_start: error in MPI_INIT'
        call MPI_COMM_RANK( MPI_COMM_WORLD, Process, MPIerror )          
          if(MPIerror/=0) print *,' P_start: error in MPI_COMM_RANK'
        call MPI_COMM_SIZE( MPI_COMM_WORLD, Numnodes, MPIerror )         
          if(MPIerror/=0) print *,' P_start: error in MPI_COMM_SIZE'
        call cpu_time(Pt0)
        call MSG('Start parallel calculation on '//Pfstr(Numnodes)//' cores    ')
        iPrint = 1000 + Process
        open(unit=iPrint,file='Process'//Pfstr(iPrint)//'.out')
   end subroutine P_start


   subroutine P_sendI(Int)                                                     ! send data to all parallel processors
        integer     :: Int
        call MPI_BCAST(Int,1,MPI_INTEGER,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendI: error in BCAST'
   end subroutine P_sendI


   subroutine P_sendL(Int)                                                     ! send data to all parallel processors
        logical     :: Int
        call MPI_BCAST(Int,1,MPI_LOGICAL,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendL: error in BCAST'
   end subroutine P_sendL


   subroutine P_sendI1(Int,N)                                                     ! send data to all parallel processors
        integer     :: Int(N)
        integer     :: N
        call MPI_BCAST(Int,N,MPI_INTEGER,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendI: error in BCAST'
   end subroutine P_sendI1


   subroutine P_sendI2(Int,N1,N2)                                                     ! send data to all parallel processors
        integer     :: Int(N1,N2)
        integer     :: N1,N2
        call MPI_BCAST(Int(1,1),N1*N2,MPI_INTEGER,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendI: error in BCAST'
   end subroutine P_sendI2


   subroutine P_sendR(R)
        real*8      :: R
        call MPI_BCAST(R,1,MPI_DOUBLE_PRECISION,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendR: error in BCAST'
   end subroutine P_sendR


   subroutine P_sendR1(R,N)
        integer     :: N
        real*8      :: R(N)
        call MPI_BCAST(R,N,MPI_DOUBLE_PRECISION,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendR: error in BCAST'
   end subroutine P_sendR1


   subroutine P_sendR2(R,N1,N2)
        integer     :: N1,N2
        real*8      :: R(N1,N2)
        call MPI_BCAST(R(1,1),N1*N2,MPI_DOUBLE_PRECISION,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendR: error in BCAST'
   end subroutine P_sendR2


   subroutine P_sendR3(R,N1,N2,N3)
        integer     :: N1,N2,N3
        real*8      :: R(N1,N2,N3)
        call MPI_BCAST(R(1,1,1),N1*N2*N3,MPI_DOUBLE_PRECISION,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendR: error in BCAST'
   end subroutine P_sendR3


   subroutine P_sendA(A)
        character(*)     :: A
        integer          :: l
        l = len(A)
        call MPI_BCAST(A,l,MPI_CHARACTER,Masternode,MPI_COMM_WORLD,MPIerror)             
          if(MPIerror/=0) print *,' P_sendA: error in BCAST'
   end subroutine P_sendA


   subroutine P_combine_R(A,R,Nn)                                         ! collect array R from parallel processors
        integer    :: Nn
        real(8)    :: A(:)
        real*8     :: R(Nn)
        call MPI_Barrier(MPI_COMM_WORLD, MPIerror)                             ! wait all processors
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_BARRIER'
        call MPI_AllGather(A,P_MPI_group,MPI_DOUBLE_PRECISION,R,P_MPI_group,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,MPIerror)   ! Masternode receives data from others cores 
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_GATHER'
   end subroutine P_combine_R


   subroutine P_sum_R(Sumi,Sum)                                         ! collect array R from parallel processors
        real(8)    :: Sumi
        real*8     :: Sum
        call MPI_Barrier(MPI_COMM_WORLD, MPIerror)                             ! wait all processors
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_BARRIER'
        call MPI_AllReduce(Sumi,Sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerror)   ! Masternode receives data from others cores 
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_GATHER'
   end subroutine P_sum_R
  

   subroutine P_max_R(Maxi,Max)                                         ! collect array R from parallel processors
        real(8)    :: Maxi
        real*8     :: Max
        call MPI_Barrier(MPI_COMM_WORLD, MPIerror)                             ! wait all processors
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_BARRIER'
        call MPI_AllReduce(Maxi,Max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,MPIerror)   ! Masternode receives data from others cores 
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_GATHER'
   end subroutine P_max_R


   subroutine P_min_R(Mini,Min)                                         ! collect array R from parallel processors
        real(8)    :: Mini
        real*8     :: Min
        call MPI_Barrier(MPI_COMM_WORLD, MPIerror)                             ! wait all processors
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_BARRIER'
        call MPI_AllReduce(Mini,Min,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,MPIerror)   ! Masternode receives data from others cores 
          if(MPIerror/=0) print *,' P_receiveR: error in MPI_GATHER'
   end subroutine P_min_R
   

   subroutine P_wait                                         ! wait all processes
        call MPI_Barrier(MPI_COMM_WORLD, MPIerror)                             ! wait all processors
          if(MPIerror/=0) print *,' P_wait: error in MPI_BARRIER'
   end subroutine P_wait


   subroutine P_calc_group(Nn)              ! calculate groups for calculation
    integer :: ist,ngr,Nn                         ! parallel calculation from ist+1 to ist+ngr
    integer :: ifin
    ngr = Nn/NumNodes                             ! number of group of points for calculation
    if(mod(Nn,NumNodes) /= 0) ngr = ngr + 1
    ist = Process*ngr                             ! starting point
    ifin = ist + ngr                              ! finil point
    if(ifin > Nn) ngr = Nn - ist                  ! correct number of group for last process
    P_MPI_group = ngr
    P_MPI_start = ist
    allocate(A1(P_MPI_group))
    allocate(A2(P_MPI_group))
    allocate(A3(P_MPI_group))
    allocate(A4(P_MPI_group))
   ! print *,'Process=',Process,' ist=',ist,' ngr=',ngr
   end subroutine P_calc_group


   subroutine MSG(comment)
    character(*)          :: comment
    integer               :: l1
    l1 = len(trim(adjustl(comment)))
    if(Process==0) print 1,trim(adjustl(comment))
 1  format(A<l1>)
   end subroutine MSG


   subroutine MSGT(comment)
    character(*)          :: comment
    integer               :: l1,l2
    call cpu_time(Ptf)
    l1 = len(trim(adjustl(comment)))
    if(l1>=60) then
     if(Process==0) print 1,trim(adjustl(comment)),(Ptf-Pt0)
    else
     l2 = 60 - l1
     if(Process==0) print 2,trim(adjustl(comment)),' ',(Ptf-Pt0)
    endif
    Pt0 = Ptf
 1  format(A<l1>,' (',F10.3,' s)')
 2  format(A<l1>,A<l2>,' (',F10.3,' s)')
   end subroutine MSGT


   subroutine P_stop                                                ! finalizing parallel calculation
     deallocate(A1)
     deallocate(A2)
     deallocate(A3)
     deallocate(A4)   
     close(unit=iPrint)  
     call MSGT('Finish parallel calculation on '//Pfstr(Numnodes)//' cores    ')
     call MPI_FINALIZE(MPIerror)
     if(MPIerror/=0) print *,' P_stop: error in MPI_FINALIZE'
   end subroutine P_stop


     character(len=4) function Pfstr(k)                             ! Convert an integer to character(5)
      integer, intent(in) :: k
      write(Pfstr,'(I4)') k
     end function Pfstr


 end module SBmpi



program mainMC
  use amp4hef
  use main_DrellYan
  implicit none
  
  integer,parameter :: eventUnit=21
  real(fltknd),parameter :: alphaWeak = 0.042227735_fltknd
  real(fltknd),parameter :: pi=3.1415926535897932384626433832795_fltknd
  
  logical :: exitLoop
  integer :: Noffshell,Nfinst,Ntotal,ii,process(13)
  integer :: id1,Nperm,NhelConf, NZ
  real(fltknd) :: Ecm,kTsq(2),sHat
  real(fltknd) :: momenta(0:3,8),directions(0:3,5),ampSquared
  real(fltknd) :: eventWeight,matElemData,psWeight,cnstWeight,totalWeight
  real(fltknd) :: partonLumi,alphaStrong,flux, renomScale
  real(fltknd) :: sumW0,sumW1,sumW2
  character(72) :: eventFile,line
  complex(fltknd) :: amp(64,120),ampStr(64,120),factor
  
  real(fltknd) :: S, P1(0:3), P2(0:3), xF, qT, M, x1, x2, k1T2, k2T2, z, fi_kappa,fi_k1, &
                  fi_k2, kT(1:2), xq, fi_k ,kT2
  real(fltknd) :: k1(0:3), k2(0:3), p3(0:3), p4(0:3), q(0:3)
! Obvious.
!  call get_command_argument(1,eventFile)
!!
!! Open event file and cycle through first few irrelevant lines.
!  open(eventUnit,file=trim(eventFile),status='old')
!  do
!    read(eventUnit,'(A)') line
!    if (trim(line).eq.'BEGINOFFILE') exit
!    write(*,*) line
!  enddo
!
!! Read the center-of-mass energy,
!! the number of off-shell partons, and the number of final-state partons.
!  read(eventUnit,*) Ecm ,Noffshell ,Nfinst, NZ
!  Ntotal = Nfinst+2
!
!! Read the process. 0=gluon ,positive=quark ,negative=anti-quark.
!! All flavors are outgoing, ie sum(process(1:Nfinst+2)) must be 0.
!  read(eventUnit,*) process(1:Ntotal)
!
!! Put the processes, and get id.
!  call put_process( id1 ,Ntotal ,Noffshell ,NZ ,process )
!! Calculate the overall constant to the event weights.
!  cnstWeight = 1 &
!!   Conversion factor from GeV to nanobarn.
!    * 389379.66_fltknd &
!!   The phase space weight "psWeight" is missing the following factor,
!!   according to the RAMBO convention.
!    * (2*pi)**(4-3*(Nfinst)) &
!!   Conversion factor from alphaStrong to g_QCD^2
!    * (4*pi)**(Nfinst) &
!!   Average over initial-state spins. The weight factor "partonLumi" includes a
!!   factor 2 for each off-shell gluon to have a uniform definition of cnstWeight.
!    / (2*2)
!!   Average over initial-state colors and symmetry factor are included in ampSquared.
!
!! Initialize statistics.
!  sumW0 = 0
!  sumW1 = 0
!  sumW2 = 0
!
!! Start the loop over events.
!  do
!!   Besides the external momenta, read_event reads the following from the event file:
!!    instWeight: the weight from the generation of initial-state variables
!!      psWeight: the weight from the generation of final-state phase space
!!    partonLumi: product of the pdfs for the initial-state partons
!!   alphaStrong: the value of the strong coupling
!!   eventWeight: the correct value of the weight determined during the
!!                creation of the event file
!    call read_event
!    if (exitLoop) exit
!
!!   The event file includes events that fell outside the cuts. These count
!!   as events where the integrand vanishes.
!    sumW0 = sumW0+1
!    if (eventWeight.eq.0) cycle
!
!!   Determine the flux factor.
!    sHat = momSquared( momenta(:,1)+momenta(:,2) )
!    kTsq(1:2) = momenta(1,1:2)**2 + momenta(2,1:2)**2
!    flux = 2*sqrt(( sHat + kTsq(1) + kTsq(2) )**2 - 4*kTsq(1)*kTsq(2))
!!   Construct directions from the energy and the z-component of the initial-state momenta.
!    directions(0,1:2) = momenta(0,1:2)
!    directions(1,1:2) = 0
!    directions(2,1:2) = 0
!    directions(3,1:2) = momenta(3,1:2)
!
!    directions(0,Ntotal) = 1
!    directions(1,Ntotal) = 0
!    directions(2,Ntotal) = 0
!    directions(3,Ntotal) = 1
!
!!   Evaluate the matrix element. Call matrix_element_b which includes factors to
!!   average over initial-state colors, and the final-state symmetry factor.
!    call put_momenta( id1 ,momenta ,directions )
!    call matrix_element_b( id1 ,ampSquared )
!    write(*,*) "matrix element ratio :", matElemData/ampSquared
!!   Determine the total weight of the event.
!    !alphaStrong = 1.
!    totalWeight = cnstWeight / flux * partonLumi &
!                * ampSquared
!
!
!!   Gather statistics.
!    sumW1 = sumW1 + totalWeight
!    sumW2 = sumW2 + totalWeight**2
!!    write(*,*) '( re-calculated weight )/( weight from file ):' &
!!              ,totalWeight/eventWeight
!!   Compare calculated weight with the number from the file.
!!    write(*,*) '( re-calculated weight )/( weight from file ):' &
!!              ,totalWeight/eventWeight
!!
!    !call matrix_element_a( id1 ,ampSquared )
!    !write(*,*) ampSquared
!  enddo
!
!! Obvious.
!  close(eventUnit)



S = 64000000.0
x1 = 0.14012238979339611
x2 = 0.26041128635406519
k1T2 =  387690.91982745461
k2T2 = 1908063.0748443501
fi_k1 = 4.6984273583699538
fi_k2 = 1.2972520490486059
qT = 461.50702546109113
xF = 2.1887153949959298E-002

M = 91.187600000000003
z = 7.9904496669769356E-002

fi_kappa = 4.6200184696872997



  call matrix_element_2x3(S, x1, x2, k1T2, k2T2, fi_k1, fi_k2,  qT, xF ,M ,z &
                               ,fi_kappa, ampSquared)
  write(*,*) "2-jet process"
  write(*,*) ampSquared/22416893515.195316
S = 64000000.0
x1 = 0.19533669948577898
x2 = 0.22782330513000509
k1T2 =  515690.65307078656
k2T2 = 1988096.8940659165
fi_k1 = 4.2053055652599260
fi_k2 = 1.5631554846008653
qT = 461.50702546109113
xF = 6.3492177845316133E-002

M = 91.187600000000003
z = 0.15797272920608535

fi_kappa = 1.5532765884057573

  call matrix_element_2x3(S, x1, x2, k1T2, k2T2, fi_k1, fi_k2,  qT, xF ,M ,z &
                               ,fi_kappa, ampSquared)
  write(*,*) ampSquared/35547311794.732552



S = 64000000.0
x1 = 9.4693064689636314E-002
x2 = 0.38839664459228551
k1T2 = 391367.86529486376
k2T2 = 2457332.9217060152
fi_k1 = 3.8171116533087472
fi_k2 = 0.85936659714704566
qT = 461.50702546109113
xF = 4.6158450045114117E-002
M = 91.187600000000003
z = 0.90467376708984437

fi_kappa = 5.9103488798390043

  call matrix_element_2x3(S, x1, x2, k1T2, k2T2, fi_k1, fi_k2,  qT, xF ,M ,z &
                               ,fi_kappa, ampSquared)
    write(*,*) ampSquared/86332404941.668747

write(*,*) "1-jet process"

kT2 =0.00
!kT2 = 391367.86529486376
fi_k =  3.8171116533087472
xq = 0.58839664459228551
  call matrix_element_2x2(S,  xF, qT, M, xq, kT2, fi_k, ampSquared)
  write(*,*) ampSquared
  contains


    subroutine read_event
    read(eventUnit,'(A)') line
    exitLoop = (trim(line).eq.'ENDOFFILE')
    if (exitLoop) return
    read(eventUnit,*) eventWeight
    if (eventWeight.eq.0) return
    read(eventUnit,*) momenta(0:3,1)
    read(eventUnit,*) momenta(0:3,2)
    do ii=1,Nfinst
      read(eventUnit,*) momenta(0:3,ii+2)
    enddo
    read(eventUnit,*) matElemData, partonLumi,alphaStrong, renomScale
    end subroutine
    function momSquared( qq ) result(rslt)
    intent(in) :: qq
    real(fltknd) :: qq(0:3),rslt
    rslt = ( qq(0)-qq(3) )*( qq(0)+qq(3) )- qq(1)*qq(1) - qq(2)*qq(2)
    end function

  
end program

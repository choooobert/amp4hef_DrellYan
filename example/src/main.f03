program mainMC
  use amp4hef
  implicit none
  
  integer,parameter :: eventUnit=21
  real(fltknd),parameter :: alphaWeak = 1./29.
  real(fltknd),parameter :: pi=3.1415926535897932384626433832795_fltknd
  
  logical :: exitLoop
  integer :: Noffshell,Nfinst,Ntotal,ii,process(13)
  integer :: id1,Nperm,NhelConf, NZ
  real(fltknd) :: Ecm,kTsq(2),sHat
  real(fltknd) :: momenta(0:3,8),directions(0:3,3),ampSquared
  real(fltknd) :: eventWeight,instWeight,psWeight,cnstWeight,totalWeight
  real(fltknd) :: partonLumi,alphaStrong,flux
  real(fltknd) :: sumW0,sumW1,sumW2
  character(72) :: eventFile,line
  complex(fltknd) :: amp(64,120),ampStr(64,120),factor
  
! Obvious.
  call get_command_argument(1,eventFile)
!
! Open event file and cycle through first few irrelevant lines.
  open(eventUnit,file=trim(eventFile),status='old')
  do
    read(eventUnit,'(A)') line
    if (trim(line).eq.'BEGINOFFILE') exit
    write(*,*) line
  enddo
  
! Read the center-of-mass energy,
! the number of off-shell partons, and the number of final-state partons.
  read(eventUnit,*) Ecm ,Noffshell ,Nfinst, NZ
  Ntotal = Nfinst+2

! Read the process. 0=gluon ,positive=quark ,negative=anti-quark.
! All flavors are outgoing, ie sum(process(1:Nfinst+2)) must be 0.
  read(eventUnit,*) process(1:Ntotal)

! Put the processes, and get id.
  call put_process( id1 ,Ntotal ,Noffshell ,NZ ,process )
! Calculate the overall constant to the event weights.
  cnstWeight = 1 &
!   Conversion factor from GeV to nanobarn.
    * 389379.66_fltknd &
!   The phase space weight "psWeight" is missing the following factor,
!   according to the RAMBO convention.
    * (2*pi)**(4-3*(Nfinst-1)) &
!   Conversion factor from alphaStrong to g_QCD^2
    * (4*pi)**(Nfinst-1) &
!   Average over initial-state spins. The weight factor "partonLumi" includes a
!   factor 2 for each off-shell gluon to have a uniform definition of cnstWeight.
    / (2*2)
!   Average over initial-state colors and symmetry factor are included in ampSquared.
  
! Initialize statistics.
  sumW0 = 0
  sumW1 = 0
  sumW2 = 0
  
! Start the loop over events.
  do
!   Besides the external momenta, read_event reads the following from the event file:
!    instWeight: the weight from the generation of initial-state variables
!      psWeight: the weight from the generation of final-state phase space
!    partonLumi: product of the pdfs for the initial-state partons
!   alphaStrong: the value of the strong coupling
!   eventWeight: the correct value of the weight determined during the
!                creation of the event file
    call read_event
    if (exitLoop) exit

!   The event file includes events that fell outside the cuts. These count
!   as events where the integrand vanishes.
    sumW0 = sumW0+1
    if (eventWeight.eq.0) cycle
  
!   Determine the flux factor.
    sHat = momSquared( momenta(:,1)+momenta(:,2) )
    kTsq(1:2) = momenta(1,1:2)**2 + momenta(2,1:2)**2
    flux = 2*sqrt(( sHat + kTsq(1) + kTsq(2) )**2 - 4*kTsq(1)*kTsq(2))
!   Construct directions from the energy and the z-component of the initial-state momenta.
    directions(0,1:2) = momenta(0,1:2)
    directions(1,1:2) = 0
    directions(2,1:2) = 0
    directions(3,1:2) = momenta(3,1:2)

!   Evaluate the matrix element. Call matrix_element_b which includes factors to
!   average over initial-state colors, and the final-state symmetry factor.
    call put_momenta( id1 ,momenta ,directions )
    call matrix_element_b( id1 ,ampSquared )
    write(*,*) "matrix element squared from data: ", instWeight
    write(*,*) "matrix element calculated by me : ", ampSquared
    write(*,*) "ratio :", instWeight/ampSquared
!   Determine the total weight of the event.
    totalWeight = cnstWeight / flux * instWeight * psWeight * partonLumi &
                * ampSquared * alphaStrong**(Nfinst-1) *alphaWeak

!   Gather statistics.
    sumW1 = sumW1 + totalWeight
    sumW2 = sumW2 + totalWeight**2

!   Compare calculated weight with the number from the file.
    write(*,*) '( re-calculated weight )/( weight from file ):' &
              ,totalWeight/eventWeight
!
    !call matrix_element_a( id1 ,ampSquared )
    !write(*,*) ampSquared 
  enddo

! Obvious.  
  close(eventUnit)

! Finalize statistics.
  write(*,*) 'cross section (in nb):',sumW1/sumW0
  write(*,*) 'error estimate (in %):',100*sqrt((sumW2*sumW0/sumW1**2-1)/(sumW0-1))
 
 
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
    read(eventUnit,*) instWeight, partonLumi,alphaStrong, psWeight
    end subroutine
    function momSquared( qq ) result(rslt)
    intent(in) :: qq
    real(fltknd) :: qq(0:3),rslt
    rslt = ( qq(0)-qq(3) )*( qq(0)+qq(3) )- qq(1)*qq(1) - qq(2)*qq(2)
    end function
  
  
end program

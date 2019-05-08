! Main module providing the routines decribed in the README file


module amp4hef

  use amp4hef_io ,only: set_path,split_line
  use amp4hef_qomentum ,only: fltknd ,NsizeProc,NsizeFlavor &
                             ,init_qomentum ,qomentum_list_type
  use amp4hef_aux ,only: NcolDof ,factorial
  use amp4hef_DrellYan

  implicit none
  private
  public :: set_path ,split_line ,fltknd ,init_amp4hef
  public :: put_process ,put_momenta
  public :: amplitude,amplitude_cpp
  public :: matrix_element_a ,matrix_element_b
  public :: all_amplitudes,all_amplitudes_cpp

  logical,private,save :: initz=.true.

  integer,save :: Nproc=0
  type(qomentum_list_type),allocatable,save :: glob(:)

contains

  subroutine increase_glob( id )
  integer,intent(out) :: id
  type(qomentum_list_type),allocatable :: tmpGlob(:)
  integer :: nn
  Nproc = Nproc+1
  if (.not.allocated(glob)) allocate(glob(2))
  nn = size(glob,1)
  if (Nproc.gt.nn) then
    allocate(tmpGlob(1:nn))
    tmpGlob(1:nn) = glob(1:nn)
    deallocate(glob)
    allocate(glob(1:Nproc*2))
    glob(1:nn) = tmpGlob(1:nn)
    deallocate(tmpGlob)
  endif
  id = Nproc
  end subroutine



  subroutine init_amp4hef
  if (initz) then
    initz=.false.
  else
    return
  endif
  write(*,'(A)') '########################################################################'
  write(*,'(A)') '#                                                                      #'
  write(*,'(A)') '#                        You are using AMP4HEF                         #'
  write(*,'(A)') '#                                                                      #'
  write(*,'(A)') '# for the evaluation of helicity amplitudes with off-shell partons.    #'
  write(*,'(A)') '#                                                                      #'
  write(*,'(A)') '# author: Andreas van Hameren <hamerenREMOVETHIS@ifj.edu.pl>           #'
  write(*,'(A)') '#   date: 07-01-2019                                                   #'
  write(*,'(A)') '#                                                                      #'
  write(*,'(A)') '# contributers: Marcin Bury <marcin.buryREMOVETHIS@ifj.edu.pl>         #'
  write(*,'(A)') '#               Kacper Bilko <kacperbilkoREMOVETHIS@gmail.com>         #'
  write(*,'(A)') '#               Hubert Milczarek <hb.milczarekREMOVETHIS@gmail.com>    #'
  write(*,'(A)') '#               Mirko Serino <mirkos.serinoREMOVETHIS@gmail.com>       #'
  write(*,'(A)') '#                                                                      #'
  write(*,'(A)') '# Please cite                                                          #'
  write(*,'(A)') '#   M. Bury and A. van Hameren, Comput.Phys.Commun. 196 (2015) 592-598 #'
  write(*,'(A)') '# in publications with results obtained with the help of this program. #'
  write(*,'(A)') '#                                                                      #'
  write(*,'(A)') '########################################################################'
  call init_qomentum
  call fill_matrices_DrellYan
  end subroutine


  subroutine put_process( id ,Ntotal ,Noffshell ,NZbos ,process )
  integer,intent(out) :: id
  integer,intent(in) :: Ntotal ,Noffshell ,NZbos ,process(*)
  integer :: ii,jj,kk,Noff2,Nsum
  integer :: NflavorFinst(-NsizeFlavor:NsizeFlavor)
!
  if (initz) call init_amp4hef
  call increase_glob( id )
!
  associate( Ntot=>glob(id)%Ntot ,Noff=>glob(id)%Noff ,NZ=>glob(id)%NZ &
            ,Nflavor=>glob(id)%Nflavor ,flavor=>glob(id)%flavor &
            ,NhelOrder=>glob(id)%NhelOrder,helOrder=>glob(id)%helOrder )
  Ntot = Ntotal
  Noff = Noffshell
  Noff2 = Noff+2
  NZ = NZbos
  glob(id)%offshell = 0
  glob(id)%onshell = 0
  do ii=1,Noff
    glob(id)%offshell(ii) = ii
  enddo

  do ii=Noff2,Ntot-1
    glob(id)%onshell(ii-Noff2+1) = ii
  enddo
  flavor =-999
  Nflavor = 0
  NflavorFinst = 0

  do ii=1,Ntot
    ! checks if the particle has a correct color index
    if (process(ii).lt.-NsizeFlavor.or.NsizeFlavor.lt.process(ii)) then
      write(*,*) 'ERROR in amp4hef: flavor',process(ii),' not defined'
      stop
    endif
    do jj=-NsizeFlavor,NsizeFlavor
      if (process(ii).eq.jj) then
        Nflavor(jj) = Nflavor(jj)+1
        flavor(Nflavor(jj),jj) = ii
        if (ii.gt.2) NflavorFinst(jj) = NflavorFinst(jj)+1
      endif
    enddo
  enddo
!

  glob(id)%symFac = NcolDof(process(1))*NcolDof(process(2))
  do ii=-NsizeFlavor,NsizeFlavor
    glob(id)%symFac = glob(id)%symFac * factorial(NflavorFinst(ii))
  enddo 

!
if(NZ.eq.1) then
    glob(id)%matrix_element => matrix_element_DrellYan
    glob(id)%all_amplitudes => all_amplitudes_DrellYan
    glob(id)%amplitude      => amplitude_DrellYan
endif

!
  if (mod(sum(process(1:Ntot)), 2).ne.0) then
    write(*,*) 'ERROR in amp4hef: process not possible'
    stop
  endif
  Nsum = 0
!  do jj=-NsizeFlavor,NsizeFlavor
!    if (Nflavor(jj).gt.0) Nsum = Nsum+1
!    if (jj.ne.0.and.Nflavor(jj).gt.1) call not_implemented
!  enddo
!  if (Nsum.gt.3) call not_implemented
!  if (Nsum.eq.1.and.Ntot.gt.7) call not_implemented
!  if (Nsum.gt.1.and.Ntot.gt.5) call not_implemented
!  if (Noff.gt.2) call not_implemented
!  if (Noff.gt.1.and.Nsum.gt.1.and.Ntot.gt.4) call not_implemented
!
! Translation array for helicity ordering:
! first on-shell gluons, then anti-quarks. Quarks automatically get the
! opposite helicity of the anti-quark. Off-shell (anti)-quarks have helicity!

!what does it do? need to rewrite it to adjust to 3 possible polarizations of massive bosons
  NhelOrder = 0
  helOrder = 0
  do jj=Noff2,Ntot
    if (process(jj).ne.0) cycle
    NhelOrder = NhelOrder+1
    helOrder(NhelOrder) = jj
  enddo
  do ii=1,NsizeFlavor
    do jj=1,Ntot
      if (process(jj).ne.-ii) cycle
      NhelOrder = NhelOrder+1
      helOrder(NhelOrder) = jj
    enddo
  enddo
!
  end associate

  contains

    subroutine not_implemented
    write(*,*) 'ERROR in amp4hef: process not implemented'
    stop
    end subroutine

  end subroutine


  subroutine put_momenta( id ,momenta ,directions )
! If there are off-shell momenta, it must be the first or the first two.
  integer,intent(in) :: id
  real(fltknd),intent(in) :: momenta(0:3,*) ,directions(0:3,*)
  integer :: ii,Noff2
  associate( Ntot=>glob(id)%Ntot ,Noff=>glob(id)%Noff ,NZ=>glob(id)%NZ)
  Noff2 = Noff+1
  do ii=1,Noff
    call glob(id)%Q(ii)%fill( momenta(0:3,ii) ,directions(0:3,ii) )
  enddo
  do ii=Noff2,(Ntot-1)
    call glob(id)%Q(ii)%fill( momenta(0:3,ii) )
  enddo
  call glob(id)%Q(Ntot)%fill( momenta(0:3,Ntot), directions(0:3,Ntot))
  do ii=1,Noff
    glob(id)%Q(ii)%kstr = glob(id)%ang(ii,ii,Noff2)/glob(id)%sqr(ii,Noff2)
    glob(id)%Q(ii)%kapp = glob(id)%ang(Noff2,ii,ii)/glob(id)%ang(Noff2,ii)
  enddo
  end associate

  end subroutine


  subroutine matrix_element_a( id ,rslt )
  integer,intent(in) :: id
  real(fltknd),intent(out) :: rslt
  rslt = glob(id)%matrix_element()
  end subroutine


  subroutine matrix_element_b( id ,rslt )
  integer,intent(in) :: id
  real(fltknd),intent(out) :: rslt
  rslt = glob(id)%matrix_element() / glob(id)%symFac
  end subroutine


  subroutine all_amplitudes( id ,NhelConf ,Nperm ,amplitude ,factor )
  integer,intent(in) :: id
  integer,intent(out) :: NhelConf,Nperm
  complex(fltknd),intent(out) :: amplitude(:,:),factor
  call glob(id)%all_amplitudes( NhelConf ,Nperm ,amplitude ,factor )
  factor = factor/glob(id)%symFac
  end subroutine

  subroutine all_amplitudes_cpp( id ,NhelConf ,Nperm ,amplitude ,factor )
  integer,intent(in) :: id
  integer,intent(out) :: NhelConf,Nperm
  complex(fltknd),intent(out) :: amplitude(64,120),factor
  call glob(id)%all_amplitudes( NhelConf ,Nperm ,amplitude ,factor )
  factor = factor/glob(id)%symFac
  end subroutine


  subroutine amplitude( id ,rslt ,helicity ,perm, type )
! Helicities should follow the process as given to put_process,
! so for example helicities refering to off-shell gluons are ignored.
  integer,intent(in) :: id ,helicity(:) ,perm(:), type
  complex(fltknd),intent(out) :: rslt
  associate( o=>glob(id) )
  rslt = o%amplitude( helicity(o%helOrder(1:o%NhelOrder)) ,perm, type )
  end associate
  end subroutine

  subroutine amplitude_cpp( id ,rslt ,helicity ,perm, type )
  integer,intent(in) :: id ,helicity(NsizeProc) ,perm(NsizeProc), type
  complex(fltknd),intent(out) :: rslt
  associate( o=>glob(id) )
  rslt = o%amplitude( helicity(o%helOrder(1:o%NhelOrder)) ,perm, type )
  end associate
  end subroutine


end module



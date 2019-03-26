program check_amplitude
  use amp4hef
  implicit none
  
  integer,parameter :: eventUnit=21
  real(fltknd),parameter :: pi=3.1415926535897932384626433832795_fltknd
  integer,parameter :: lineLen=144
  
  integer :: Noffshell,Nfinst,Ntotal,ii,process(13),helicity(13),perm(13)
  integer :: id1,nn
  real(fltknd) :: momenta(0:3,8),directions(0:3,2)
  character(lineLen) :: eventFile,line,frmt
  complex(fltknd) :: xAmp,amp
  character(1),parameter :: arab(0:9)=['0','1','2','3','4','5','6','7','8','9']
  
  call get_command_argument(1,eventFile)

  open(eventUnit,file=trim(eventFile),status='old')
  do
    read(eventUnit,'(A)',end=999) line
    if (line(1:10).eq.'BEGIN FILE') exit
  enddo
   
  do
    do
      read(eventUnit,'(A)',end=999) line
      if (line(1:11).eq.'BEGIN EVENT') exit
      if (line(1:11).eq.'END FILE') goto 999
    enddo
  
    read(eventUnit,*) Noffshell ,Nfinst
    Ntotal = Nfinst+2

! Read the process. 0=gluon ,positive=quark ,negative=anti-quark.
! All flavors are outgoing, ie sum(process(1:Nfinst+2)) must be 0.
    read(eventUnit,*) process(1:Ntotal)
    frmt = '('//arab(Ntotal)//'i3,2x,'//arab(Ntotal-2)//'i0,2e23.15)'

    call put_process( id1 ,Ntotal ,Noffshell ,process )

    do ii=1,Ntotal
      read(eventUnit,*) momenta(0:3,ii)
    enddo
    directions(0    ,1:2) = momenta(0,1:2)
    directions( 1:2 ,1:2) = 0
    directions(    3,1:2) = momenta(3,1:2)
    call put_momenta( id1 ,momenta ,directions )

    do
      read(eventUnit,'(A)',end=999) line
      if (line(1:9).eq.'END EVENT') exit
      if (line(1:11).eq.'END FILE') goto 999
      if (line(1:8).eq.'COMMENT:') then
        write(*,'(A)') line(9:80)
      else
        read(line,*) helicity(1:Ntotal),perm(1:Ntotal-2)
        call amplitude( id1 ,amp ,helicity ,perm )
        write(*,trim(frmt)) helicity(1:Ntotal),perm(1:Ntotal-2),amp
      endif
    enddo

  enddo

  999 continue
  close(eventUnit)

end program



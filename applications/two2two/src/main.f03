program main
  use amp4hef
  implicit none
  integer,parameter :: readUnit=21,writeUnit=22,lineLen=256
  character(lineLen) :: line,inputFile,outputFile
  integer :: Nwords,iBgn(lineLen),iEnd(lineLen)
  real(fltknd) :: Ecm,k1T,phi1,k3T,y3,phi3 ,y4,phi4,k4T ,rslt
  real(fltknd) :: DyMin,DyMax,Dy ,DphiMin,DphiMax,Dphi ,DyStp,DphiStp
  real(fltknd) :: kk(0:3,1:4),pp(0:3,1:2)
  integer :: NDy,NDphi,ii,jj,process(4),procID
  integer,parameter :: O=0,X=1,Y=2,Z=3

  if (command_argument_count().le.0) then
    write(*,*)
    write(*,*) 'Execute as:'
    write(*,*) '$ ./main.out <inputfile>'
    write(*,*)
    stop
  endif

  call get_command_argument( 1 ,inputFile )
  open(readUnit,file=trim(inputFile),status='old')  
  do
    read(readUnit,'(A)',end=111) line
    call split_line( Nwords ,iBgn ,iEnd ,line )
    if     (line(iBgn(1):iEnd(1)).eq.'k1T')  then ;read(line(iBgn(3):iEnd(3)),*) k1T
    elseif (line(iBgn(1):iEnd(1)).eq.'phi1') then ;read(line(iBgn(3):iEnd(3)),*) phi1
    elseif (line(iBgn(1):iEnd(1)).eq.'k3T')  then ;read(line(iBgn(3):iEnd(3)),*) k3T
    elseif (line(iBgn(1):iEnd(1)).eq.'y3')   then ;read(line(iBgn(3):iEnd(3)),*) y3
    elseif (line(iBgn(1):iEnd(1)).eq.'phi3') then ;read(line(iBgn(3):iEnd(3)),*) phi3
    elseif (line(iBgn(1):iEnd(1)).eq.'k4T')  then ;read(line(iBgn(3):iEnd(3)),*) k4T
    elseif (line(iBgn(1):iEnd(1)).eq.'DyMin') then
      read(line(iBgn(5):iEnd(7)),*) DyMin,DyMax,NDy
    elseif (line(iBgn(1):iEnd(1)).eq.'DphiMin') then
      read(line(iBgn(5):iEnd(7)),*) DphiMin,DphiMax,NDphi
    elseif (line(iBgn(1):iEnd(1)).eq.'output') then
      read(line(iBgn(4):iEnd(4)),*) outputFile
    elseif (line(iBgn(1):iEnd(1)).eq.'process') then
      read(line(iBgn(3):iEnd(6)),*) process 
    endif
  enddo
  111 continue
  close(readUnit)

  call put_process( procID ,4,2 ,process )

  open(writeUnit,file=trim(outputFile),status='new')
  DyStp = (DyMax-DyMin)/(NDy-1)
  DphiStp = (DphiMax-DphiMin)/(NDphi-1)
  do ii=0,NDy-1
    do jj=0,NDphi-1
      Dy = DyMin + ii*DyStp
      Dphi = DphiMin + jj*DphiStp
      if (Dy.eq.0.and.Dphi.eq.0) then
        Dy = Dy + DyStp/2
        Dphi = Dphi + DphiStp/2
      endif
      y4 = y3 + Dy
      phi4 = phi3 + Dphi
      kk(O:Z,3) = k3T*[cosh(y3),cos(phi3),sin(phi3),sinh(y3)]
      kk(O:Z,4) = k4T*[cosh(y4),cos(phi4),sin(phi4),sinh(y4)]
      kk(O,1) =-( (kk(O,3)+kk(Z,3)) + (kk(O,4)+kk(Z,4)) )/2
      kk(O,2) =-( (kk(O,3)-kk(Z,3)) + (kk(O,4)-kk(Z,4)) )/2
      kk(Z,1) = kk(O,1)
      kk(Z,2) =-kk(O,2)
      kk(X:Y,1) =-k1T*[cos(phi1),sin(phi1)]
      kk(X:Y,2) =-kk(X:Y,1)-kk(X:Y,3)-kk(X:Y,4)
      pp(O    ,1:2) = kk(O,1:2)
      pp(    Z,1:2) = kk(Z,1:2)
      pp( X:Y ,1:2) = 0
      call put_momenta( procID ,kk ,pp )
      call matrix_element_a( procID ,rslt )
      write(writeUnit,'(3e16.8)') Dy,Dphi,rslt
    enddo
  enddo
  close(writeUnit)

end program

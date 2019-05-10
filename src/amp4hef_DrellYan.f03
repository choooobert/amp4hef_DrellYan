module amp4hef_DrellYan
  use amp4hef_io
  use amp4hef_aux
  use amp4hef_qomentum
  implicit none
  private
  public :: fill_matrices_DrellYan,matrix_element_DrellYan ,amplitude_DrellYan ,all_amplitudes_DrellYan

  real, parameter :: sqrt_2 = 1.41421356237_fltknd
  real(fltknd) :: MZ
  real(fltknd) :: MZ_sq

  real, parameter :: cV = 0.203666_fltknd
  real, parameter :: cA = 0.5_fltknd

  integer,parameter :: gluon=0 ,quark=1 ,antiq=-1, Zboson=2
     integer,parameter :: helTable_DrellYan(3,6)=reshape(&
     [ -1,1,-1,  -1, 1, 0,  -1, 1, 1,&
       1,-1,-1,   1,-1, 0,   1,-1, 1 &
     ], [3,6])

  integer,allocatable,save :: mtx_4_sqr(:,:)

contains

  function matrix_element_DrellYan(Tin) result(rslt)
  class(qomentum_list_type),intent(in) :: Tin
  real(fltknd) :: rslt
  complex(fltknd) :: amp(12, 2, 2)
  integer :: ii,NhelSum,Nminus2, NhelConf, Nperm, jj,kk
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )

    MZ_sq = square(Tin%Q(Ntot)%k)
    MZ    = sqrt(MZ_sq)
    if (MZ_sq.lt.1) then
        MZ = 0
        MZ_sq =0
    end if

    NhelConf = 6
    NhelSum = 3
    if(Ntot.eq.5) then
        NPerm = 2
    else if(Ntot.eq.4) then
        NPerm = 1
    end if

  rslt = 0
    do ii=1, (NhelConf/2)
      do jj=1, Nperm
        if (Ntot.eq.4) then
          amp(ii, jj, 1) = amplitude_DrellYan(Tin ,helTable_DrellYan(:,ii), [1], 1)
        else if (Ntot.eq.5) then
          do kk=1,2
            amp(ii, kk, jj) = amplitude_DrellYan(Tin ,helTable_DrellYan(:,ii), perTable(1:NPerm, jj), kk)
          enddo
        endif
      enddo
    rslt = rslt + colorSum(Ntot, amp(ii,:, :))
  enddo
  do ii=1, Noff
    rslt=rslt*(Tin%Q(ii)%kapp*Tin%Q(ii)%kstr)
  enddo
  rslt = 2*rslt
  end associate
  end function 


  subroutine all_amplitudes_DrellYan(Tin ,NhelConf ,Nperm ,amplitude ,factor ) ! not relevant
  class(qomentum_list_type),intent(in) :: Tin
  integer,intent(out) :: NhelConf, Nperm
  complex(fltknd),intent(out) :: amplitude(:,:),factor
  integer :: ii,jj,NhelSum
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )
    NhelSum = 3
    NhelConf= 6
    if(Ntot.eq.5) then
        NPerm = 2
    else if(Ntot.eq.4) then
        NPerm = 1
    end if
  do ii=1,(NhelConf/2)
  do jj=1,Nperm
!       amplitude(ii,jj) =(cV+cA)*amplitude_DrellYan(Tin ,helTable_DrellYan(1:NhelSum,ii) ,perTable(1:Nperm,jj) )
!        amplitude(NhelConf-ii+1,jj) =(cV-cA)*conjg(amplitude(ii,jj))
  enddo
  enddo

  factor = 2 ! only half of the helicity configurations is returned
  end associate
  end subroutine


 function amplitude_DrellYan(Tin ,helicity ,perm, type ) result(rslt)
  class(qomentum_list_type),intent(in) :: Tin
  integer,intent(in) :: helicity(:), perm(:), type
  complex(fltknd) :: rslt
  integer ::  i2, i3, i4, j2, j3, j4
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff ,offshell=>Tin%offshell )
    !
    rslt = 0
    j2=1 ;j3=3; j4=2 ! indexing for helicity configuration
    i2=4 ;i3=5; i4=3 ! indexing for qomentum class
    if(Ntot.eq.5) then
      if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.-1) then
        call Tin%set_direction(i3,i4)
        if(type.eq.1) then
          rslt = amp_201(Tin, perm) + amp_207(Tin, perm)
        else if(type.eq.2) then
          rslt = amp_237(Tin, perm)
        end if
      else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.0) then
!        call Tin%set_direction(i3,i4)
        rslt = 0
      else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.1) then
!        call Tin%set_direction(i3,i2)
!        if(type.eq.1) then
!          rslt = amp_203(Tin, perm) + amp_215(Tin, perm)
!        else if(type.eq.2) then
!          rslt = amp_245(Tin, perm)
!        end if
        rslt = 0
      end if

    else if(Ntot.eq.4) then
      if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.-1) then
!        rslt = amp_101(Tin)
      else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.0) then
        rslt = 0
      else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.1) then
        rslt = amp_103(Tin)
      end if

!        if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.-1) then
!        rslt = amp_201(Tin)
!        else if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.0) then
!        rslt = 0
!        else if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.1) then
!        rslt = amp_203(Tin)
!        end if
    end if

  end associate
end function

!! !amplitudes for 1-jet process
    function amp_101(T) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,vv ,xx,yy, zz
    integer :: i1, i2, i3, i4
    !
    i1=1 ;i2=2 ;i3=4; i4=3
    !
    rslt = 0
    call T%set_direction(i3,i4)
    vv =T%sqr(i4,i1)*T%sqr(i4,i1)*T%ang(i3,i2)
    yy = (MZ_sq+T%ang(i2,i3,i2))*T%Q(i1)%kapp*T%sqr(i4,i3)
    rslt = 2*vv/yy
    end function

    function amp_102(T) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,vv ,xx,yy, zz
    integer :: i1, i2, i3, i4
    !
    i1=1 ;i2=2 ;i3=4; i4=3
    !
    rslt = 0
    call T%set_direction(i3,i1)
    xx = (T%sqr(i4,i3)*T%ang(i3,i4)/MZ+MZ*T%ang(i1,i2)*T%sqr(i2,i1)/(T%ang(i3,i1)*T%sqr(i1,i3)))/(MZ*MZ+T%ang(i4,i3,i4))
    yy = (T%sqr(i2,i3)*T%ang(i3,i2)/MZ+MZ*T%ang(i1,i4)*T%sqr(i4,i1)/(T%ang(i3,i1)*T%sqr(i1,i3)))/(MZ*MZ+T%ang(i2,i3,i2))
    zz = T%sqr(i4,i1)*T%ang(i1,i2)/(T%Q(i1)%kapp*T%Q(i1)%kstr)
    rslt = sqrt_2*(-xx+yy)/zz
    end function


    function amp_103(T) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,vv ,xx,yy, zz
    integer :: i1, i2, i3, i4
    !
    i1=1 ;i2=2 ;i3=4; i4=3
    !
    rslt = 0
    call T%set_direction(i3,i2)
    vv =T%ang(i2,i1)*T%ang(i2,i1)*T%sqr(i3,i4)
    yy = (MZ_sq+T%ang(i4,i3,i4))*T%Q(i1)%kstr*T%ang(i2,i3)
    rslt = 2*vv/yy
    end function


!!!!! amplitudes for 2-jet process
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ! 1st diagram amplitudes
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    function amp_201(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    rslt = 0
    write(*,*) "amp_201"
    xx= T%sqr(i4,i5)*T%sqr(i4,i5)* T%ang(i1,i2)/T%sqr(i4,i3)
    yy = T%ang(i3,i1) + T%ang(i3,i2)*T%sqr(i2,i1)/ T%Q(i1)%kapp
    zz = (square(T%Q(i1)%k)+T%ang(i2,i1,i2))*(square(T%Q(i5)%k)+T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
    rslt = 2*sqrt_2*xx*yy/zz
    end function
!
!   function amp_202(T, perm) result(rslt)
!   type(qomentum_list_type),intent(in) :: T
!   integer, intent(in) :: perm(:)
!   complex(fltknd) :: rslt,tt,uu,vv,xx,ww,yy,zz
!   integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!   rslt = 0
!   call T%set_direction(i3, i4)
!   tt = T%sqr(i4,i5)*T%ang(i1,i2)
!   uu = (-T%Q(i1)%kapp*T%Q(i1)%kstr+T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
!   ww = T%sqr(i5,i3) + T%ang(i5,i4)*T%sqr(i4,i3)/T%Q(i5)%kstr
!   xx = T%ang(i3,i1) + T%ang(i3,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
!   yy = T%sqr(i5,i4)/(T%ang(i3,i4)*T%sqr(i4,i3))
!   zz = T%ang(i4,i1) + T%ang(i4,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
!   rslt = -2*tt/uu*(-ww*xx/MZ+MZ*yy*zz)
!   end function
!
    function amp_203(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    rslt = 0
    write(*,*) "amp_203"

    xx = T%sqr(i4,i5)*T%ang(i1,i2)*T%ang(i1,i2)/T%ang(i2,i3)
    yy = T%sqr(i5,i3)+T%ang(i5,i4)*T%sqr(i4,i3)/T%Q(i5)%kstr
    zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr+T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
    rslt = 2*sqrt_2*xx*yy/zz
    end function
!
!
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 2nd diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    function amp_207(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz, vv
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    rslt = 0
    write(*,*) "amp_207"
    vv = T%sqr(i4,i5)*T%ang(i3,i2)
    xx = T%sqr(i5,i1) + T%ang(i5,i4)*T%sqr(i4,i1)/ T%Q(i5)%kstr
    yy = T%ang(i1,i3) + T%ang(i1,i2)*T%sqr(i2,i4)/T%sqr(i3,i4)
    zz = (MZ_sq + T%ang(i2, i3, i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i4,i5,i4)) &
       * T%Q(i1)%kapp*T%Q(i1)%kstr*T%Q(i5)%kapp
    rslt = 2*sqrt_2*vv*xx*yy/zz
    end function
!
!    function amp_208(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,ww,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    call T%set_direction(i3, i2)
!    rslt = 0
!    call T%set_direction(i3, i4)
!    ww = T%ang(i1,i2)*T%sqr(i4,i5)
!    xx = T%sqr(i5,i1) - T%ang(i5,i4)*T%sqr(i4,i1)/T%Q(i5)%kstr
!    yy = MZ + T%sqr(i2,i3)*T%ang(i3,i2)/MZ
!    zz = (MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i5)%kapp*T%Q(i1)%kapp*T%Q(i1)%kstr
!    if (ww.ne.0.and.xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*ww*xx*yy/zz
!    end function
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 3rd diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    function amp_214(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,ww,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    call T%set_direction(i3, i4)
!    rslt = 0
!    call T%set_direction(i3, i4)
!    ww = T%ang(i2,i1)*T%sqr(i4,i5)
!    xx = T%ang(i5,i1) - T%ang(i5,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
!    yy = MZ + T%sqr(i4,i3)*T%ang(i3,i4)/MZ
!    zz = (MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i5)%kapp*T%Q(i5)%kstr*T%Q(i1)%kstr
!    if (ww.ne.0.and.xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*ww*xx*yy/zz
!    end function
!
    function amp_215(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,yy,zz, vv
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    rslt = 0
    write(*,*) "amp_215"
    vv = T%sqr(i3,i4)*T%ang(i1,i2)
    xx = T%sqr(i3,i5) + T%ang(i2,i4)*T%sqr(i4,i5)/T%ang(i2,i3)
    yy = T%ang(i5,i1) + T%ang(i5,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
    zz = (MZ_sq + T%ang(i4, i3, i4))*(-T%Q(i1)%kapp*T%Q(i1)%kstr+T%ang(i2,i1,i2))*T%Q(i5)%kapp*T%Q(i5)%kstr*T%Q(i1)%kstr
    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*vv*xx*yy/zz
    end function

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 7th diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    function amp_237(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx, uu,vv, bb, dd, ff, i
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    i=(0,1)
    rslt = 0
    write(*,*) "amp_237"
    xx = T%ang(i3,i2)/((MZ_sq+ T%ang(i2,i3,i2)+ T%ang(i4,i3,i4) +T%sqr(i2,i4)*T%ang(i4,i2))*(MZ_sq+T%ang(i2,i3,i2))&
        * T%Q(i5)%kapp*T%Q(i5)%kstr*T%Q(i1)%kapp*T%Q(i1)%kstr)
    uu = T%ang(i1,i5)*T%sqr(i5,i1)/2
    vv = T%sqr(i4,i1,i3) + T%sqr(i2,i4)/T%sqr(i3,i4)*T%sqr(i4,i1,i2)
    bb = T%ang(i5,i1,i5) - T%ang(i1,i5)*T%sqr(i5,i1)*T%Q(i1)%kapp*T%Q(i1)%kstr/T%ang(i1,i5,i1)
    dd = T%ang(i1,i3) + T%ang(i1,i2)*T%sqr(i2,i4)/T%sqr(i3,i4)
    ff = T%sqr(i4,i1)
    rslt =i*2*xx*(uu*vv - ff*bb*dd)
    end function
!
!    function amp_238(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    call T%set_direction(i3, i2)
!    xx= T%ang(i1,i5)*T%sqr(i5,i1)*T%sqr(i4,i5)*T%ang(i5,i2)
!    yy =T%ang(i3,i2)*T%sqr(i2,i3)/MZ+MZ
!    zz = T%Q(i1)%kapp*T%Q(i1)%kstr*(MZ*MZ+T%ang(i2,i3,i2))*T%ang(i5,i1,i5) &
!       *(MZ*MZ+T%ang(i2,i4)*T%sqr(i4,i2)+T%ang(i2,i3,i2)+T%ang(i4,i3,i4))
!    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = sqrt_2*xx*yy/zz
!    end function
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 8th diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    function amp_244(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    call T%set_direction(i3, i4)
!    xx= T%ang(i1,i5)*T%sqr(i5,i1)*T%sqr(i4,i5)*T%ang(i5,i2)
!    yy =T%ang(i3,i4)*T%sqr(i4,i3)/MZ+MZ
!    zz = T%Q(i1)%kapp*T%Q(i1)%kstr*(MZ*MZ+T%ang(i4,i3,i4))*T%ang(i5,i1,i5) &
!       *(MZ*MZ+T%ang(i2,i4)*T%sqr(i4,i2)+T%ang(i2,i3,i2)+T%ang(i4,i3,i4))
!    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = -sqrt_2*xx*yy/zz
!    end function
!
  function amp_245(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx, uu,vv, bb, dd, ff, i
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    i=(0,1)
    rslt = 0
    write(*,*) "amp_245"
    xx = T%sqr(i4,i3)/((MZ_sq + T%ang(i2,i3,i2)+ T%ang(i4,i3,i4)+T%ang(i2,i4)*T%sqr(i4,i2))*(MZ_sq + T%ang(i4,i3,i4)) &
       * T%Q(i5)%kapp*T%Q(i5)%kstr*T%Q(i1)%kapp*T%Q(i1)%kstr)

    uu = T%ang(i1,i5)*T%sqr(i5,i1)/2
    vv = T%sqr(i3,i1,i2) + T%ang(i2,i4)/T%ang(i3,i4)*T%ang(i4,i1,i2)

    bb = -T%ang(i1,i5,i1) + T%ang(i1,i5)*T%sqr(i5,i1)*T%Q(i5)%kapp*T%Q(i5)%kstr/T%ang(i5,i1,i5)
    dd = T%sqr(i3,i5) + T%ang(i2,i4)*T%sqr(i4,i5)/T%ang(i3,i4)
    ff = T%ang(i5,i2)
    rslt = i*2*sqrt_2*xx*(uu*vv + ff*bb*dd)
  end function


  subroutine fill_matrices_DrellYan
  integer :: rUnit
  call check_file('amp4hef.tbl')
  open(newunit=rUnit,file=trim(get_path('amp4hef.tbl')),status='old')
  allocate(mtx_4_sqr(36,24))
  call read_matrix( mtx_4_sqr ,'BEGIN qq4g square' )
  close(rUnit)
!
  contains
!
    subroutine read_matrix( matrix ,tag )
    integer,intent(out) :: matrix(:,:)
    character(*),intent(in) :: tag
    character(144) :: line
    integer :: ii
    do
      read(rUnit,'(A)') line
      if (line(1:len(tag)).eq.tag) exit
    enddo
    do ii=1,size(matrix,2)
      read(rUnit,*) matrix(1:size(matrix,1),ii)
    enddo
    end subroutine
!
  end subroutine


function colorSum( Ntot ,t ) result(rslt)
  integer,intent(in) :: Ntot
  complex(fltknd),intent(in) :: t(2,2)
  complex(fltknd) :: rslt ,z(2,2), a, b, c, az, bz, cz
  real(fltknd) :: AA, BB, CC, AB, AC, BC
  integer :: Nadj,i
!
  Nadj = Ncolor(2)-1
!
  AA = 256./3. ; BB = 256./3. ; CC = 2*384.
  AB = -32./3. ; AC = -192.   ; BC = 192.

  select case (Ntot)
  case (4)
    rslt = conjg(t(1,1))*t(1, 1) * Nadj

  case (5)
    z(1:2, 1:2) = conjg(t(1:2, 1:2))
    a = t(1,1)
    az= z(1,1)
    b = t(1,2)
    bz= z(1,2)
    c = t(2,1)-t(2,2)
    cz= z(2,1)-z(2,2)
    rslt = AA*a*az + BB*b*bz + CC*c*cz &
         + AB*(a*bz + b*az) &
         + AC*(a*cz + c*az) &
         + BC*(b*cz + c*bz)
  case default
    rslt = 0
  end select
end function


end module




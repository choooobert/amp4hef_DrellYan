module amp4hef_DrellYan
  use amp4hef_io
  use amp4hef_aux
  use amp4hef_qomentum
  implicit none
  private
  public :: fill_matrices_DrellYan,matrix_element_DrellYan ,amplitude_DrellYan ,all_amplitudes_DrellYan, impact_factors_all

  real, parameter :: sqrt_2 = 1.41421356237_fltknd
  real(fltknd) :: MZ
  real(fltknd) :: MZ_sq

  real, parameter :: cV = 0.203666_fltknd
  real, parameter :: cA = 0.5_fltknd

  integer,parameter :: gluon=0 ,quark=1 ,antiq=-1, Zboson=2
     integer,parameter :: helTable_DrellYan(3,6)=reshape(&
     [-1, 1,-1,  -1, 1, 1,   -1, 1, 0,&
       1,-1, 1,   1,-1,-1,    1,-1, 0 &
     ], [3,6])

  integer,allocatable,save :: mtx_4_sqr(:,:)

contains

  function impact_factors_all(Tin) result(rslt)
      class(qomentum_list_type),intent(in) :: Tin
      real(fltknd) :: rslt(9)
  complex(fltknd) :: amp(12, 2)
  integer :: ii,Nminus2, NhelConf, Nperm, jj,kk
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )

    MZ_sq = square(Tin%Q(Ntot)%k)
    MZ    = sqrt(MZ_sq)

    if((MZ.lt.1E-3).and.(MZ.gt.-1E-3)) then
      NhelConf = 4
      rslt(5:9) = [0, 0, 0, 0, 0]
    else
      NhelConf = 6
    endif

    if(Ntot.eq.5) then
        NPerm = 2
    else if(Ntot.eq.4) then
        NPerm = 1
    end if

    rslt = 0
    do ii=1, (NhelConf/2)
      do jj=1, Nperm
          amp(ii, jj) = amplitude_DrellYan(Tin ,helTable_DrellYan(:,ii), perTable(1:NPerm, jj))
      enddo
  enddo

  kk =1
  do ii=1, (Nhelconf/2)
    do jj=1, (Nhelconf/2)
      rslt(kk) = 2*colorSum(Ntot, amp(ii,:), amp(jj,:))
    enddo
  enddo
  end associate
  end function



  function matrix_element_DrellYan(Tin) result(rslt)
  class(qomentum_list_type),intent(in) :: Tin
  real(fltknd) :: rslt
  complex(fltknd) :: amp(12, 2)
  integer :: ii,Nminus2, NhelConf, Nperm, jj,kk
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )

    MZ_sq = square(Tin%Q(Ntot)%k)
    MZ    = sqrt(MZ_sq)

    if((MZ.lt.1E-3).and.(MZ.gt.-1E-3)) then
      NhelConf =4
    else
      NhelConf = 6
    endif

    if(Ntot.eq.5) then
        NPerm = 2
    else if(Ntot.eq.4) then
        NPerm = 1
    end if

    rslt = 0
    do ii=1, (NhelConf/2)
      do jj=1, Nperm
          amp(ii, jj) = amplitude_DrellYan(Tin ,helTable_DrellYan(:,ii), perTable(1:NPerm, jj))
      enddo
    rslt = rslt + colorSum(Ntot, amp(ii,:), amp(ii,:))
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


 function amplitude_DrellYan(Tin ,helicity ,perm ) result(rslt)
  class(qomentum_list_type),intent(in) :: Tin
  integer,intent(in) :: helicity(:), perm(:)
  complex(fltknd) :: rslt, rslt_X, rslt_Y, i
  integer ::  i2, i3, i4, j2, j3, j4, iX, iY, iZ
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff ,offshell=>Tin%offshell )
    !
    i = (0._fltknd,1._fltknd)
    rslt = 0
    j2=1 ;j3=3; j4=2 ! indexing for helicity configuration
    i2=4 ;i3=5; i4=3 ! indexing for qomentum class
    iX = Ntot + 1
    iY = Ntot + 2
    iZ = Ntot + 3
    if(Ntot.eq.5) then

    if (helicity(j3).ne.0) then
      rslt_X =  amp_202(Tin, perm, iX) + amp_205(Tin, perm, iX) + amp_208(Tin, perm, iX) &
              + amp_211(Tin, perm, iX)+  amp_214(Tin, perm, iX)
      rslt_Y =  amp_202(Tin, perm, iY) + amp_205(Tin, perm, iY) + amp_208(Tin, perm, iY) &
              + amp_211(Tin, perm, iY)+  amp_214(Tin, perm, iY)
      if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.-1) then
        rslt = (rslt_X - i*rslt_Y)/sqrt_2
      else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.1) then
        rslt = -(rslt_X + i*rslt_Y)/sqrt_2
      endif
    else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.0) then
      rslt =  amp_202(Tin, perm, iZ) + amp_205(Tin, perm, iZ) + amp_208(Tin, perm, iZ) &
            + amp_211(Tin, perm, iZ)+  amp_214(Tin, perm, iZ)
    end if

    else if(Ntot.eq.4) then
      if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.-1) then
        rslt_X = amp_102(Tin, iX)
        rslt_Y = amp_102(Tin, iY)
        rslt = (rslt_X - i*rslt_Y)/sqrt_2
!        write(*,*) rslt_X, rslt_Y

      else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.0) then
        rslt = amp_102(Tin, iZ)

      else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.1) then
        rslt_X = amp_102(Tin, iX)
        rslt_Y = amp_102(Tin, iY)
        rslt = -(rslt_X + i*rslt_Y)/sqrt_2
!        write(*,*) rslt_X, rslt_Y
      end if
    end if


  end associate
end function

!! !amplitudes for 1-jet process



    function amp_102(T, Z) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: Z
    complex(fltknd) :: rslt,vv ,xx,yy, zz, r1, r2, l2
    complex(fltknd) :: l1, r3, r4, r5, r6, r7, r8, j
    integer :: i1, i2, i3, i4, X, Y
    !
    X = 5
    Y = 6
    j=(0.,1.)
    i1=1 ;i2=2 ;i3=4; i4=3
    !
    rslt = 0
    l2 = square(T%Q(i3)%k) + twodot(T%Q(i2)%p, T%Q(i3)%k)

    l1 = -square(T%Q(i3)%k) - twodot(T%Q(i4)%p, T%Q(i3)%k)
    r3 = T%sqr(i4, X,[i3,i4], Z, i4)
    r4 = T%sqr(i4, Z,[i3,i4], X, i4)

    r5 = T%ang(i1, i2, X, [i2, i3], i1)
    r6 = T%ang(i1, i2, Z, [i2, i3], i1)
    rslt = r3/l1+r4/l2

end function



!!!!! amplitudes for 2-jet process
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ! 1st diagram amplitudes
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function amp_202(T, perm, U) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:), U
    complex(fltknd) :: rslt, xx, zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    rslt = 0
    xx = T%sqr(i4,i5)*T%ang(i5,[i5,i4],U,[i1,i2], i1)* T%ang(i1,i2)
    zz = (square(T%Q(i1)%k)+T%ang(i2,i1,i2))*(square(T%Q(i5)%k)+T%ang(i4,i5,i4))
    rslt = -xx/zz
    end function

!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 2nd diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

    function amp_205(T, perm,U) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:), U
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    rslt = 0
    xx = T%sqr(i4,i5)*T%ang(i5,[i5,i4],i1)*T%ang(i1,[i2,i3],U, i2)
    zz = (MZ_sq+T%ang(i2,i3,i2))*(square(T%Q(i5)%k)+T%ang(i4,i5,i4))
    rslt = -xx/zz
    end function

!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 3rd diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

    function amp_208(T, perm, U) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:), U
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    rslt = 0
    xx = T%sqr(i4, U,[i3,i4], i5)*T%ang(i5,[i1,i2],i1) *T%ang(i1,i2)
    zz = (square(T%Q(i1)%k)+T%ang(i2,i1,i2))*(MZ_sq+T%ang(i4,i3,i4))
    rslt = -xx/zz
    end function

!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 7th diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

    function amp_211(T, perm, U) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx, uu,vv, bb, dd, ff, gg
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:), U
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    rslt = 0
    xx = 1/((square(T%Q(i1)%k) +square(T%Q(i5)%k) +twodot(T%Q(i1)%k, T%Q(i5)%k))*(MZ_sq+T%ang(i2,i3,i2)))
    uu =T%ang(i1,i5)*T%sqr(i5,i1)/2
    vv = T%sqr(i4,i5,[i2,i3], U, i2) - T%sqr(i4,i1,[i2,i3], U, i2)

    bb = T%ang(i5,i1,i5)&
       + T%ang(i1,i5)*T%sqr(i5,i1)*square(T%Q(i1)%k)/T%ang(i1,i5,i1)
    dd = T%sqr(i4,i1)*T%ang(i1,[i2,i3], U, i2)

    ff = T%ang(i1,i5,i1)&
       + T%ang(i1,i5)*T%sqr(i5,i1)*square(T%Q(i5)%k)/T%ang(i5,i1,i5)
    gg = T%sqr(i4,i5)*T%ang(i5,[i2,i3], U, i2)

    rslt = xx*(uu*vv + bb*dd - ff*gg)
    end function

!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 8th diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


    function amp_214(T, perm, U) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx, uu,vv, bb, dd, ff, gg
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:), U
    !
    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
    !
    rslt = 0
    xx = 1/((square(T%Q(i1)%k) +square(T%Q(i5)%k) +twodot(T%Q(i1)%k, T%Q(i5)%k))*(MZ_sq+T%ang(i4,i3,i4)))
    uu = T%ang(i1,i5)*T%sqr(i5,i1)/2
    vv = T%sqr(i4,U,[i3,i4], i1, i2) - T%sqr(i4,U,[i3,i4], i5, i2)

    bb = T%ang(i5,i1,i5)&
       + T%ang(i1,i5)*T%sqr(i5,i1)*square(T%Q(i1)%k)/T%ang(i1,i5,i1)
    dd = T%sqr(i4,U,[i3,i4], i1)*T%ang(i1, i2)

    ff = T%ang(i1,i5,i1)&
       + T%ang(i1,i5)*T%sqr(i5,i1)*square(T%Q(i5)%k)/T%ang(i5,i1,i5)
    gg = T%sqr(i4,U,[i3,i4], i5)*T%ang(i5, i2)

    rslt =xx*(uu*vv - bb*dd + ff*gg)
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


function colorSum( Ntot ,t, u ) result(rslt)
  integer,intent(in) :: Ntot
  complex(fltknd),intent(in) :: t(2), u(2)
  complex(fltknd) :: rslt ,z(2), uu, vv
  integer :: Nadj,i
!
  Nadj = Ncolor(2)-1
!
  select case (Ntot)
  case (4)
    rslt = t(1)*conjg(u(1)) * Nadj
  case (5)
    z = conjg(u(1:2))
    uu = z(1)*t(1) + z(2)*t(2)
    vv =-z(1)*t(2) - z(2)*t(1)
    rslt = ( vv + Nadj*uu ) * Nadj/Ncolor(1)
  end select
end function


end module

!    function amp_101(T) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,vv ,xx,yy, zz
!    integer :: i1, i2, i3, i4
!    !
!    i1=1 ;i2=2 ;i3=4; i4=3
!    !
!    rslt = 0
!    call T%set_direction(i3,i4)
!    vv =T%sqr(i4,i1)*T%sqr(i4,i1)*T%ang(i3,i2)
!    yy = (MZ_sq+T%ang(i2,i3,i2))*T%Q(i1)%kapp*T%sqr(i4,i3)
!    rslt = 2*vv/yy
!    end function
!
!    function amp_102(T) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,vv ,xx,yy, zz
!    integer :: i1, i2, i3, i4, Z
!    !
!    i1=1 ;i2=2 ;i3=4; i4=3; Z=5
!    !
!    rslt = 0
!    call T%set_direction(i3,i1)
!    xx = T%sqr(i4, Z, [i3, i4], i1)*T%ang(i1,i2)
!    yy = T%sqr(i4,i1)*T%ang(i1,[i3,i2], Z, i2)
!    vv = square(T%Q(i3)%k) + T%ang(i4,i3,i4)
!    zz = square(T%Q(i3)%k) + T%ang(i2,i3,i2)
!    rslt = xx/vv - yy/zz
!    end function
!
!
!    function amp_103(T) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,vv ,xx,yy, zz
!    integer :: i1, i2, i3, i4
!    !
!    i1=1 ;i2=2 ;i3=4; i4=3
!    !
!    rslt = 0
!    call T%set_direction(i3,i2)
!    vv =T%ang(i2,i1)*T%ang(i2,i1)*T%sqr(i3,i4)
!    yy = (MZ_sq+T%ang(i4,i3,i4))*T%Q(i1)%kstr*T%ang(i2,i3)
!    rslt = 2*vv/yy
!    end function




!
!    function amp_201(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,xx,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    integer, intent(in) :: perm(:)
!    !
!    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
!    !
!    rslt = 0
!    xx= T%sqr(i4,i5)*T%sqr(i4,i5)* T%ang(i1,i2)/T%sqr(i4,i3)
!    yy = T%ang(i3,i1) + T%ang(i3,i2)*T%sqr(i2,i1)/ T%Q(i1)%kapp
!    zz = (square(T%Q(i1)%k)+T%ang(i2,i1,i2))*(square(T%Q(i5)%k)+T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
!    rslt = xx*yy/zz
!    end function
!!
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
!    function amp_203(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
!    !
!    rslt = 0
!    xx = T%sqr(i4,i5)*T%ang(i1,i2)*T%ang(i1,i2)/T%ang(i2,i3)
!    yy = T%sqr(i5,i3)+T%ang(i5,i4)*T%sqr(i4,i3)/T%Q(i5)%kstr
!    zz = (square(T%Q(i1)%k)+T%ang(i2,i1,i2))*(square(T%Q(i5)%k)+T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
!    rslt = xx*yy/zz
!    end function
!
!
!    function amp_204(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,xx,yy,zz, vv
!    integer :: i1, i2, i3, i4, i5
!    integer, intent(in) :: perm(:)
!    !
!    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
!    !
!    rslt = 0
!    vv = T%sqr(i4,i5)*T%ang(i3,i2)
!    xx = T%sqr(i5,i1) + T%ang(i5,i4)*T%sqr(i4,i1)/ T%Q(i5)%kstr
!    yy = T%ang(i1,i3) + T%ang(i1,i2)*T%sqr(i2,i4)/T%sqr(i3,i4)
!    zz = (MZ_sq + T%ang(i2, i3, i2))*(square(T%Q(i5)%k)+T%ang(i4,i5,i4)) &
!       * T%Q(i1)%kapp*T%Q(i1)%kstr*T%Q(i5)%kapp
!    rslt = vv*xx*yy/zz
!    end function
!
!    function amp_205(T, perm) result(rslt)
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

!    function amp_208(T, perm) result(rslt)
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
!    function amp_209(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,yy,zz, vv
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
!    !
!    rslt = 0
!    vv = T%sqr(i3,i4)*T%ang(i1,i2)
!    xx = T%sqr(i3,i5) + T%ang(i2,i4)*T%sqr(i4,i5)/T%ang(i2,i3)
!    yy = T%ang(i5,i1) + T%ang(i5,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
!    zz = (MZ_sq + T%ang(i4, i3, i4))*(square(T%Q(i1)%k)+T%ang(i2,i1,i2))*T%Q(i5)%kapp*T%Q(i5)%kstr*T%Q(i1)%kstr
!    rslt = vv*xx*yy/zz
!    end function
!
!    function amp_210(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,xx, uu,vv, bb, dd, ff, gg
!    integer :: i1, i2, i3, i4, i5
!    integer, intent(in) :: perm(:)
!    !
!    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
!    !
!    rslt = 0
!    xx = T%ang(i3,i2)/((MZ_sq+ T%ang(i2,i3,i2)+ T%ang(i4,i3,i4) &
!       + T%sqr(i2,i4)*T%ang(i4,i2))*(MZ_sq+T%ang(i2,i3,i2)))
!    uu = T%ang(i1,i5)*T%sqr(i5,i1)/(2*square(T%Q(i1)%k)*square(T%Q(i5)%k))
!    vv = (T%sqr(i4,i1,i3)-T%sqr(i4,i5,i3))&
!       + T%sqr(i2,i4)/T%sqr(i3,i4)*(T%sqr(i4,i1,i2)-T%sqr(i4,i5,i2))
!
!    bb = -T%ang(i5,i1,i5)/(square(T%Q(i1)%k)*square(T%Q(i5)%k))&
!    - T%ang(i1,i5)*T%sqr(i5,i1)/(T%ang(i1,i5,i1)*square(T%Q(i5)%k))
!    dd = T%sqr(i4,i1)*T%ang(i1,i3) + T%sqr(i4,i1)*T%ang(i1,i2)*T%sqr(i2,i4)/T%sqr(i3,i4)
!
!    ff = +T%ang(i1,i5,i1)/(square(T%Q(i1)%k)*square(T%Q(i5)%k))&
!    + T%ang(i1,i5)*T%sqr(i5,i1)/(T%ang(i5,i1,i5)*square(T%Q(i1)%k))
!    gg = T%sqr(i4,i5)*T%ang(i5,i3) + T%sqr(i4,i5)*T%ang(i5,i2)*T%sqr(i2,i4)/T%sqr(i3,i4)
!
!    rslt = xx*(uu*vv + bb*dd + ff*gg)
!    end function

!
!    function amp_211(T, perm) result(rslt)
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



!    function amp_214(T, perm) result(rslt)
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
!  function amp_215(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,xx, uu,vv, bb, dd, ff, gg
!    integer :: i1, i2, i3, i4, i5
!    integer, intent(in) :: perm(:)
!    !
!    i1=perm(2) ;i2=4 ;i3=5; i4=3; i5=perm(1)
!    !
!    rslt = 0
!    xx = T%sqr(i4,i3)/((MZ_sq+ T%ang(i2,i3,i2)+ T%ang(i4,i3,i4) &
!       + T%sqr(i2,i4)*T%ang(i4,i2))*(MZ_sq+T%ang(i4,i3,i4)))
!    uu = T%ang(i1,i5)*T%sqr(i5,i1)/(2*square(T%Q(i1)%k)*square(T%Q(i5)%k))
!    vv = (T%sqr(i3,i1,i2)-T%sqr(i3,i5,i2))&
!       + T%ang(i2,i4)/T%ang(i2,i3)*(T%sqr(i4,i1,i2)-T%sqr(i4,i5,i2))
!    bb = -T%ang(i5,i1,i5)/(square(T%Q(i1)%k)*square(T%Q(i5)%k))&
!    - T%ang(i1,i5)*T%sqr(i5,i1)/(T%ang(i1,i5,i1)*square(T%Q(i5)%k))
!    dd = T%sqr(i3,i1)*T%ang(i1,i2) + T%sqr(i4,i1)*T%ang(i1,i2)*T%ang(i2,i4)/T%ang(i2,i3)
!
!    ff = T%ang(i1,i5,i1)/(square(T%Q(i1)%k)*square(T%Q(i5)%k))&
!    + T%ang(i1,i5)*T%sqr(i5,i1)/(T%ang(i5,i1,i5)*square(T%Q(i1)%k))
!    gg = T%sqr(i3,i5)*T%ang(i5,i2) + T%sqr(i4,i5)*T%ang(i5,i2)*T%ang(i2,i4)/T%ang(i2,i3)
!
!    rslt = xx*(uu*vv + bb*dd + ff*gg)
!    end function



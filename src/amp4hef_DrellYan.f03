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

    if((MZ_sq.lt.1E-3).and.(MZ_sq.gt.-1E-3)) then
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

    write(*,*) "MZ ", MZ_sq

    if((MZ_sq.lt.1E-3).and.(MZ_sq.gt.-1E-3)) then
      NhelConf = 4
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
    MZ_sq = square(Tin%Q(Ntot)%k)
    i = (0._fltknd,1._fltknd)
    rslt = 0
    j2=1 ;j3=3; j4=2 ! indexing for helicity configuration
    i2=4 ;i3=5; i4=3 ! indexing for qomentum class
    iX = Ntot + 1
    iY = Ntot + 2
    iZ = Ntot + 3
    if(Ntot.eq.5) then
    if (helicity(j3).ne.0) then
      rslt_X =  amp_2jet_1(Tin, perm, iX) + amp_2jet_2(Tin, perm, iX) + amp_2jet_3(Tin, perm, iX) &
              + amp_2jet_7(Tin, perm, iX)+  amp_2jet_8(Tin, perm, iX)
      rslt_Y =  amp_2jet_1(Tin, perm, iY) + amp_2jet_2(Tin, perm, iY) + amp_2jet_3(Tin, perm, iY) &
              + amp_2jet_7(Tin, perm, iY)+  amp_2jet_8(Tin, perm, iY)
      if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.-1) then
        if((MZ_sq.lt.1E-3).and.(MZ_sq.gt.-1E-3)) then
          write(*,*) "I'm in minus"
          rslt = amp_2jet_1m(tin, perm) + amp_2jet_2m(tin, perm) + amp_2jet_7m(Tin, perm)
        else
          write(*,*) "I'm in minus"
          rslt = (rslt_X - i*rslt_Y)/sqrt_2
        endif
      else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.1) then
        if((MZ_sq.lt.1E-3).and.(MZ_sq.gt.-1E-3)) then
          write(*,*) "I'm in plus "
          rslt = amp_2jet_1p(tin, perm) + amp_2jet_3p(tin, perm) + amp_2jet_8p(Tin, perm)
        else
          write(*,*) "I'm in plus "
          rslt = -(rslt_X + i*rslt_Y)/sqrt_2
        endif
      endif
    else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.0) then
       if((MZ_sq.lt.1E-3).and.(MZ_sq.gt.-1E-3)) then
         write(*,*) "I'm in zero  "
         rslt = 0
       else
         write(*,*) "I'm in zero  "
         rslt =  amp_2jet_1(Tin, perm, iZ) + amp_2jet_2(Tin, perm, iZ) + amp_2jet_3(Tin, perm, iZ) &
            + amp_2jet_7(Tin, perm, iZ)+  amp_2jet_8(Tin, perm, iZ)
       endif
    end if

    else if(Ntot.eq.4) then
    if (helicity(j3).ne.0) then
      rslt_X =  amp_1jet(Tin, iX)
      rslt_Y =  amp_1jet(Tin, iY)
      if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.-1) then
        if((MZ_sq.lt.1E-3).and.(MZ_sq.gt.-1E-3)) then
          !rslt = amp_1jet_m(Tin)
        else
          rslt = (rslt_X - i*rslt_Y)/sqrt_2
        endif
      else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.1) then
        if((MZ_sq.lt.1E-3).and.(MZ_sq.gt.-1E-3)) then
          rslt = amp_1jet_p(Tin)
        else
          rslt = -(rslt_X + i*rslt_Y)/sqrt_2
        endif
      endif
    else if (helicity(j2).eq.-1.and.helicity(j4).eq.1.and.helicity(j3).eq.0) then
      rslt =  amp_1jet(Tin, iZ)
    end if

    end if


  end associate
end function

!!!! amplitudes for 1-jet process

    function amp_1jet(T, Z) result(rslt)
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

    function amp_1jet_m(T) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,vv ,xx,yy, zz
    integer :: i1, i2, i3, i4
    !
    i1=1 ;i2=2 ;i3=4; i4=3
    !
    rslt = 0
    vv =T%sqr(i4,i1)*T%sqr(i4,i1)*T%Q(i1)%kstr
    yy = T%sqr(i2, i3)*T%sqr(i3,i4)
    rslt = -vv/yy
    write(*,*) "spinor ", abs(T%sqr(i4,i1))**2

    end function



    function amp_1jet_p(T) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,vv ,xx,yy, zz
    integer :: i1, i2, i3, i4
    !
    i1=1 ;i2=2 ;i3=4; i4=3
    !
    rslt = 0
    vv =T%ang(i1,i2)*T%ang(i1,i2)*T%Q(i1)%kapp
    yy = T%ang(i3,i4)*T%ang(i2,i3)
    rslt = vv/yy
    write(*,*) "spinor ", abs(T%sqr(i4,i1))**2
    end function

!!!! amplitudes for 2-jet process

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 1st diagram amplitudes
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function amp_2jet_1(T, perm, U) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:), U
    complex(fltknd) :: rslt, xx, zz
    integer :: i1, i2, i3, i4, i5
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    xx = T%sqr(i5,i1)*T%ang(i1,[i1,i5],U,[i2,i4], i2)* T%ang(i2,i4)
    zz = (square(T%Q(i2)%k)+T%ang(i4,i2,i4))*(square(T%Q(i1)%k)+T%ang(i5,i1,i5))
    rslt = -xx/zz
    end function

    function amp_2jet_1m(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    xx= T%sqr(i1,i5)*T%sqr(i1,i5)* T%ang(i2,i4)/T%sqr(i5,i3)*T%Q(i1)%kstr
    yy = T%Q(i2)%kapp*T%ang(i3,i2) + T%ang(i3,i4)*T%sqr(i4,i2)
    zz = (square(T%Q(i1)%k)+T%ang(i5,i1,i5))*(square(T%Q(i2)%k)+T%ang(i4,i2,i4))
    rslt = xx*yy/zz
    end function


    function amp_2jet_1p(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    xx = T%sqr(i5,i1)*T%ang(i2,i4)*T%ang(i2,i4)*T%Q(i2)%kapp/T%ang(i4,i3)
    yy = T%Q(i1)%kstr*T%sqr(i1,i3)+T%ang(i1,i5)*T%sqr(i5,i3)
    zz = (square(T%Q(i1)%k)+T%ang(i5,i1,i5))*(square(T%Q(i2)%k)+T%ang(i4,i2,i4))
    rslt = xx*yy/zz
    end function

!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 2nd diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

    function amp_2jet_2(T, perm,U) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:), U
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    xx = T%sqr(i5,i1)*T%ang(i1,[i1,i5],i2)*T%ang(i2,[i4,i3],U, i4)
    zz = (MZ_sq+T%ang(i4,i3,i4))*(square(T%Q(i1)%k)+T%ang(i5,i1,i5))
    rslt = -xx/zz
    end function

    function amp_2jet_2m(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz, vv
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    vv = T%sqr(i1,i5)*T%ang(i3,i2)/T%sqr(i4, i3)
    xx = T%Q(i1)%kstr*T%sqr(i1,i2) + T%ang(i1,i5)*T%sqr(i5,i2)
    yy = T%ang(i2,i3) + T%ang(i2,i4)*T%sqr(i4,i5)/T%sqr(i3,i5)
    zz = (square(T%Q(i1)%k)+T%ang(i5,i1,i5))
    rslt = vv*xx*yy/zz

    end function

!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 3rd diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

    function amp_2jet_3(T, perm, U) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:), U
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    xx = T%sqr(i5, U,[i3,i5], i1)*T%ang(i1,[i2,i4],i2) *T%ang(i2,i4)
    zz = (square(T%Q(i2)%k)+T%ang(i4,i2,i4))*(MZ_sq+T%ang(i5,i3,i5))
    rslt = -xx/zz
    end function


    function amp_2jet_3p(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,yy,zz, vv
    integer :: i1, i2, i3, i4, i5
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    vv = T%ang(i2,i4)/T%ang(i5, i3)
    xx = T%sqr(i3,i1) + T%ang(i2,i4)*T%sqr(i5,i1)/T%ang(i4,i3)
    yy = T%Q(i2)%kapp*T%ang(i1,i2) + T%ang(i1,i4)*T%sqr(i4,i2)
    zz = (square(T%Q(i2)%k)+T%ang(i4,i2,i4))
    rslt = vv*xx*yy/zz
    end function
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 7th diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

    function amp_2jet_7(T, perm, U) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx, uu,vv, bb, dd, ff, gg
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:), U
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    xx = 1/((square(T%Q(i2)%k) +square(T%Q(i1)%k) +twodot(T%Q(i2)%k, T%Q(i1)%k))*(MZ_sq+T%ang(i4,i3,i4)))
    uu =T%ang(i2,i1)*T%sqr(i1,i2)/2
    vv = T%sqr(i5,i1,[i4,i3], U, i4) - T%sqr(i5,i2,[i4,i3], U, i4)

    bb = T%ang(i1,i2,i1)&
       + T%ang(i2,i1)*T%sqr(i1,i2)*square(T%Q(i2)%k)/T%ang(i2,i1,i2)
    dd = T%sqr(i5,i2)*T%ang(i2,[i4,i3], U, i4)

    ff = T%ang(i2,i1,i2)&
       + T%ang(i2,i1)*T%sqr(i1,i2)*square(T%Q(i1)%k)/T%ang(i1,i2,i1)
    gg = T%sqr(i5,i1)*T%ang(i1,[i4,i3], U, i4)

    rslt = xx*(uu*vv + bb*dd - ff*gg)
    end function

    function amp_2jet_7m(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx, uu,vv, bb, dd, ff, gg
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    xx = T%ang(i3,i4)/((MZ_sq+ T%ang(i4,i3,i4)+ T%ang(i5,i3,i5) &
       + T%sqr(i4,i5)*T%ang(i5,i4))*(MZ_sq+T%ang(i4,i3,i4)))
    uu = T%ang(i2,i1)*T%sqr(i1,i2)/(2*square(T%Q(i2)%k)*square(T%Q(i1)%k))
    vv = (T%sqr(i5,i2,i3)-T%sqr(i5,i1,i3))&
       + T%sqr(i4,i5)/T%sqr(i3,i5)*(T%sqr(i5,i2,i4)-T%sqr(i5,i1,i4))

    bb = -T%ang(i1,i2,i1)/(square(T%Q(i2)%k)*square(T%Q(i1)%k))&
    - T%ang(i2,i1)*T%sqr(i1,i2)/(T%ang(i2,i1,i2)*square(T%Q(i1)%k))
    dd = T%sqr(i5,i2)*T%ang(i2,i3) + T%sqr(i5,i2)*T%ang(i2,i4)*T%sqr(i4,i5)/T%sqr(i3,i5)

    ff = +T%ang(i2,i1,i2)/(square(T%Q(i2)%k)*square(T%Q(i1)%k))&
    + T%ang(i2,i1)*T%sqr(i1,i2)/(T%ang(i1,i2,i1)*square(T%Q(i2)%k))
    gg = T%sqr(i5,i1)*T%ang(i1,i3) + T%sqr(i5,i1)*T%ang(i1,i4)*T%sqr(i4,i5)/T%sqr(i3,i5)

    rslt = xx*(uu*vv + bb*dd + ff*gg)*(T%Q(i1)%kapp*T%Q(i1)%kstr*T%Q(i2)%kapp*T%Q(i2)%kstr)
    end function


!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 8th diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


    function amp_2jet_8(T, perm, U) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx, uu,vv, bb, dd, ff, gg
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:), U
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    xx = 1/((square(T%Q(i2)%k) +square(T%Q(i1)%k) +twodot(T%Q(i2)%k, T%Q(i1)%k))*(MZ_sq+T%ang(i5,i3,i5)))
    uu = T%ang(i2,i1)*T%sqr(i1,i2)/2
    vv = T%sqr(i5,U,[i3,i5], i2, i4) - T%sqr(i5,U,[i3,i5], i1, i4)

    bb = T%ang(i1,i2,i1)&
       + T%ang(i2,i1)*T%sqr(i1,i2)*square(T%Q(i2)%k)/T%ang(i2,i1,i2)
    dd = T%sqr(i5,U,[i3,i5], i2)*T%ang(i2, i4)

    ff = T%ang(i2,i1,i2)&
       + T%ang(i2,i1)*T%sqr(i1,i2)*square(T%Q(i1)%k)/T%ang(i1,i2,i1)
    gg = T%sqr(i5,U,[i3,i5], i1)*T%ang(i1, i4)

    rslt =xx*(uu*vv - bb*dd + ff*gg)
    end function

  function amp_2jet_8p(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx, uu,vv, bb, dd, ff, gg
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i2=perm(2) ;i5=4 ;i3=5; i4=3; i1=perm(1)
    !
    rslt = 0
    xx = T%sqr(i5,i3)/((MZ_sq+ T%ang(i4,i3,i4)+ T%ang(i5,i3,i5) &
       + T%sqr(i4,i5)*T%ang(i5,i4))*(MZ_sq+T%ang(i5,i3,i5)))
    uu = T%ang(i2,i1)*T%sqr(i1,i2)/(2*square(T%Q(i2)%k)*square(T%Q(i1)%k))
    vv = (T%sqr(i3,i2,i4)-T%sqr(i3,i1,i4))&
       + T%ang(i4,i5)/T%ang(i4,i3)*(T%sqr(i5,i2,i4)-T%sqr(i5,i1,i4))
    bb = -T%ang(i1,i2,i1)/(square(T%Q(i2)%k)*square(T%Q(i1)%k))&
    - T%ang(i2,i1)*T%sqr(i1,i2)/(T%ang(i2,i1,i2)*square(T%Q(i1)%k))
    dd = T%sqr(i3,i2)*T%ang(i2,i4) + T%sqr(i5,i2)*T%ang(i2,i4)*T%ang(i4,i5)/T%ang(i4,i3)

    ff = T%ang(i2,i1,i2)/(square(T%Q(i2)%k)*square(T%Q(i1)%k))&
    + T%ang(i2,i1)*T%sqr(i1,i2)/(T%ang(i1,i2,i1)*square(T%Q(i2)%k))
    gg = T%sqr(i3,i1)*T%ang(i1,i4) + T%sqr(i5,i1)*T%ang(i1,i4)*T%ang(i4,i5)/T%ang(i4,i3)

    rslt = xx*(uu*vv + bb*dd + ff*gg)*(T%Q(i2)%kapp*T%Q(i2)%kstr*T%Q(i1)%kapp*T%Q(i1)%kstr)

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


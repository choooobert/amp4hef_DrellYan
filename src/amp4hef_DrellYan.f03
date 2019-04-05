module amp4hef_DrellYan
  use amp4hef_io
  use amp4hef_aux
  use amp4hef_qomentum
  implicit none
  private
  public :: fill_matrices_DrellYan,matrix_element_DrellYan ,amplitude_DrellYan ,all_amplitudes_DrellYan

  real, parameter :: sqrt_2 = 1.41421356237_fltknd
  real, parameter :: MZ = 0.8315287819239998E+04
  real, parameter :: cV = 1._fltknd
  real, parameter :: cA = 1.26_fltknd

  integer,parameter :: gluon=0 ,quark=1 ,antiq=-1, Zboson=2
	 integer,parameter :: helTable_DrellYan(3,6)=reshape(&
	 [ -1,1,-1,	 -1, 1, 0,	-1, 1, 1,&
	   1,-1,-1,   1,-1, 0,	 1,-1, 1 &
	 ], [3,6])

  integer,allocatable,save :: mtx_4_sqr(:,:)

contains

  function matrix_element_DrellYan(Tin) result(rslt)
  class(qomentum_list_type),intent(in) :: Tin
  real(fltknd) :: rslt
  complex(fltknd) :: amp(12,4)
  integer :: ii,NhelSum,Nminus2, NhelConf, Nperm, jj
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )

    NhelConf = 6
	NhelSum = 3
	NPerm = 2
!
  rslt = 0
	do ii=1, (NhelConf/2)
		do jj=1, Nperm
      amp(ii, jj) = amplitude_DrellYan(Tin ,helTable_DrellYan(:,ii), perTable(1:NPerm,jj))
      amp(NhelConf-ii+1,jj) =conjg(amp(ii,jj))
      rslt = (cV+cA)*amp(ii, jj)+ (cV-cA)*amp(NhelConf-ii+1,jj)
!      write(*,'(2e16.8,99i3)') amp(ii),helTable_DrellYan(1:NhelSum,ii) !DEBUG
    enddo
!    write(*,*) !DEBUG
  enddo
  rslt = rslt*2
  end associate
  end function 


  subroutine all_amplitudes_DrellYan(Tin ,NhelConf ,Nperm ,amplitude ,factor )
  class(qomentum_list_type),intent(in) :: Tin
  integer,intent(out) :: NhelConf, Nperm
  complex(fltknd),intent(out) :: amplitude(:,:),factor
  integer :: ii,jj,NhelSum
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )
	NhelSum = 3
	NhelConf= 6
	Nperm = 2
  do ii=1,(NhelConf/2)
  do jj=1,Nperm
		amplitude(ii,jj) =(cV+cA)*amplitude_DrellYan(Tin ,helTable_DrellYan(1:NhelSum,ii) ,perTable(1:Nperm,jj) )
        amplitude(NhelConf-ii+1,jj) =(cV-cA)*conjg(amplitude(ii,jj))

  enddo
	enddo
	
  factor = 2 ! only half of the helicity configurations is returned
  end associate
  end subroutine


 function amplitude_DrellYan(Tin ,helicity ,perm ) result(rslt)
! The first entries of helicity refer to the on-shell gluons. The
! next helicity to the anti-quark. Off-shell (anti)-quarks have helicity,
! and the quark always gets the opposite helicity of the anti-quark.
! If 1(2) gluon(s) is(are) off-shell, it must be the first (2) gluon(s).
  class(qomentum_list_type),intent(in) :: Tin
  integer,intent(in) :: helicity(:), perm(:)
  complex(fltknd) :: rslt
 ! type(qomentum_list_type) :: T
  integer :: hel(-1:NsizeProc)
  integer ::  i2, i3, i4
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff ,offshell=>Tin%offshell )
  hel(3:5) = helicity
	!
	 i2=1 ;i3=3; i4=2
	!
!  T%Ntot = Ntot

!what is it for?
!  T%Q(Ntot-1)   = Tin%Q(Tin%flavor(    1        ,antiq))
!  T%Q(Ntot)     = Tin%Q(Tin%flavor(    1        ,quark))
	rslt = 0
	if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.-1) then
	rslt = 24*( amp_01(Tin, perm) + amp_07(Tin, perm) + amp_13(Tin, perm) &
	          + amp_19(Tin, perm) + amp_25(Tin, perm) + amp_31(Tin, perm) )
	else if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.0) then
	rslt = 24*( amp_02(Tin, perm) + amp_08(Tin, perm) + amp_14(Tin, perm) &
		      + amp_20(Tin, perm) + amp_26(Tin, perm) + amp_32(Tin, perm) )
	else if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.1) then
	rslt = 24*( amp_03(Tin, perm) + amp_09(Tin, perm) + amp_15(Tin, perm) &
	          + amp_21(Tin, perm) + amp_27(Tin, perm) + amp_33(Tin, perm) )
	end if






  end associate
end function

!!!! amplitude case part,++ not relevant now

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 1st diagram amplitudes
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	function amp_01(T, perm) result(rslt)
	type(qomentum_list_type),intent(in) :: T
	complex(fltknd) :: rslt,xx,yy,zz
	integer :: i1, i2, i3, i4, i5
	integer, intent(in) :: perm(:)

	!
	i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
	!
	rslt = 0
	call T%set_direction(i3, i4)
	xx= T%sqr(i4,i5)*T%sqr(i4,i5)* T%ang(i1,i2)
	yy = T%ang(i3,i1) - T%ang(i3,i2)*T%sqr(i2,i1)/ T%Q(i1)%kapp
	zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
	if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
	end function
	
	function amp_02(T, perm) result(rslt)
	type(qomentum_list_type),intent(in) :: T
	integer, intent(in) :: perm(:)
	complex(fltknd) :: rslt,tt,uu,vv,xx,ww,yy,zz
	integer :: i1, i2, i3, i4, i5
	!
	i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
	!
	rslt = 0
	call T%set_direction(i3, i4)
	tt = T%sqr(i4,i5)*T%ang(i1,i2)
	uu = (-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i1)%kstr * T%Q(i5)%kapp
	vv = T%sqr(i5,i4)
	ww = T%ang(i4,i3) *T%sqr(i3,i4)
	xx = T%ang(i4,i1) - T%ang(i4,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
	yy = T%sqr(i5,i3) - T%ang(i5,i4)*T%sqr(i4,i3)/T%Q(i5)%kstr
	zz = T%ang(i3,i1) - T%ang(i3,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
	if (tt.ne.0.and.uu.ne.0.and.ww.ne.0) rslt = 4*tt/uu*(MZ*vv/ww*xx-yy*zz/MZ)
	end function
	
	function amp_03(T, perm) result(rslt)
	type(qomentum_list_type),intent(in) :: T
	integer, intent(in) :: perm(:)
	complex(fltknd) :: rslt,xx,yy,zz
	integer :: i1, i2, i3, i4, i5
	!
	i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
	!
	rslt = 0
	call T%set_direction(i3, i2)

	xx = T%ang(i1,i2) * T%ang(i1,i2) * T%sqr(i4,i5) 
	yy = T%sqr(i5,i3) - T%ang(i5,i4)*T%sqr(i4,i3)/T%Q(i5)%kstr
	zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
	if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
	end function
	


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2nd diagram amplitudes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function amp_07(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)

    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i2)
    xx= T%sqr(i4,i5)*T%ang(i2,i3)* T%ang(i1,i3)
    yy = T%sqr(i5,i1) - T%ang(i5,i4)*T%sqr(i4,i1)/ T%Q(i5)%kstr
    zz =(MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i1)%kapp*T%Q(i1)%kstr*T%Q(i5)%kapp

    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
    end function

    function amp_08(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,ww,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    call T%set_direction(i3, i2)
    rslt = 0
    call T%set_direction(i3, i4)
    ww = T%ang(i1,i2)*T%sqr(i4,i5)
    xx = T%sqr(i5,i1) - T%ang(i5,i4)*T%sqr(i4,i1)/T%Q(i5)%kstr
    yy = MZ + T%sqr(i2,i3)*T%ang(i3,i2)/MZ
    zz = (MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i5)%kapp*T%Q(i1)%kapp*T%Q(i1)%kstr
    if (ww.ne.0.and.xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*ww*xx*yy/zz
    end function

    function amp_09(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i1)

    xx = T%sqr(i5,i4) * T%ang(i3,i2) * T%ang(i1,i3)
    yy = T%sqr(i5,i1) - T%ang(i5,i4)*T%sqr(i4,i2)/T%Q(i5)%kstr
    zz = (MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i5)%kapp*T%Q(i1)%kstr*T%Q(i1)%kapp
    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 3rd diagram amplitudes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function amp_13(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i5)
    xx= T%sqr(i4,i5)*T%sqr(i4,i5)*T%ang(i1,i2)* T%ang(i3,i4)
    yy = T%ang(i5,i1) - T%ang(i5,i2)*T%sqr(i2,i1)/ T%Q(i1)%kapp
    zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr - T%ang(i2,i1,i2))*(MZ*MZ+T%ang(i4,i3,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp*T%Q(i5)%kstr
    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
    end function

    function amp_14(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,ww,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    call T%set_direction(i3, i4)
    rslt = 0
    call T%set_direction(i3, i4)
    ww = T%ang(i2,i1)*T%sqr(i4,i5)
    xx = T%ang(i5,i1) - T%ang(i5,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
    yy = MZ + T%sqr(i4,i3)*T%ang(i3,i4)/MZ
    zz = (MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i5)%kapp*T%Q(i5)%kstr*T%Q(i1)%kstr
    if (ww.ne.0.and.xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*ww*xx*yy/zz
    end function

    function amp_15(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i4)

    xx = T%sqr(i4,i3) * T%sqr(i3,i5) * T%ang(i1,i2)
    yy = T%ang(i5,i1) - T%ang(i5,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
    zz = (MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i5)%kapp*T%Q(i5)%kstr*T%Q(i1)%kstr
    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 4th diagram amplitudes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function amp_19(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i4)
    xx= T%sqr(i4,i1)*T%sqr(i4,i1)* T%ang(i2,i5)
    yy = T%ang(i3,i5) - T%ang(i3,i2)*T%sqr(i2,i5)/ T%Q(i5)%kapp
    zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i4,i1,i4))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i2,i5,i2))*T%Q(i5)%kstr*T%Q(i1)%kapp

    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
    end function

    function amp_20(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,tt,uu,vv,xx,ww,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    call T%set_direction(i3, i4)
    rslt = 0
    call T%set_direction(i3, i4)
    tt = T%sqr(i4,i1)*T%ang(i5,i2)
    uu = (-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i4,i1,i4))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i2,i5,i2))*T%Q(i5)%kstr * T%Q(i1)%kapp
    vv = T%sqr(i1,i4)
    ww = T%ang(i4,i3) *T%sqr(i3,i4)

    xx = T%ang(i4,i5) - T%ang(i4,i2)*T%sqr(i2,i5)/T%Q(i5)%kapp
    yy = T%sqr(i1,i3) - T%ang(i1,i4)*T%sqr(i4,i3)/T%Q(i1)%kstr
    zz = T%ang(i3,i5) - T%ang(i3,i2)*T%sqr(i2,i5)/T%Q(i5)%kapp
    if (tt.ne.0.and.uu.ne.0.and.ww.ne.0) rslt = 4*tt/uu*(MZ*vv*xx/ww-yy*zz/MZ)
    end function
    function amp_21(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i2)
    xx = T%ang(i2,i5) * T%ang(i2,i5) * T%sqr(i1,i4)
    yy = T%sqr(i1,i3) - T%ang(i1,i4)*T%sqr(i4,i3)/T%Q(i1)%kstr
    zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i4,i1,i4))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i2,i5,i2))*T%Q(i5)%kstr*T%Q(i1)%kapp
    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 5th diagram amplitudes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function amp_25(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)

    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i2)
    xx= T%sqr(i4,i1)*T%ang(i3,i2)* T%ang(i5,i3)
    yy = T%sqr(i1,i5) - T%ang(i1,i4)*T%sqr(i4,i5)/ T%Q(i1)%kstr
    zz =(MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i4,i1,i4))*T%Q(i5)%kapp*T%Q(i5)%kstr*T%Q(i1)%kapp

    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
    end function

    function amp_26(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,ww,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i2)
    ww = T%ang(i1,i4)*T%sqr(i5,i2)
    xx = T%sqr(i1,i5) - T%ang(i1,i4)*T%sqr(i4,i5)/T%Q(i1)%kstr
    yy = MZ + T%sqr(i2,i3)*T%ang(i3,i2)/MZ
    zz = (MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i4,i1,i4))*T%Q(i1)%kapp*T%Q(i5)%kapp*T%Q(i5)%kstr
    if (ww.ne.0.and.xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*ww*xx*yy/zz
    end function

    function amp_27(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i2)

    xx = T%sqr(i4,i1) * T%ang(i3,i2) * T%ang(i5,i3)
    yy = T%sqr(i1,i5) - T%ang(i1,i4)*T%sqr(i4,i5)/T%Q(i1)%kstr
    zz = (MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i4,i1,i4))*T%Q(i1)%kapp*T%Q(i5)%kstr*T%Q(i5)%kapp
    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 6th diagram amplitudes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function amp_31(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    integer, intent(in) :: perm(:)
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i1)
    xx= T%sqr(i4,i1)*T%sqr(i4,i1)*T%ang(i5,i2)* T%ang(i3,i4)
    yy = T%ang(i5,i1) - T%ang(i1,i2)*T%sqr(i2,i5)/ T%Q(i5)%kapp
    zz = (-T%Q(i5)%kapp*T%Q(i5)%kstr - T%ang(i2,i5,i2))*(MZ*MZ+T%ang(i4,i3,i4))*T%Q(i5)%kstr*T%Q(i1)%kapp*T%Q(i1)%kstr
    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
    end function

    function amp_32(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,ww,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i4)
    ww = T%ang(i5,i2)*T%sqr(i1,i4)
    xx = T%ang(i1,i5) - T%ang(i1,i2)*T%sqr(i2,i5)/T%Q(i5)%kapp
    yy = MZ + T%sqr(i3,i4)*T%ang(i4,i3)/MZ
    zz = (MZ*MZ + T%ang(i4, i3, i4))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i2,i5,i2))*T%Q(i1)%kapp*T%Q(i1)%kstr*T%Q(i5)%kstr
    if (ww.ne.0.and.xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*ww*xx*yy/zz
    end function

    function amp_33(T, perm) result(rslt)
    type(qomentum_list_type),intent(in) :: T
    integer, intent(in) :: perm(:)
    complex(fltknd) :: rslt,xx,yy,zz
    integer :: i1, i2, i3, i4, i5
    !
    i1=perm(2) ;i2=3 ;i3=4; i4=5; i5=perm(1)
    !
    rslt = 0
    call T%set_direction(i3, i4)

    xx = T%sqr(i4,i3) * T%sqr(i3,i1) * T%ang(i5,i2)
    yy = T%ang(i1,i5) - T%ang(i1,i2)*T%sqr(i2,i5)/T%Q(i5)%kapp
    zz = (MZ*MZ + T%ang(i4, i3, i4))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i2,i5,i2))*T%Q(i1)%kapp*T%Q(i1)%kstr*T%Q(i5)%kstr
    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
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


end module




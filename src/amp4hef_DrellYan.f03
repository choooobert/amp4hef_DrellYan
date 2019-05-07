module amp4hef_DrellYan
  use amp4hef_io
  use amp4hef_aux
  use amp4hef_qomentum
  implicit none
  private
  public :: fill_matrices_DrellYan,matrix_element_DrellYan ,amplitude_DrellYan ,all_amplitudes_DrellYan

  real, parameter :: sqrt_2 = 1.41421356237_fltknd
  real, parameter :: MZ = 91.1882
!  real, parameter :: MZ = 0
  real(fltknd), parameter:: MZ_sq = 8315.2878192399512
!  real(fltknd), parameter:: MZ_sq = 0

  real, parameter :: cV = 0.203666_fltknd
  real, parameter :: cA = 0.5_fltknd

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
    if(Ntot.eq.5) then
	    NPerm = 2
    else if(Ntot.eq.4) then
        NPerm = 1
    end if
  rslt = 0
    do ii=1, (NhelConf/2)
        do jj=1, Nperm
      amp(ii, jj) = amplitude_DrellYan(Tin ,helTable_DrellYan(:,ii), perTable(1:NPerm,jj))
      rslt = rslt + amp(ii, jj)*conjg(amp(ii,jj))
!      write(*,'(2e16.8,99i3)') amp(ii),helTable_DrellYan(1:NhelSum,ii) !DEBUG
    enddo
!    write(*,*) !DEBUG
  enddo
  do ii=1, Noff
      rslt = rslt*(Tin%Q(ii)%kapp*Tin%Q(ii)%kstr)
  enddo
  rslt= 2*rslt
  !strong factor
  rslt = 8*rslt
  !weak factor
  rslt = 2*(cV*cV+cA*cA)*rslt
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
!		amplitude(ii,jj) =(cV+cA)*amplitude_DrellYan(Tin ,helTable_DrellYan(1:NhelSum,ii) ,perTable(1:Nperm,jj) )
!        amplitude(NhelConf-ii+1,jj) =(cV-cA)*conjg(amplitude(ii,jj))
  enddo
  enddo
	
  factor = 2 ! only half of the helicity configurations is returned
  end associate
  end subroutine


 function amplitude_DrellYan(Tin ,helicity ,perm ) result(rslt)
  class(qomentum_list_type),intent(in) :: Tin
  integer,intent(in) :: helicity(:), perm(:)
  complex(fltknd) :: rslt
  integer ::  i2, i3, i4
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff ,offshell=>Tin%offshell )
	!
	rslt = 0
	i2=1 ;i3=3; i4=2 !indexing only for helicity configuration not for qomentum class
    if(Ntot.eq.5) then
!	if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.-1) then
!	rslt =  amp_01(Tin, perm) + amp_07(Tin, perm) &
!	       + amp_19(Tin, perm) + amp_25(Tin, perm) &
!	       + amp_37(Tin, perm))
!
!	else if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.0) then
!	rslt = ( amp_02(Tin, perm) + amp_08(Tin, perm) + amp_14(Tin, perm) &
!		      + amp_20(Tin, perm) + amp_26(Tin, perm) + amp_32(Tin, perm) &
!		      + amp_38(Tin, perm) + amp_44(Tin, perm))
!
!	else if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.1) then
!	rslt = ( amp_03(Tin, perm) + amp_15(Tin, perm) &
!	       + amp_21(Tin, perm) + amp_33(Tin, perm) &
!	       + amp_45(Tin, perm))
!	end if
    else if(Ntot.eq.4) then
        if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.-1) then
!        write(*,*) "helicity", helicity(i2), helicity(i4), helicity(i3)
!        write(*,*)  "amp 101"
        rslt = amp_101(Tin)
        else if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.0) then
!        write(*,*) "helicity", helicity(i2), helicity(i4), helicity(i3)
!        write(*,*)  "amp 102"
        rslt = amp_102(Tin)
        else if (helicity(i2).eq.-1.and.helicity(i4).eq.1.and.helicity(i3).eq.1) then
!        write(*,*) "helicity", helicity(i2), helicity(i4), helicity(i3)
!        write(*,*)  "amp 103"
        rslt = amp_103(Tin)
        end if
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
    complex(fltknd) :: rslt, vv ,xx,yy, zz, uu, ww
    integer :: i1, i2, i3, i4
    !
    i1=1 ;i2=2 ;i3=4; i4=3
    !
    rslt = 0
!    write(*,*) "102"
!    call T%set_direction(i3,i1)
!    vv = MZ/(T%ang(i3,i1)*T%sqr(i1,i3))
!    xx = (T%sqr(i4,i3)*T%ang(i3,i4)/MZ + vv*T%ang(i1,i2)*T%sqr(i2,i1))/(MZ_sq+T%ang(i4,i3,i4))
!    yy = (T%sqr(i2,i3)*T%ang(i3,i2)/MZ + vv*T%ang(i1,i4)*T%sqr(i4,i1))/(MZ_sq+T%ang(i2,i3,i2))
!    zz = T%sqr(i4,i1)*T%ang(i1,i2)/(T%Q(i1)%kapp*T%Q(i1)%kstr)
!    rslt = sqrt_2*zz*(-xx+yy)


    call T%set_direction(i3,i2)
    uu = T%sqr(i4,i2)*T%ang(i2,i1)*T%ang(i2,i1)/(T%ang(i3,i2)*T%sqr(i2,i3))
    ww = 1/(MZ*T%Q(i1)%kapp*T%Q(i1)%kstr)
    xx = T%sqr(i4,i3)*T%ang(i1,i2)*(T%ang(i3,i1)*T%Q(i1)%kapp+T%ang(i3,i2)*T%sqr(i2,i1)) &
       /(MZ_sq+T%ang(i4,i3,i4))
    yy = T%sqr(i4,i1)*T%ang(i3,i2)*(T%sqr(i1,i3)*T%Q(i1)%kstr+T%ang(i1,i4)*T%sqr(i4,i3)) &
       /(MZ_sq+T%ang(i2,i3,i2))
    rslt = sqrt_2*(ww*uu+ww*(xx-yy))
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
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	! 1st diagram amplitudes
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	function amp_01(T, perm) result(rslt)
!	type(qomentum_list_type),intent(in) :: T
!	complex(fltknd) :: rslt,xx,yy,zz
!	integer :: i1, i2, i3, i4, i5
!	integer, intent(in) :: perm(:)
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!	rslt = 0
!	xx= T%sqr(i4,i5)*T%sqr(i4,i5)* T%ang(i1,i2)
!	yy = T%ang(i3,i1) + T%ang(i3,i2)*T%sqr(i2,i1)/ T%Q(i1)%kapp
!	zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr+T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
!	rslt = 2*sqrt_2*xx*yy/zz
!	end function
!
!	function amp_02(T, perm) result(rslt)
!	type(qomentum_list_type),intent(in) :: T
!	integer, intent(in) :: perm(:)
!	complex(fltknd) :: rslt,tt,uu,vv,xx,ww,yy,zz
!	integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!	rslt = 0
!	call T%set_direction(i3, i4)
!	tt = T%sqr(i4,i5)*T%ang(i1,i2)
!	uu = (-T%Q(i1)%kapp*T%Q(i1)%kstr+T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
!	ww = T%sqr(i5,i3) + T%ang(i5,i4)*T%sqr(i4,i3)/T%Q(i5)%kstr
!	xx = T%ang(i3,i1) + T%ang(i3,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
!	yy = T%sqr(i5,i4)/(T%ang(i3,i4)*T%sqr(i4,i3))
!	zz = T%ang(i4,i1) + T%ang(i4,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
!	rslt = -2*tt/uu*(-ww*xx/MZ+MZ*yy*zz)
!	end function
!
!	function amp_03(T, perm) result(rslt)
!	type(qomentum_list_type),intent(in) :: T
!	integer, intent(in) :: perm(:)
!	complex(fltknd) :: rslt,xx,yy,zz
!	integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!	rslt = 0
!
!	xx = T%sqr(i4,i5)*T%ang(i1,i2)*T%ang(i1,i2)/T%ang(i2,i3)
!	yy = T%sqr(i5,i3)+t%ang(i5,i4)*T%sqr(i4,i3)/T%Q(i5)%kstr
!	zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr+T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
!	rslt = 2*sqrt_2*xx*yy/zz
!	end function
!
!
!
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 2nd diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    function amp_07(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,xx,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    integer, intent(in) :: perm(:)
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    vv = T%sqr(i4,i5)*T%ang(i3,i2)
!    xx = T%sqr(i5,i1) + T%ang(i5,i4)*T%sqr(i4,i1)/ T%Q(i5)%kstr
!    yy = T%ang(i1,i3) + T%ang(i1,i2)*T%sqr(i2,i4)/T%sqr(i3,i4)
!    zz =(MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i4,i5,i4))*T%Q(i1)%kapp*T%Q(i1)%kstr*T%Q(i5)%kapp
!    rslt = 2*sqrt_2*vv*xx*yy/zz
!    end function
!
!    function amp_08(T, perm) result(rslt)
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
!    function amp_14(T, perm) result(rslt)
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
!    function amp_15(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    vv = T%sqr(i3,i4)* T%ang(i1,i2)
!    xx = T%sqr(i3,i5)+T%ang(i2,i4)*T%sqr(i4,i5)/T%ang(i2,i3)
!    yy = T%ang(i5,i1) + T%ang(i5,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
!    zz = (MZ*MZ + T%ang(i4, i3, i4))*(-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i2,i1,i2))*T%Q(i5)%kapp*T%Q(i5)%kstr*T%Q(i1)%kstr
!    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*vv*xx*yy/zz
!    end function
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 4th diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    function amp_19(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,xx,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    integer, intent(in) :: perm(:)
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    xx = T%sqr(i4,i1)*T%sqr(i4,i1)* T%ang(i5,i2)/T%sqr(i4,i3)
!    yy = T%ang(i3,i5)+T%ang(i3,i2)*T%sqr(i2,i5)/T%Q(i5)%kapp
!    zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr+T%ang(i4,i1,i4))*(-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i2,i5,i2))*T%Q(i5)%kstr*T%Q(i1)%kapp
!
!    rslt = 2*sqrt_2*xx*yy/zz
!    end function
!
!    function amp_20(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,tt,uu,vv,xx,ww,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    call T%set_direction(i3, i4)
!    tt = T%sqr(i4,i1)*T%ang(i5,i2)
!    uu = (-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i2,i5,i2))*(-T%Q(i1)%kapp*T%Q(i1)%kstr+T%ang(i4,i1,i4))*T%Q(i5)%kstr*T%Q(i1)%kapp
!    ww = T%sqr(i1,i3) + T%ang(i1,i4)*T%sqr(i4,i3)/T%Q(i1)%kstr
!    xx = T%ang(i3,i5) + T%ang(i3,i2)*T%sqr(i2,i5)/T%Q(i5)%kapp
!    yy = T%sqr(i1,i4)/(T%ang(i3,i4)*T%sqr(i4,i3))
!    zz = T%ang(i4,i5) + T%ang(i4,i2)*T%sqr(i2,i5)/T%Q(i5)%kapp
!    rslt = -2*tt/uu*(-ww*xx/MZ+MZ*yy*zz)
!    end function
!
!    function amp_21(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    xx = T%sqr(i4,i1)*T%ang(i5,i2)*T%ang(i5,i2)/T%ang(i2,i3)
!    yy = T%sqr(i1,i3) + T%ang(i1,i4)*T%sqr(i4,i3)/T%Q(i1)%kstr
!    zz = (-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i2,i5,i2))*(-T%Q(i1)%kapp*T%Q(i1)%kstr+T%ang(i1,i4,i1))*T%Q(i5)%kstr*T%Q(i1)%kapp
!    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
!    end function
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 5th diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    function amp_25(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,xx,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    integer, intent(in) :: perm(:)
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    vv = T%sqr(i4,i1)*T%ang(i3,i2)
!    xx = T%sqr(i1,i5)+T%ang(i1,i4)*T%sqr(i4,i5)/T%Q(i1)%kstr
!    yy = T%ang(i5,i3)+T%ang(i5,i2)*T%sqr(i2,i4)/T%sqr(i3,i4)
!    zz = (MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i1)%kapp*T%Q(i1)%kstr+T%ang(i4,i1,i4))*T%Q(i5)%kapp*T%Q(i5)%kstr*T%Q(i1)%kapp
!
!    if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*vv*xx*yy/zz
!    end function
!
!    function amp_26(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,ww,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    call T%set_direction(i3, i2)
!    ww = T%ang(i1,i4)*T%sqr(i5,i2)
!    xx = T%sqr(i1,i5) - T%ang(i1,i4)*T%sqr(i4,i5)/T%Q(i1)%kstr
!    yy = MZ + T%sqr(i2,i3)*T%ang(i3,i2)/MZ
!    zz = (MZ*MZ + T%ang(i2, i3, i2))*(-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i4,i1,i4))*T%Q(i1)%kapp*T%Q(i5)%kapp*T%Q(i5)%kstr
!    if (ww.ne.0.and.xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*ww*xx*yy/zz
!    end function
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 6th diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    function amp_32(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,ww,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    call T%set_direction(i3, i4)
!    ww = T%ang(i5,i2)*T%sqr(i1,i4)
!    xx = T%ang(i1,i5) - T%ang(i1,i2)*T%sqr(i2,i5)/T%Q(i5)%kapp
!    yy = MZ + T%sqr(i3,i4)*T%ang(i4,i3)/MZ
!    zz = (MZ*MZ + T%ang(i4, i3, i4))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i2,i5,i2))*T%Q(i1)%kapp*T%Q(i1)%kstr*T%Q(i5)%kstr
!    if (ww.ne.0.and.xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*ww*xx*yy/zz
!    end function
!
!    function amp_33(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    integer, intent(in) :: perm(:)
!    complex(fltknd) :: rslt,xx,yy,zz
!    integer :: i1, i2, i3, i4, i5
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    vv = T%sqr(i3,i4)*T%ang(i5,i2)
!    xx = T%sqr(i3,i1) + T%ang(i2,i4)*T%sqr(i4,i1)/T%ang(i2,i3)
!    yy = T%ang(i1,i5) + T%ang(i1,i2)*T%sqr(i2,i5)/T%Q(i5)%kapp
!    zz = (MZ*MZ + T%ang(i4, i3, i4))*(-T%Q(i5)%kapp*T%Q(i5)%kstr+T%ang(i2,i5,i2))*T%Q(i1)%kapp*T%Q(i1)%kstr*T%Q(i5)%kstr
!    rslt = 2*sqrt_2*vv*xx*yy/zz
!    end function
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! 7th diagram amplitudes
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    function amp_37(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,xx, zz
!    integer :: i1, i2, i3, i4, i5
!    integer, intent(in) :: perm(:)
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    call T%set_direction(i3, i2)
!    xx= T%ang(i1,i5)*T%sqr(i5,i1)*T%sqr(i4,i5)*T%ang(i5,i3)*T%ang(i2,i3)
!    zz = T%Q(i1)%kapp*T%Q(i1)%kstr*(MZ*MZ+T%ang(i2,i3,i2))*T%ang(i5,i1,i5) &
!       *(MZ*MZ+T%ang(i2,i4)*T%sqr(i4,i2)+T%ang(i2,i3,i2)+T%ang(i4,i3,i4))
!    if (xx.ne.0.and.zz.ne.0) rslt =-2*xx/zz
!    end function
!
!    function amp_38(T, perm) result(rslt)
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
!    function amp_44(T, perm) result(rslt)
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
!    function amp_45(T, perm) result(rslt)
!    type(qomentum_list_type),intent(in) :: T
!    complex(fltknd) :: rslt,xx, zz
!    integer :: i1, i2, i3, i4, i5
!    integer, intent(in) :: perm(:)
!    !
!    i1=2 ;i2=4 ;i3=5; i4=3; i5=1
!    !
!    rslt = 0
!    call T%set_direction(i3, i4)
!    xx= T%ang(i1,i5)*T%sqr(i5,i1)*T%sqr(i4,i3)*T%sqr(i3,i5)*T%ang(i5,i2)
!    zz = T%Q(i1)%kapp*T%Q(i1)%kstr*(MZ*MZ+T%ang(i4,i3,i4))*T%ang(i5,i1,i5) &
!       *(MZ*MZ+T%ang(i2,i4)*T%sqr(i4,i2)+T%ang(i2,i3,i2)+T%ang(i4,i3,i4))
!    if (xx.ne.0.and.zz.ne.0) rslt =2*xx/zz
!    end function

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




module amp4hef_DrellYan
  use amp4hef_io
  use amp4hef_aux
  use amp4hef_qomentum
! TODO: fix problem with helicity table
! TODO: make outcome of DrellYan amplitudes non-zero
! TODO: check if set_direction work correctly
  implicit none
  private
  public :: fill_matrices_DrellYan,matrix_element_DrellYan ,amplitude_DrellYan ,all_amplitudes_DrellYan

  real, parameter :: sqrt_2 = 1.41421356237_fltknd
  integer,parameter :: gluon=0 ,quark=1 ,antiq=-1, Zboson=2
	 integer,parameter :: helTable_DrellYan(3,12)=reshape(&
	 [-1,-1,-1,	 -1,-1, 1,	-1, 0,-1,	-1, 0, 1,	-1, 1,-1,	 -1, 1,1,&
	   1,-1,-1,   1,-1, 1,	 1, 0,-1,	 1, 0, 1,	 1, 1,-1,	  1, 1, 1&
	 ], [3,12])

  integer,allocatable,save :: mtx_4_sqr(:,:)

contains

  function matrix_element_DrellYan(Tin) result(rslt)
  class(qomentum_list_type),intent(in) :: Tin
  real(fltknd) :: rslt
  complex(fltknd) :: amp(12,4)
  integer :: ii,NhelSum,Nminus2, NhelConf, Nperm, jj
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )
    NhelConf = 12
	NhelSum = 3
	NPerm = 2
!
  rslt = 0
	do jj=1, Nperm
		do ii=1, NhelConf
      amp(ii, jj) = amplitude_DrellYan(Tin ,helTable_DrellYan(:,ii), perTable(1:NPerm,jj))
      rslt = rslt + amp(ii, jj)
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
	NhelConf =12
	Nperm = 2
  do ii=1,NhelConf
  do jj=1,Nperm
		amplitude(ii,jj) = amplitude_DrellYan(Tin ,helTable_DrellYan(1:NhelSum,ii) ,perTable(1:Nperm,jj) )
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
  type(qomentum_list_type) :: T
  integer :: hel(-1:NsizeProc)
  integer :: i1, i2, i3, i4, i5
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff ,offshell=>Tin%offshell )

  hel(3:5) = helicity
	!
	i1=2 ;i2=3 ;i3=4; i4=5; i5=1
	!
  T%Ntot = Ntot

!what is it for?
!  T%Q(Ntot-1)   = Tin%Q(Tin%flavor(    1        ,antiq))
!  T%Q(Ntot)     = Tin%Q(Tin%flavor(    1        ,quark))
    write (*,*) "helicity : ", helicity(i3)
	rslt = 0
	if (helicity(i2).eq.-1.and.helicity(i3).eq.-1.and.helicity(i4).eq.1) then 
	rslt = amp_1(T) !+ amp_7(T) + amp_13(T)
	else if (helicity(i2).eq.-1.and.helicity(i3).eq.0.and.helicity(i4).eq.1) then 
	rslt = amp_2(T) !+ amp_8(T) + amp_14(T)
	else if (helicity(i2).eq.-1.and.helicity(i3).eq.1.and.helicity(i4).eq.1) then 
	rslt = amp_3(T) !+ amp_9(T) + amp_15(T)
	
	else if (helicity(i2).eq.1.and.helicity(i3).eq.-1.and.helicity(i4).eq.-1) then 
	rslt = amp_1(T) !+ amp_7(T) + amp_13(T)
	rslt = conjg(rslt)
	else if (helicity(i2).eq.1.and.helicity(i3).eq.0.and.helicity(i4).eq.-1) then 
	rslt = amp_2(T) !+ amp_8(T) + amp_14(T)
	rslt = conjg(rslt)
	else if (helicity(i2).eq.1.and.helicity(i3).eq.1.and.helicity(i4).eq.-1) then 
	rslt = amp_3(T) !+ amp_9(T) + amp_15(T)
	rslt = conjg(rslt)
	end if
  end associate
end function

!!!! amplitude case part,++ not relevant now
	
	function amp_1(T) result(rslt)
	type(qomentum_list_type),intent(in) :: T
	complex(fltknd) :: rslt,xx,yy,zz
	integer :: i1, i2, i3, i4, i5
	!
	i1=2 ;i2=3 ;i3=4; i4=5; i5=1
	!
	rslt = 0
	xx= T%sqr(i4,i5)*T%sqr(i4,i5)* T%ang(i1,i2)
	yy = T%ang(i3,i1) - T%ang(i3,i2)*T%sqr(i2,i1)/ T%Q(i1)%kapp
	zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
	if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
	write (*,*) "In amp1"
	end function
	
	function amp_2(T) result(rslt)
	type(qomentum_list_type),intent(in) :: T
	complex(fltknd) :: rslt,tt,uu,vv,xx,ww,yy,zz,m
	integer :: i1, i2, i3, i4, i5
	!
	i1=2 ;i2=3 ;i3=4; i4=5; i5=1
	!
	m=1 ! mass of Z bozon 
	rslt = 0
	tt = T%sqr(i4,i5)*T%ang(i1,i2)
	uu = (-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i1)%kstr * T%Q(i5)%kapp
	vv = T%sqr(i5,i4)
	ww = T%ang(i4,i3) *T%sqr(i3,i4)
	xx = T%ang(i4,i1) - T%ang(i4,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
	yy = T%sqr(i5,i3) - T%ang(i5,i4)*T%sqr(i4,i3)/T%Q(i5)%kstr
	zz = T%ang(i3,i1) - T%ang(i3,i2)*T%sqr(i2,i1)/T%Q(i1)%kapp
	if (tt.ne.0.and.uu.ne.0.and.ww.ne.0) rslt = 4*tt/uu*(m*vv/ww*xx-yy*zz/m)
    write (*,*) "In amp2"
	end function
	
	function amp_3(T) result(rslt)
	type(qomentum_list_type),intent(in) :: T
	complex(fltknd) :: rslt,xx,yy,zz
	integer :: i1, i2, i3, i4, i5
	!
	i1=2 ;i2=3 ;i3=4; i4=5; i5=1
	!
	rslt = 0
	xx = T%ang(i1,i2) * T%ang(i1,i2) * T%sqr(i4,i5) 
	yy = T%sqr(i5,i3) - T%ang(i5,i4)*T%sqr(i4,i3)/T%Q(i5)%kstr
	zz = (-T%Q(i1)%kapp*T%Q(i1)%kstr-T%ang(i2,i1,i2))*(-T%Q(i5)%kapp*T%Q(i5)%kstr-T%ang(i4,i5,i4))*T%Q(i1)%kstr*T%Q(i5)%kapp
	if (xx.ne.0.and.yy.ne.0.and.zz.ne.0) rslt = 2*sqrt_2*xx*yy/zz
    write (*,*) "In amp3"
	end function
	

  
	function colorSum( Nminus2 ,a ) result(rslt)
  integer,intent(in) :: Nminus2
  complex(fltknd),intent(in) :: a(factorial(Nminus2))
  complex(fltknd) :: rslt ,z(factorial(5)) ,uu,vv,ww,xx,yy
  integer :: Nadj,i
!
  Nadj = Ncolor(2)-1
!
  select case (Nminus2)

  case (1)
    rslt = conjg(a(1))*a(1) * Nadj

  case (2)
    z(1:2) = conjg(a(1:2))
    uu = z(1)*a(1) + z(2)*a(2)
    vv =-z(1)*a(2) - z(2)*a(1)
    rslt = ( vv + Nadj*uu ) * Nadj/Ncolor(1)
   
  case (3)
    z(1:6) = conjg(a(1:6))
    uu = z(1)*a(1) + z(2)*a(2) + z(3)*a(3) + z(4)*a(4) + z(5)*a(5) + z(6)*a(6)
    vv = z(1)*( a(4) - a(2) - a(6) ) &
       + z(2)*( a(5) - a(1) - a(3) ) &
       + z(3)*( a(6) - a(2) - a(4) ) &
       + z(4)*( a(1) - a(3) - a(5) ) &
       + z(5)*( a(2) - a(4) - a(6) ) &
       + z(6)*( a(3) - a(1) - a(5) )
    ww = z(1)*( 2*a(4) + a(3) + a(5) ) &
       + z(2)*( 2*a(5) + a(4) + a(6) ) &
       + z(3)*( 2*a(6) + a(1) + a(5) ) &
       + z(4)*( 2*a(1) + a(2) + a(6) ) &
       + z(5)*( 2*a(2) + a(1) + a(3) ) &
       + z(6)*( 2*a(3) + a(2) + a(4) )
    rslt = ( ww + Nadj*( vv + Nadj*uu ) ) * Nadj/Ncolor(2)
   
  case (4)
    associate(m=>mtx_4_sqr)
    z(1:24) = conjg(a(1:24))
    uu=0 ;vv=0 ;ww=0 ;xx=0
    do i=1,24
      uu = uu + z(i)*( a(m(1,i)) )
      vv = vv + z(i)*( sum(a(m(2:7,i))) - sum(a(m(8:10,i))) )
      ww = ww + z(i)*( 2*sum(a(m(11:12,i))) + sum(a(m(13:18,i))) &
                             - sum(a(m(19:22,i))) - 3*a(m(23,i)) )
      xx = xx - z(i)*( 4*a(m(24,i)) + 2*sum(a(m(25:31,i))) + sum(a(m(32:36,i))) )
    enddo
    rslt = ( xx + Nadj*( ww + Nadj*( vv + Nadj*uu ) ) ) * Nadj/Ncolor(3)
    end associate

  case default
    rslt = 0

  end select
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




module amp4hef_qomentum
  implicit none
  private
  public :: fltknd ,imag ,init_qomentum ,zeroSlash
  public :: slashed_type ,qomentum_type ,qomentum_list_type
  public :: square ,twodot
  public :: operator(*),operator(+),operator(-)
  public :: NsizeProc,NsizeFlavor

  integer,parameter :: fltknd=kind(1d0)
  integer,parameter :: vecPerm(3)=[3,1,2]
  integer,parameter :: NsizeProc=12
  integer,parameter :: NsizeFlavor=3

  complex(fltknd),save,protected :: imag
  logical,save :: initd=.false.

  type :: sqrL_type
    private
    complex(fltknd) :: x1,x2
  end type
  type :: Rsqr_type
    private
    complex(fltknd) :: x1,x2
  end type
  type :: angL_type
    private
    complex(fltknd) :: y1,y2
  end type
  type :: Rang_type
    private
    complex(fltknd) :: y1,y2
  end type

  type :: slashed_type
    private
    complex(fltknd) :: c11,c12,c21,c22
  end type

  type(slashed_type),save,protected :: zeroSlash
  
  type :: qomentum_type
    type(sqrL_type) :: sqrL
    type(Rsqr_type) :: Rsqr
    type(angL_type) :: angL
    type(Rang_type) :: Rang
    type(slashed_type) :: p ! direction
    type(slashed_type) :: k ! momentum
    real(fltknd) :: momenta(0:3),direction(0:3)
    complex(fltknd) :: kapp,kstr
    logical :: lightlike
  contains
    procedure :: qom_fill_r
    procedure :: qom_fill_sl
    generic :: fill=>qom_fill_r,qom_fill_sl
    procedure :: minus=>qom_minus
    procedure :: prnt=>qom_prnt
  end type

  type :: qomentum_list_type
    integer :: Ntot,Noff, NZ
    integer :: onshell(NsizeProc),offshell(2), Z
    integer :: Nflavor(-NsizeFlavor:NsizeFlavor)
    integer :: flavor(NsizeProc,-NsizeFlavor:NsizeFlavor)
    integer :: symFac,sizeHel,sizePerm
    integer :: NhelOrder,helOrder(NsizeProc)
    type(qomentum_type) :: Q(NsizeProc)
    procedure(me_interface),pointer :: matrix_element
    procedure(all_amp_interface),pointer :: all_amplitudes
    procedure(amp_interface),pointer :: amplitude
  contains
    procedure :: set_direction
    procedure :: list_ang
    procedure :: list_ang_k
    procedure :: list_ang_sl
    procedure :: list_ang_l
    procedure :: list_ang_kk
    procedure :: list_ang_kl
    procedure :: list_ang_lk
    procedure :: list_ang_ll
    procedure :: list_sqr
    procedure :: list_sqr_k
    procedure :: list_sqr_sl
    procedure :: list_sqr_l
    procedure :: list_sqr_kk
    procedure :: list_sqr_kl
    procedure :: list_sqr_lk
    procedure :: list_sqr_ll
    generic :: ang=>list_ang,list_ang_k,list_ang_l ,list_ang_sl &
                   ,list_ang_kk,list_ang_kl,list_ang_lk,list_ang_ll
    generic :: sqr=>list_sqr,list_sqr_k,list_sqr_l ,list_sqr_sl &
                   ,list_sqr_kk,list_sqr_kl,list_sqr_lk,list_sqr_ll
    procedure :: sinv=>list_sinv_2
    procedure :: transvec=>list_transvec
    procedure :: sumk=>list_sumk
    procedure :: shift_ang=>list_shift_ang
    procedure :: shift_sqr=>list_shift_sqr
    procedure :: shift_kapp=>list_shift_kapp
    procedure :: shift_kstr=>list_shift_kstr
    procedure :: cterm=>list_cterm
    procedure :: dterm=>list_dterm
  end type


  interface operator (*)
    module procedure sqrL_Rsqr ,angL_Rang ,Rsqr_angL ,Rang_sqrL
    module procedure sqrL_slash ,angL_slash ,slash_Rsqr ,slash_Rang
    module procedure c_sqrL ,c_Rsqr ,c_angL ,c_Rang ,c_slash
    module procedure sqrL_c ,Rsqr_c ,angL_c ,Rang_c ,slash_c
  end interface
  interface operator (+)
    module procedure plus_sl_sl
  end interface
  interface operator (-)
    module procedure minus_sl_sl
  end interface

  interface square
    module procedure square_sl
  end interface

  interface twodot
    module procedure twodot_sl_sl
  end interface


  abstract interface
    function me_interface( obj ) result(rslt)
    import
    class(qomentum_list_type),intent(in) :: obj
    real(fltknd) :: rslt
    end function
  end interface 

  abstract interface
    subroutine all_amp_interface( obj ,NhelConf ,Nperm ,amplitude ,factor )
    import
    class(qomentum_list_type),intent(in) :: obj
    integer,intent(out) :: NhelConf,Nperm
    complex(fltknd),intent(out) :: amplitude(:,:),factor
    end subroutine
  end interface

  abstract interface
    function amp_interface( obj ,helicity ,perm ) result(rslt)
    import
    class(qomentum_list_type),intent(in) :: obj
    integer,intent(in) :: helicity(:),perm(:)
    complex(fltknd) :: rslt
    end function
  end interface 


contains


  subroutine init_qomentum
! initialize module
  if (initd) return ;initd=.true.
  imag = cmplx(0,1,kind=kind(imag))
  zeroSlash%c11 = 0
  zeroSlash%c12 = 0
  zeroSlash%c21 = 0
  zeroSlash%c22 = 0
  end subroutine


  subroutine qom_fill_r( obj ,momentum ,direction )
  class(qomentum_type) :: obj
  real(fltknd),intent(in) :: momentum(0:3),direction(0:3)
  optional :: direction
  complex(fltknd) :: cmpnnt2ima
  if (present(direction)) then
    obj%lightlike = .false.
    cmpnnt2ima = direction(vecPerm(2))*imag
    obj%p%c11 = direction(        0 ) + direction(vecPerm(3))
    obj%p%c12 = direction(vecPerm(1)) - cmpnnt2ima
    obj%p%c21 = direction(vecPerm(1)) + cmpnnt2ima
    call finish_p
    cmpnnt2ima = momentum(vecPerm(2))*imag
    obj%k%c11 = momentum(        0 ) + momentum(vecPerm(3))
    obj%k%c12 = momentum(vecPerm(1)) - cmpnnt2ima
    obj%k%c21 = momentum(vecPerm(1)) + cmpnnt2ima
    obj%k%c22 = momentum(        0 ) - momentum(vecPerm(3))
  else
    obj%lightlike = .true.
    cmpnnt2ima = momentum(vecPerm(2))*imag
    obj%p%c11 = momentum(        0 ) + momentum(vecPerm(3))
    obj%p%c12 = momentum(vecPerm(1)) - cmpnnt2ima
    obj%p%c21 = momentum(vecPerm(1)) + cmpnnt2ima
    call finish_p
    obj%k = obj%p
    obj%kstr = 0
    obj%kapp = 0
  endif
  contains
    subroutine finish_p
    obj%p%c22 = obj%p%c12*obj%p%c21/obj%p%c11 !p(l0)-p(l3) !
    call put_spinors( obj )
    end subroutine
  end subroutine

  subroutine qom_fill_sl( obj ,pp )
  class(qomentum_type) :: obj
  type(slashed_type),intent(in) :: pp
  obj%lightlike = .true.
  obj%k = pp
  obj%p = pp
  call put_spinors( obj )
  end subroutine

  subroutine put_spinors( obj )
  class(qomentum_type) :: obj
  complex(fltknd) :: hh
  hh = sqrt(abs(obj%p%c11))
  obj%sqrL%x1 = obj%p%c11/hh
  obj%sqrL%x2 = obj%p%c12/hh
  obj%Rsqr%x1 =-obj%sqrL%x2
  obj%Rsqr%x2 = obj%sqrL%x1
  obj%Rang%y1 = hh
  obj%Rang%y2 = obj%p%c21*hh/obj%p%c11
  obj%angL%y1 =-obj%Rang%y2
  obj%angL%y2 = hh
  end subroutine


  subroutine qom_prnt( obj ,iunit )
  class(qomentum_type) :: obj
  integer,intent(in),optional :: iunit
  integer :: wunit
  wunit = 0
  if (present(iunit)) wunit = iunit
  if (wunit.le.0) return
  write(wunit,*) 'p%c11',obj%p%c11
  write(wunit,*) 'p%c12',obj%p%c12
  write(wunit,*) 'p%c21',obj%p%c21
  write(wunit,*) 'p%c22',obj%p%c22
  write(wunit,*) 'sqrL%x1',obj%sqrL%x1
  write(wunit,*) 'sqrL%x2',obj%sqrL%x2
  write(wunit,*) 'Rsqr%x1',obj%Rsqr%x1
  write(wunit,*) 'Rsqr%x2',obj%Rsqr%x2
  write(wunit,*) 'Rang%y1',obj%Rang%y1
  write(wunit,*) 'Rang%y2',obj%Rang%y2
  write(wunit,*) 'angL%y1',obj%angL%y1
  write(wunit,*) 'angL%y2',obj%angL%y2
  write(wunit,*) 'k%c11',obj%k%c11
  write(wunit,*) 'k%c12',obj%k%c12
  write(wunit,*) 'k%c21',obj%k%c21
  write(wunit,*) 'k%c22',obj%k%c22
  write(wunit,*) 'kapp',obj%kapp
  write(wunit,*) 'kstr',obj%kstr
  write(wunit,*) 'lightlike',obj%lightlike
  end subroutine


  function sqrL_Rsqr( a,b ) result(c)
! [ab]
  type(sqrL_type),intent(in) :: a
  type(Rsqr_type),intent(in) :: b
  complex(fltknd) :: c
  c = a%x1*b%x1 + a%x2*b%x2
  end function

  function angL_Rang( a,b ) result(c)
! <ab>
  type(angL_type),intent(in) :: a
  type(Rang_type),intent(in) :: b
  complex(fltknd) :: c
  c = a%y1*b%y1 + a%y2*b%y2
  end function

  function Rsqr_angL( a,b ) result(c)
! |a]<b|
  type(Rsqr_type),intent(in) :: a
  type(angL_type),intent(in) :: b
  type(slashed_type) :: c
  c%c11 = a%x2*b%y2
  c%c12 =-a%x1*b%y2
  c%c21 =-a%x2*b%y1
  c%c22 = a%x1*b%y1
  end function

  function Rang_sqrL( a,b ) result(c)
! |a>[b|
  type(Rang_type),intent(in) :: a
  type(sqrL_type),intent(in) :: b
  type(slashed_type) :: c
  c%c11 = a%y1*b%x1
  c%c12 = a%y1*b%x2
  c%c21 = a%y2*b%x1
  c%c22 = a%y2*b%x2
  end function

  function sqrL_slash( a,b ) result(c)
! [a|b_slash
  type(sqrL_type),intent(in) :: a
  type(slashed_type),intent(in) :: b
  type(angL_type) :: c
  c%y1 = a%x1*b%c22 - a%x2*b%c21
  c%y2 =-a%x1*b%c12 + a%x2*b%c11
  end function

  function angL_slash( a,b ) result(c)
! <a|b_slash
  type(angL_type),intent(in) :: a
  type(slashed_type),intent(in) :: b
  type(sqrL_type) :: c
  c%x1 = a%y1*b%c11 + a%y2*b%c21
  c%x2 = a%y1*b%c12 + a%y2*b%c22
  end function

  function slash_Rsqr( a,b ) result(c)
! a_slash|b]
  type(slashed_type),intent(in) :: a
  type(Rsqr_type),intent(in) :: b
  type(Rang_type) :: c
  c%y1 = a%c11*b%x1 + a%c12*b%x2
  c%y2 = a%c21*b%x1 + a%c22*b%x2
  end function

  function slash_Rang( a,b ) result(c)
! a_slash|b>
  type(slashed_type),intent(in) :: a
  type(Rang_type),intent(in) :: b
  type(Rsqr_type) :: c
  c%x1 = a%c22*b%y1 - a%c12*b%y2
  c%x2 =-a%c21*b%y1 + a%c11*b%y2
  end function

  function plus_sl_sl( a,b ) result(c)
! a_slash + b_slash
  type(slashed_type),intent(in) :: a,b
  type(slashed_type) :: c
  c%c11 = a%c11+b%c11
  c%c12 = a%c12+b%c12
  c%c21 = a%c21+b%c21
  c%c22 = a%c22+b%c22
  end function
  
  function minus_sl_sl( a,b ) result(c)
! a_slash - b_slash
  type(slashed_type),intent(in) :: a,b
  type(slashed_type) :: c
  c%c11 = a%c11-b%c11
  c%c12 = a%c12-b%c12
  c%c21 = a%c21-b%c21
  c%c22 = a%c22-b%c22
  end function
  
  function c_slash( a,b ) result(c)
! scalar * b_slash
  complex(fltknd),intent(in) :: a
  type(slashed_type),intent(in) :: b
  type(slashed_type) :: c
  c%c11 = a*b%c11
  c%c12 = a*b%c12
  c%c21 = a*b%c21
  c%c22 = a*b%c22
  end function
  
  function slash_c( b,a ) result(c)
! b_slash * scalar
  complex(fltknd),intent(in) :: a
  type(slashed_type),intent(in) :: b
  type(slashed_type) :: c
  c%c11 = a*b%c11
  c%c12 = a*b%c12
  c%c21 = a*b%c21
  c%c22 = a*b%c22
  end function
  
  function c_sqrL( a,b ) result(c)
! scalar * [b|
  complex(fltknd),intent(in) :: a
  type(sqrL_type),intent(in) :: b
  type(sqrL_type) :: c
  c%x1 = a*b%x1
  c%x2 = a*b%x2
  end function
  function sqrL_c( b,a ) result(c)
! [b| * scalar
  complex(fltknd),intent(in) :: a
  type(sqrL_type),intent(in) :: b
  type(sqrL_type) :: c
  c%x1 = a*b%x1
  c%x2 = a*b%x2
  end function
  
  function c_Rsqr( a,b ) result(c)
! scalar * |b]
  complex(fltknd),intent(in) :: a
  type(Rsqr_type),intent(in) :: b
  type(Rsqr_type) :: c
  c%x1 = a*b%x1
  c%x2 = a*b%x2
  end function
  function Rsqr_c( b,a ) result(c)
! |b] * scalar
  complex(fltknd),intent(in) :: a
  type(Rsqr_type),intent(in) :: b
  type(Rsqr_type) :: c
  c%x1 = a*b%x1
  c%x2 = a*b%x2
  end function
  
  function c_angL( a,b ) result(c)
! scalar * <b|
  complex(fltknd),intent(in) :: a
  type(angL_type),intent(in) :: b
  type(angL_type) :: c
  c%y1 = a*b%y1
  c%y2 = a*b%y2
  end function
  function angL_c( b,a ) result(c)
! <b| * scalar
  complex(fltknd),intent(in) :: a
  type(angL_type),intent(in) :: b
  type(angL_type) :: c
  c%y1 = a*b%y1
  c%y2 = a*b%y2
  end function
  
  function c_Rang( a,b ) result(c)
! scalar * |b>
  complex(fltknd),intent(in) :: a
  type(Rang_type),intent(in) :: b
  type(Rang_type) :: c
  c%y1 = a*b%y1
  c%y2 = a*b%y2
  end function
  function Rang_c( b,a ) result(c)
! |b> * scalar
  complex(fltknd),intent(in) :: a
  type(Rang_type),intent(in) :: b
  type(Rang_type) :: c
  c%y1 = a*b%y1
  c%y2 = a*b%y2
  end function



  function square_sl( a ) result(rslt)
  type(slashed_type),intent(in) :: a
  complex(fltknd) :: rslt
  rslt = a%c11*a%c22 - a%c12*a%c21
  end function


  function twodot_sl_sl( a,b ) result(c)
  type(slashed_type),intent(in) :: a
  type(slashed_type),intent(in) :: b
  complex(fltknd) :: c
  c = a%c11*b%c22 + b%c11*a%c22 &
    - a%c12*b%c21 - b%c12*a%c21
  end function


! Below the type-bound routines for the qomentum_list_type, that
! ask for integer arguments, 
! eg. l%ang(1,2)= l%p(1)%angL*l%p(2)%Rang etc.

  subroutine set_direction( obj ,i3, i1 )
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i3
  if (.not.obj%Q(i3)%lightlike) then
    obj%Q(i3)%p = obj%Q(i1)%p
  end if
  end subroutine



  function list_ang( obj ,i1,i2 ) result(rslt)
! <12>
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2
  complex(fltknd) :: rslt
  rslt = obj%Q(i1)%angL*obj%Q(i2)%Rang
  end function
  
  function list_ang_k( obj ,i1 ,i3 ,i2 ) result(rslt)
! <1|k3|2]
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3
  complex(fltknd) :: rslt
  rslt = obj%Q(i1)%angL * obj%Q(i3)%k * obj%Q(i2)%Rsqr
  end function
  
  function list_ang_sl( obj ,i1 ,k3 ,i2 ) result(rslt)
! <1|k3|2]
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2
  type(slashed_type),intent(in) :: k3
  complex(fltknd) :: rslt
  rslt = obj%Q(i1)%angL * k3 * obj%Q(i2)%Rsqr
  end function
  
  function list_ang_l( obj ,i1 ,i3 ,i2 ) result(rslt)
! <1|k3(1)+k3(2)+...|2]
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3(:)
  complex(fltknd) :: rslt
  type(slashed_type) :: kk
  integer :: ii
  kk = zeroSlash
  do ii=1,size(i3)
    kk = kk + obj%Q(i3(ii))%k
  enddo
  rslt = obj%Q(i1)%angL * kk * obj%Q(i2)%Rsqr
  end function
  
  function list_ang_kk( obj ,i1 ,i3,i4 ,i2 ) result(rslt)
! <1|k3*k4|2>
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3,i4
  complex(fltknd) :: rslt
  rslt = obj%Q(i1)%angL * obj%Q(i3)%k * obj%Q(i4)%k * obj%Q(i2)%Rang
  end function
  
  function list_ang_kl( obj ,i1 ,i3,i4 ,i2 ) result(rslt)
! <1|k3*(k4(1)+k4(2)+...)|2>
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3,i4(:)
  complex(fltknd) :: rslt
  type(slashed_type) :: kk
  integer :: ii
  kk = zeroSlash
  do ii=1,size(i4)
    kk = kk + obj%Q(i4(ii))%k
  enddo
  rslt = obj%Q(i1)%angL * obj%Q(i3)%k * kk * obj%Q(i2)%Rang
  end function
  
  function list_ang_lk( obj ,i1 ,i3,i4 ,i2 ) result(rslt)
! <1|(k3(1)+k3(2)+...)*k4|2>
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3(:),i4
  complex(fltknd) :: rslt
  type(slashed_type) :: kk
  integer :: ii
  kk = zeroSlash
  do ii=1,size(i3)
    kk = kk + obj%Q(i3(ii))%k
  enddo
  rslt = obj%Q(i1)%angL * kk * obj%Q(i4)%k * obj%Q(i2)%Rang
  end function
  
  function list_ang_ll( obj ,i1 ,i3,i4 ,i2 ) result(rslt)
! <1|(k3(1)+k3(2)+...)*(k4(1)+k4(2)+...)|2>
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3(:),i4(:)
  complex(fltknd) :: rslt
  type(slashed_type) :: k3,k4
  integer :: ii
  k3 = zeroSlash
  do ii=1,size(i3)
    k3 = k3 + obj%Q(i3(ii))%k
  enddo
  k4 = zeroSlash
  do ii=1,size(i4)
    k4 = k4 + obj%Q(i4(ii))%k
  enddo
  rslt = obj%Q(i1)%angL * k3 * k4 * obj%Q(i2)%Rang
  end function
  


  function list_sqr( obj ,i1,i2 ) result(rslt)
! [12]
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2
  complex(fltknd) :: rslt
  rslt = obj%Q(i1)%sqrL*obj%Q(i2)%Rsqr
  end function
  
  function list_sqr_k( obj ,i1 ,i3 ,i2 ) result(rslt)
! [1|k3|2>
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3
  complex(fltknd) :: rslt
  rslt = obj%Q(i1)%sqrL * obj%Q(i3)%k * obj%Q(i2)%Rang
  end function
  
  function list_sqr_sl( obj ,i1 ,k3 ,i2 ) result(rslt)
! [1|k3|2>
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2
  type(slashed_type),intent(in) :: k3
  complex(fltknd) :: rslt
  rslt = obj%Q(i1)%sqrL * k3 * obj%Q(i2)%Rang
  end function
  
  function list_sqr_l( obj ,i1 ,i3 ,i2 ) result(rslt)
! [1|k3(1)+k3(2)+...|2>
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3(:)
  complex(fltknd) :: rslt
  type(slashed_type) :: kk
  integer :: ii
  kk = zeroSlash
  do ii=1,size(i3)
    kk = kk + obj%Q(i3(ii))%k
  enddo
  rslt = obj%Q(i1)%sqrL * kk * obj%Q(i2)%Rang
  end function
  
  function list_sqr_kk( obj ,i1 ,i3,i4 ,i2 ) result(rslt)
! [1|k3*k4|2]
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3,i4
  complex(fltknd) :: rslt
  rslt = obj%Q(i1)%sqrL * obj%Q(i3)%k * obj%Q(i4)%k * obj%Q(i2)%Rsqr
  end function
  
  function list_sqr_kl( obj ,i1 ,i3,i4 ,i2 ) result(rslt)
! [1|k3*(k4(1)+k4(2)+...)|2]
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3,i4(:)
  complex(fltknd) :: rslt
  type(slashed_type) :: kk
  integer :: ii
  kk = zeroSlash
  do ii=1,size(i4)
    kk = kk + obj%Q(i4(ii))%k
  enddo
  rslt = obj%Q(i1)%sqrL * obj%Q(i3)%k * kk * obj%Q(i2)%Rsqr
  end function
  
  function list_sqr_lk( obj ,i1 ,i3,i4 ,i2 ) result(rslt)
! [1|(k3(1)+k3(2)+...)*k4|2]
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3(:),i4
  complex(fltknd) :: rslt
  type(slashed_type) :: kk
  integer :: ii
  kk = zeroSlash
  do ii=1,size(i3)
    kk = kk + obj%Q(i3(ii))%k
  enddo
  rslt = obj%Q(i1)%sqrL * kk * obj%Q(i4)%k * obj%Q(i2)%Rsqr
  end function
  
  function list_sqr_ll( obj ,i1 ,i3,i4 ,i2 ) result(rslt)
! [1|(k3(1)+k3(2)+...)*(k4(1)+k4(2)+...)|2]
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2,i3(:),i4(:)
  complex(fltknd) :: rslt
  type(slashed_type) :: k3,k4
  integer :: ii
  k3 = zeroSlash
  do ii=1,size(i3)
    k3 = k3 + obj%Q(i3(ii))%k
  enddo
  k4 = zeroSlash
  do ii=1,size(i4)
    k4 = k4 + obj%Q(i4(ii))%k
  enddo
  rslt = obj%Q(i1)%sqrL * k3 * k4 * obj%Q(i2)%Rsqr
  end function
  

  function list_sinv_2( obj ,i1,i2 ) result(rslt)
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2
  complex(fltknd) :: rslt ,aa,bb,cc,dd
  aa = obj%Q(i1)%k%c11 + obj%Q(i2)%k%c11
  bb = obj%Q(i1)%k%c22 + obj%Q(i2)%k%c22
  cc = obj%Q(i1)%k%c12 + obj%Q(i2)%k%c12
  dd = obj%Q(i1)%k%c21 + obj%Q(i2)%k%c21
  rslt = aa*bb - cc*dd
  end function


  function list_sumk( obj ,lst ) result(rslt)
  class(qomentum_list_type) :: obj
  integer,intent(in) :: lst(:)
  type(slashed_type) :: rslt
  integer :: ii
  rslt = zeroSlash
  do ii=1,size(lst)
    rslt = rslt + obj%Q(lst(ii))%k
  enddo
  end function


  function qom_minus( obj ) result(rslt)
! multiply qomentum_type with -1, also affects the spinors
  class(qomentum_type) :: obj
  type(qomentum_type) :: rslt
  rslt%lightlike = obj%lightlike
  rslt%k%c11 =-obj%k%c11
  rslt%k%c12 =-obj%k%c12
  rslt%k%c21 =-obj%k%c21
  rslt%k%c22 =-obj%k%c22
  rslt%p%c11 =-obj%p%c11
  rslt%p%c12 =-obj%p%c12
  rslt%p%c21 =-obj%p%c21
  rslt%p%c22 =-obj%p%c22
  rslt%sqrL%x1 =-obj%sqrL%x1
  rslt%sqrL%x2 =-obj%sqrL%x2
  rslt%Rsqr%x1 =-obj%Rsqr%x1
  rslt%Rsqr%x2 =-obj%Rsqr%x2
  rslt%Rang%y1 = obj%Rang%y1
  rslt%Rang%y2 = obj%Rang%y2
  rslt%angL%y1 = obj%angL%y1
  rslt%angL%y2 = obj%angL%y2
  end function

  
  function list_transvec( obj ,i1,i2 ) result(rslt)
! <1|gamma_mu|2] / 2
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2
  type(qomentum_type) :: rslt
  rslt%lightlike = .true.
  rslt%angL = obj%Q(i1)%angL
  rslt%Rang = obj%Q(i1)%Rang
  rslt%sqrL = obj%Q(i2)%sqrL
  rslt%Rsqr = obj%Q(i2)%Rsqr
  rslt%p = rslt%Rsqr * rslt%angL
  rslt%k = rslt%p
  end function


  function list_shift_sqr( obj ,i1,i2 ,zz ) result(rslt)
! construct qomentum_type obj by shifting its square spinors
! with ee, so  |1> -> |1>, |1] -> |1] + z*|2]
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2
  complex(fltknd),intent(in) :: zz
  type(qomentum_type) :: rslt
  rslt%lightlike = .true.
  rslt%sqrL%x1 = obj%Q(i1)%sqrL%x1 + zz*obj%Q(i2)%sqrL%x1
  rslt%sqrL%x2 = obj%Q(i1)%sqrL%x2 + zz*obj%Q(i2)%sqrL%x2
  rslt%Rsqr%x1 =-rslt%sqrL%x2
  rslt%Rsqr%x2 = rslt%sqrL%x1
  rslt%Rang%y1 = obj%Q(i1)%Rang%y1
  rslt%Rang%y2 = obj%Q(i1)%Rang%y2
  rslt%angL%y1 = obj%Q(i1)%angL%y1
  rslt%angL%y2 = obj%Q(i1)%angL%y2
  rslt%p = rslt%Rsqr*rslt%angL
  rslt%k = rslt%p
  end function


  function list_shift_ang( obj ,i1,i2 ,zz ) result(rslt)
! construct qomentum_type obj by shifting its angular spinors
! with ee, so  |1] -> |1], |1> -> |1> + z*|2>
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2
  complex(fltknd),intent(in) :: zz
  type(qomentum_type) :: rslt
  rslt%lightlike = .true.
  rslt%Rang%y1 = obj%Q(i1)%Rang%y1 + zz*obj%Q(i2)%Rang%y1
  rslt%Rang%y2 = obj%Q(i1)%Rang%y2 + zz*obj%Q(i2)%Rang%y2
  rslt%angL%y1 =-rslt%Rang%y2
  rslt%angL%y2 = rslt%Rang%y1
  rslt%sqrL%x1 = obj%Q(i1)%sqrL%x1
  rslt%sqrL%x2 = obj%Q(i1)%sqrL%x2
  rslt%Rsqr%x1 = obj%Q(i1)%Rsqr%x1
  rslt%Rsqr%x2 = obj%Q(i1)%Rsqr%x2
  rslt%p = rslt%Rsqr*rslt%angL
  rslt%k = rslt%p
  end function


  function list_shift_kapp( obj ,i1,i2 ,zz ) result(rslt)
! only shift momentum
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2
  complex(fltknd),intent(in) :: zz
  type(qomentum_type) :: rslt
  rslt%lightlike = .false.
  rslt%Rang = obj%Q(i1)%Rang
  rslt%angL = obj%Q(i1)%angL
  rslt%sqrL = obj%Q(i1)%sqrL
  rslt%Rsqr = obj%Q(i1)%Rsqr
  rslt%p    = obj%Q(i1)%p
  rslt%k    = obj%Q(i1)%k + zz * obj%Q(i1)%Rang*obj%Q(i2)%sqrL
  rslt%kapp = obj%Q(i1)%kapp - zz * obj%Q(i1)%sqrL*obj%Q(i2)%Rsqr
  rslt%kstr = obj%Q(i1)%kstr
  end function


  function list_shift_kstr( obj ,i1,i2 ,zz ) result(rslt)
! only shift momentum
  class(qomentum_list_type) :: obj
  integer,intent(in) :: i1,i2
  complex(fltknd),intent(in) :: zz
  type(qomentum_type) :: rslt
  rslt%lightlike = .false.
  rslt%Rang = obj%Q(i1)%Rang
  rslt%angL = obj%Q(i1)%angL
  rslt%sqrL = obj%Q(i1)%sqrL
  rslt%Rsqr = obj%Q(i1)%Rsqr
  rslt%p    = obj%Q(i1)%p
  rslt%k    = obj%Q(i1)%k + zz * obj%Q(i2)%Rang*obj%Q(i1)%sqrL
  rslt%kstr = obj%Q(i1)%kstr - zz * obj%Q(i2)%angL*obj%Q(i1)%Rang
  rslt%kapp = obj%Q(i1)%kapp
  end function


  subroutine list_cterm( obj ,i1,i2 ,xx )
! Off-shell obj%Q(i1) goes on-shell 
  class(qomentum_list_type) :: obj
  complex(fltknd)     ,intent(out) :: xx
  integer,intent(in) :: i1,i2
  complex(fltknd) :: sqr12,sqrtx
  associate( Q1=>obj%Q(i1) ,Q2=>obj%Q(i2) )
  sqr12 = Q1%sqrL * Q2%Rsqr
  xx = twodot(Q2%p,Q1%k)/twodot(Q2%p,Q1%p)
  sqrtx = sqrt(xx)
!
  if (Q2%lightlike) then
    Q2%Rang = ((Q2%p+Q1%k) * Q1%Rsqr) * (-1/sqr12)
    Q2%angL%y1 =-Q2%Rang%y2
    Q2%angL%y2 = Q2%Rang%y1
    Q2%p = Q2%Rsqr * Q2%angL
    Q2%k = Q2%p
    Q2%kstr = 0
    Q2%kapp = 0
  else
    Q2%kstr = (Q2%angL*(Q2%k+Q1%k)*Q1%Rsqr)*(-1/sqr12)
    Q2%kapp = Q2%kapp
    Q2%k    = Q2%k - (Q1%kapp/sqr12)*(Q2%Rsqr*Q1%angL)
  endif
!
  Q1%lightlike = .true.
  Q1%Rang = (Q1%k * Q2%Rsqr) * (1/(sqrtx*sqr12))
  Q1%sqrL = Q1%sqrL * sqrtx
  Q1%angL%y1 =-Q1%Rang%y2
  Q1%angL%y2 = Q1%Rang%y1
  Q1%Rsqr%x1 =-Q1%sqrL%x2
  Q1%Rsqr%x2 = Q1%sqrL%x1
  Q1%p = Q1%Rsqr * Q1%angL
  Q1%k = Q1%p
  Q1%kstr = 0
  Q1%kapp = 0
  end associate
  end subroutine


  subroutine list_dterm( obj ,i1,i2 ,xx )
! Off-shell obj%Q(i2) goes on-shell
  class(qomentum_list_type) :: obj
  complex(fltknd)     ,intent(out) :: xx
  integer,intent(in) :: i1,i2
  complex(fltknd) :: ang21,sqrtx
  associate( Q1=>obj%Q(i1) ,Q2=>obj%Q(i2) )
  ang21 = Q2%angL * Q1%Rang
  xx = twodot(Q1%p,Q2%k)/twodot(Q1%p,Q2%p)
  sqrtx = sqrt(xx)
!
  if (Q1%lightlike) then
    Q1%sqrL = (1/ang21)*Q2%angL*(Q1%p+Q2%k)
    Q1%Rsqr%x1 =-Q1%sqrL%x2
    Q1%Rsqr%x2 = Q1%sqrL%x1
    Q1%p = Q1%Rsqr*Q1%angL
    Q1%k = Q1%p
    Q1%kstr = 0
    Q1%kapp = 0
  else
    Q1%kstr = Q1%kstr
    Q1%kapp = (1/ang21)*(Q2%angL*(Q1%k+Q2%k)*Q1%Rsqr)
    Q1%k    = Q1%k + (Q2%kstr/ang21)*(Q2%Rsqr*Q1%angL)
  endif
 !
  Q2%lightlike = .true.
  Q2%sqrL = (-1/(sqrtx*ang21))*Q1%angL*Q2%k
  Q2%Rang = sqrtx*Q2%Rang
  Q2%angL%y1 =-Q2%Rang%y2
  Q2%angL%y2 = Q2%Rang%y1
  Q2%Rsqr%x1 =-Q2%sqrL%x2
  Q2%Rsqr%x2 = Q2%sqrL%x1
  Q2%p = Q2%Rsqr*Q2%angL
  Q2%k = Q2%p
  Q2%kstr = 0
  Q2%kapp = 0
  end associate
  end subroutine


end module



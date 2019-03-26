module amp4hef_ng
  use amp4hef_io
  use amp4hef_aux
  use amp4hef_qomentum

  implicit none
  private
  public :: fill_matrices_ng,matrix_element_ng ,amplitude_ng ,all_amplitudes_ng
  public :: print_color_matrix


  integer,allocatable,save :: mtx_4_sqr(:,:),mtx_5_sqr(:,:) &
                             ,mtx_6_sqr(:,:),mtx_7_sqr(:,:)


contains


  function matrix_element_ng( Tin ) result(rslt)
  class(qomentum_list_type),intent(in) :: Tin
  real(fltknd) :: rslt
  complex(fltknd) :: amp(120)
  integer :: ii,jj,Nons,Nminus2
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )
  Nons = Ntot-Noff
  Nminus2 = Ntot-2
  rslt = 0
  do ii=1,twoPower(Nons-1)
    do jj=1,factorial(Nminus2)
      amp(jj) = amplitude_ng( Tin ,helTable(1:Nons,ii) ,perTable(1:Nminus2,jj) )
    enddo
    rslt = rslt + colorSum( Nminus2 ,amp(1:factorial(Nminus2)) )
  enddo
  rslt = rslt*2
  do ii=1,Noff
    rslt = rslt * (Tin%Q(Tin%offshell(ii))%kapp*Tin%Q(Tin%offshell(ii))%kstr)
  enddo
  end associate
  end function


  subroutine all_amplitudes_ng( Tin ,NhelConf ,Nperm ,amplitude ,factor )
  class(qomentum_list_type),intent(in) :: Tin
  integer,intent(out) :: NhelConf,Nperm
  complex(fltknd),intent(out) :: amplitude(:,:),factor
  integer :: ii,jj,Nons,Nminus2
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )
  Nons = Ntot-Noff
  Nminus2 = Ntot-2
  NhelConf = twoPower(Nons-1)
  Nperm = factorial(Nminus2)
  do ii=1,NhelConf
  do jj=1,Nperm
    amplitude(ii,jj) = amplitude_ng( Tin ,helTable(1:Nons,ii) ,perTable(1:Nminus2,jj) )
  enddo
  enddo
  factor = 2 ! only half of the helicity configurations is returned
  do ii=1,Noff
    factor = factor * (Tin%Q(Tin%offshell(ii))%kapp*Tin%Q(Tin%offshell(ii))%kstr)
  enddo
  end associate
  end subroutine


  function amplitude_ng( Tin ,helicity ,perm ) result(rslt)
! helicity(1:Ntot-Noff) refer to the on-shell gluons.
! perm(1:Ntot-2) permutes gluon 3 to Ntot.
  class(qomentum_list_type),intent(in) :: Tin
  integer,intent(in) :: helicity(:),perm(:)
  type(qomentum_list_type) :: T
  complex(fltknd) :: rslt
  integer :: hel(NsizeProc) ,mhv_type,tags(2),ii
  T%Ntot = Tin%Ntot
  T%Noff = Tin%Noff
  do ii=1,T%Noff
    hel(ii) = 0
  enddo
  do ii=T%Noff+1,2
    hel(ii) = helicity(ii-T%Noff)
  enddo
  do ii=3,T%Ntot
    hel(ii) = helicity(perm(ii-2)+2-T%Noff)
  enddo
  call determine_mhv( mhv_type ,tags ,T%Ntot ,hel )
  if (mhv_type.eq.0) then
    rslt = 0
  else
    do ii=1,T%Noff
      T%Q(ii) = Tin%Q(Tin%offshell(ii))
    enddo
    do ii=T%Noff+1,2
      T%Q(ii) = Tin%Q(Tin%onshell(ii-T%Noff))
    enddo
    do ii=3,T%Ntot
      T%Q(ii) = Tin%Q(Tin%onshell(perm(ii-2)+2-T%Noff))
    enddo
    if     (mhv_type.eq. 1) then ; rslt = amp_mostly_plus( tags ,T )
    elseif (mhv_type.eq.-1) then ; rslt = amp_mostly_minus( tags ,T )
    elseif (T%Noff.eq.0) then ; rslt = amp_0( mhv_type ,tags ,hel ,T )
    elseif (T%Noff.eq.1) then ; rslt = amp_1( mhv_type ,tags ,hel ,T )
    elseif (T%Noff.eq.2) then ; rslt = amp_2( mhv_type ,tags ,hel ,T )
    else
      write(*,*) 'ERROR in amp4hef_ng: Noff=',T%Noff,' not implemented, returning 0'
      rslt = 0
    endif
  endif
  end function 


  subroutine determine_mhv( mhv_type ,tags ,Ntot ,helicity )
! Determine if input leads to a vanishing amplitude, or an MHV amplitude.
  integer,intent(out) :: mhv_type ,tags(2)
  integer,intent(in) :: Ntot
  integer,intent(in) :: helicity(Ntot)
  integer :: ii,Npls,Nmin,Noff,plus(NsizeProc),mins(NsizeProc),offs(NsizeProc)
  Npls = 0
  Nmin = 0
  Noff = 0
  do ii=1,Ntot
    if (helicity(ii).gt.0) then
      Npls = Npls+1
      plus(Npls) = ii
    elseif (helicity(ii).lt.0) then
      Nmin = Nmin+1
      mins(Nmin) = ii
    else
      Noff = Noff+1
      offs(Noff) = ii
    endif
  enddo
  if     ((Noff.le.1).and.((Nmin.eq.0).or.(Npls.eq.0))) then !vanishing
    mhv_type = 0
  elseif ((Noff+Nmin.le.1).and.(Npls.ge.3)) then !vanishing
    mhv_type = 0
  elseif ((Noff+Npls.le.1).and.(Nmin.ge.3)) then !vanishing
    mhv_type = 0
  elseif ((Noff+Nmin).eq.2) then !mostly_plus
    mhv_type = 1
    tags(1:Noff) = offs(1:Noff)
    tags(Noff+1:Noff+Nmin) = mins(1:Nmin)
  elseif ((Noff+Npls).eq.2) then !mostly_minus
    mhv_type =-1
    tags(1:Noff) = offs(1:Noff)
    tags(Noff+1:Noff+Npls) = plus(1:Npls)
  else
    mhv_type = 999
  endif
  end subroutine


  recursive function amp_0( mhv_type ,tags ,helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: mhv_type ,tags(2) ,helicity(T%Ntot)
  complex(fltknd) :: rslt ,Ksqr,zz,Urslt,Vrslt
  type(slashed_type) :: Kslash
  type(qomentum_list_type) :: U,V
  type(qomentum_type) :: ee,Khat 
  integer :: cyc(NsizeProc),jj,ii,ihel
  integer :: U_mhv,U_tags(2),U_h(NsizeProc)
  integer :: V_mhv,V_tags(2),V_h(NsizeProc)
  associate( N=>T%Ntot )
  if     (mhv_type.eq. 0) then ; rslt = 0
  elseif (mhv_type.eq. 1) then ; rslt = amp_mostly_plus( tags ,T )
  elseif (mhv_type.eq.-1) then ; rslt = amp_mostly_minus( tags ,T )
  else
    rslt = 0
! Rotate until appropriate helicities for leg 1 and N are found.
    cyc(1:N) = unity(1:N)
    do while (helicity(cyc(1)).eq.helicity(cyc(N)).or.helicity(cyc(1)).gt.0)
      jj = cyc(1)
      cyc(1:N-1) = cyc(2:N)
      cyc(N) = jj
    enddo
    ee = T%transvec( cyc(1) ,cyc(N) )
! Sum over partitions of external gluons
    do ii=2,N-2
      jj = N-ii
      if (ii.lt.jj) then
        Kslash = T%sumk(cyc(1:ii))
      else
        Kslash = zeroSlash-T%sumk(cyc(ii+1:N))
      endif
      Ksqr = square(Kslash)
      zz =-Ksqr/T%ang(cyc(1),Kslash,cyc(N))
      call Khat%fill( zz*ee%p + Kslash )
! Sum over internal helicities
      do ihel=-1,1,2
        U%Ntot = ii+1
        U_h(1:ii) = helicity(cyc(1:ii))
        U_h(ii+1) = ihel
        call determine_mhv( U_mhv ,U_tags ,U%Ntot ,U_h )
        if (U_mhv.eq.0) cycle
        V%Ntot = jj+1
        V_h(1) =-ihel
        V_h(2:jj+1) = helicity(cyc(ii+1:N))
        call determine_mhv( V_mhv ,V_tags ,V%Ntot ,V_h )
        if (V_mhv.eq.0) cycle
        U%Q(   1) = T%shift_sqr(cyc(1),cyc(N), zz)
        V%Q(jj+1) = T%shift_ang(cyc(N),cyc(1),-zz)
        U%Q(2:ii) = T%Q(cyc(2:ii))
        U%Q(ii+1) = Khat%minus()
        V%Q(1) = Khat
        V%Q(2:jj) = T%Q(cyc(ii+1:N-1))
        U%Noff = 0
        V%Noff = 0
        Urslt = amp_0( U_mhv ,U_tags ,U_h ,U )
        Vrslt = amp_0( V_mhv ,V_tags ,V_h ,V )
        rslt = rslt + Urslt*Vrslt/Ksqr
      enddo
    enddo
  endif
  end associate
  end function


  recursive function amp_1( mhv_type ,tags ,helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: mhv_type ,tags(2) ,helicity(T%Ntot)
  complex(fltknd) :: rslt ,Ksqr,zz,Urslt,Vrslt
  type(slashed_type) :: Kslash
  type(qomentum_list_type) :: U,V
  type(qomentum_type) :: ee,Khat 
  complex(fltknd) :: xx
  logical :: minusPlus
  integer :: cyc(NsizeProc),jj,ii,ihel
  integer :: U_mhv,U_tags(2),U_h(NsizeProc)
  integer :: V_mhv,V_tags(2),V_h(NsizeProc)
  associate( N=>T%Ntot )
  if     (mhv_type.eq. 0) then ; rslt = 0
  elseif (mhv_type.eq. 1) then ; rslt = amp_mostly_plus( tags ,T )
  elseif (mhv_type.eq.-1) then ; rslt = amp_mostly_minus( tags ,T )
  elseif (  T%Ntot.eq. 5) then ; rslt = amp_5_1( helicity ,T )
  else
    rslt = 0
! Rotate until leg 1 is off-shell.
    cyc(1:N) = unity(1:N)
    do while (helicity(cyc(1)).ne.0)
      jj = cyc(N)
      cyc(2:N) = cyc(1:N-1)
      cyc(1) = jj
    enddo
    minusPlus = (helicity(cyc(N)).eq.1)
    if (minusPlus) then ; ee = T%transvec( cyc(1) ,cyc(N) )
    else                ; ee = T%transvec( cyc(N) ,cyc(1) )
    endif
! Sum over partitions of external gluons
    do ii=2,N-2
      jj = N-ii
      if (ii.lt.jj) then
        Kslash = T%sumk(cyc(1:ii))
      else
        Kslash = zeroSlash-T%sumk(cyc(ii+1:N))
      endif
      Ksqr = square(Kslash)
      zz =-Ksqr/twodot(Kslash,ee%p)
      call Khat%fill( zz*ee%p + Kslash )
! Sum over internal helicities
      do ihel=-1,1,2
        U%Ntot = ii+1
        U_h(1:ii) = helicity(cyc(1:ii))
        U_h(ii+1) = ihel
        call determine_mhv( U_mhv ,U_tags ,U%Ntot ,U_h )
        if (U_mhv.eq.0) cycle
        V%Ntot = jj+1
        V_h(1) =-ihel
        V_h(2:jj+1) = helicity(cyc(ii+1:N))
        call determine_mhv( V_mhv ,V_tags ,V%Ntot ,V_h )
        if (V_mhv.eq.0) cycle
        if (minusPlus) then
          U%Q(   1) = T%shift_kapp(cyc(1),cyc(N), zz)
          V%Q(jj+1) = T%shift_ang( cyc(N),cyc(1),-zz)
        else
          U%Q(   1) = T%shift_kstr(cyc(1),cyc(N), zz)
          V%Q(jj+1) = T%shift_sqr( cyc(N),cyc(1),-zz)
        endif
        U%Q(2:ii) = T%Q(cyc(2:ii))
        U%Q(ii+1) = Khat%minus()
        V%Q(1) = Khat
        V%Q(2:jj) = T%Q(cyc(ii+1:N-1))
        !write(*,*) 'kapp',U%Q(1)%kapp !DEBUG
        !write(*,*) U_h(1:U%Ntot) ,'|  ',V_h(1:V%Ntot) ,'|  ',U_mhv!DEBUG
        U%Noff = 1
        V%Noff = 0
        Urslt = amp_1( U_mhv ,U_tags ,U_h ,U )
        Vrslt = amp_0( V_mhv ,V_tags ,V_h ,V )
        !write(*,'(99e16.8)') Urslt,Vrslt,Ksqr !DEBUG
        rslt = rslt + Urslt*Vrslt/Ksqr
      enddo
    enddo
!
    U%Ntot = N
    U_h(2:N) = helicity(cyc(2:N))
    if (minusPlus) then
      U_h(1) = 1
      call determine_mhv( U_mhv ,U_tags ,U%Ntot ,U_h )
      if (U_mhv.ne.0) then
        U%Q(1:N) = T%Q(cyc(1:N))
        call U%cterm( 1,N ,xx )
        U%Noff = 0
        Urslt = amp_0( U_mhv ,U_tags ,U_h ,U )
        rslt = rslt + Urslt/( xx*T%Q(cyc(1))%kapp )
        !write(*,'(a6,2e16.8,99i3)') 'amp1 C',Urslt,U_h(1:U%Ntot) !DEBUG
      endif
    else
      U_h(1) =-1
      call determine_mhv( U_mhv ,U_tags ,U%Ntot ,U_h )
      if (U_mhv.ne.0) then
        U%Q(1:N) = T%Q(cyc(1:N))
        call U%dterm( N,1 ,xx )
        U%Noff = 0
        Urslt = amp_0( U_mhv ,U_tags ,U_h ,U )
        rslt = rslt + Urslt/( xx*T%Q(cyc(1))%kstr )
        !write(*,'(a6,2e16.8,99i3)') 'amp1 D',Urslt,U_h(1:U%Ntot) !DEBUG
      endif
    endif
!
  endif
  end associate
  end function


  function amp_2( mhv_type ,tags ,helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: mhv_type ,tags(2) ,helicity(T%Ntot)
  complex(fltknd) :: rslt ,Ksqr,zz
  type(slashed_type) :: Kslash
  type(qomentum_list_type) :: U,V
  type(qomentum_type) :: ee,Khat 
  complex(fltknd) :: xx,Urslt,Vrslt
  integer :: cyc(NsizeProc),jj,ii,ihel
  integer :: U_mhv,U_tags(2),U_h(NsizeProc)
  integer :: V_mhv,V_tags(2),V_h(NsizeProc)
  associate( N=>T%Ntot )
  if     (mhv_type.eq. 0) then ; rslt = 0
  elseif (mhv_type.eq. 1) then ; rslt = amp_mostly_plus( tags ,T )
  elseif (mhv_type.eq.-1) then ; rslt = amp_mostly_minus( tags ,T )
  elseif (  T%Ntot.eq. 4) then ; rslt = amp_4_2( helicity ,T )
  elseif (  T%Ntot.eq. 5) then ; rslt = amp_5_2( helicity ,T )
  else
    rslt = 0
! Rotate until leg 1 and N are off-shell.
    cyc(1:N) = unity(1:N)
    do while (helicity(cyc(1)).ne.0.or.helicity(cyc(N)).ne.0)
      jj = cyc(1)
      cyc(1:N-1) = cyc(2:N)
      cyc(N) = jj
    enddo
    ee = T%transvec( cyc(1) ,cyc(N) )
! Sum over partitions of external gluons
    do ii=2,N-2
      jj = N-ii
      if (ii.lt.jj) then
        Kslash = T%sumk(cyc(1:ii))
      else
        Kslash = zeroSlash-T%sumk(cyc(ii+1:N))
      endif
      Ksqr = square(Kslash)
      zz =-Ksqr/T%ang(cyc(1),Kslash,cyc(N))
      call Khat%fill( zz*ee%p + Kslash )
! Sum over internal helicities
      do ihel=-1,1,2
        U%Ntot = ii+1
        U_h(1:ii) = helicity(cyc(1:ii))
        U_h(ii+1) = ihel
        call determine_mhv( U_mhv ,U_tags ,U%Ntot ,U_h )
        if (U_mhv.eq.0) cycle
        V%Ntot = jj+1
        V_h(1) =-ihel
        V_h(2:jj+1) = helicity(cyc(ii+1:N))
        call determine_mhv( V_mhv ,V_tags ,V%Ntot ,V_h )
        if (V_mhv.eq.0) cycle
        U%Q(   1) = T%shift_kapp(cyc(1),cyc(N), zz)
        U%Q(2:ii) = T%Q(cyc(2:ii))
        U%Q(ii+1) = Khat%minus()
        V%Q(   1) = Khat
        V%Q(2:jj) = T%Q(cyc(ii+1:N-1))
        V%Q(jj+1) = T%shift_kstr(cyc(N),cyc(1),-zz)
        U%Noff = 1
        V%Noff = 1
        Urslt = amp_1( U_mhv ,U_tags ,U_h ,U )
        Vrslt = amp_1( V_mhv ,V_tags ,V_h ,V )
        rslt = rslt + Urslt*Vrslt/Ksqr
        !write(*,'(a6,2e16.8,99i3)') 'amp2 U',Urslt,U_h(1:U%Ntot) !DEBUG
        !write(*,'(a6,2e16.8,99i3)') 'amp2 V',Vrslt,V_h(1:V%Ntot) !DEBUG
      enddo
    enddo
!
    U%Ntot = N
    U_h(2:N-1) = helicity(cyc(2:N-1))
    U_h(1) = 1
    U_h(N) = helicity(cyc(N))
    call determine_mhv( U_mhv ,U_tags ,U%Ntot ,U_h )
    if (U_mhv.ne.0) then
      U%Q(1:N) = T%Q(cyc(1:N))
      call U%cterm( 1,N ,xx )
      U%Noff = 1
      Urslt = amp_1( U_mhv ,U_tags ,U_h ,U )
      rslt = rslt + Urslt/( xx*T%Q(cyc(1))%kapp )
      !write(*,'(a6,2e16.8,99i3)') 'amp2 C',Urslt,U_h(1:U%Ntot) !DEBUG
    endif
    U_h(1) = helicity(cyc(1))
    U_h(N) =-1
    call determine_mhv( U_mhv ,U_tags ,U%Ntot ,U_h )
    if (U_mhv.ne.0) then
      U%Q(1:N) = T%Q(cyc(1:N))
      call U%dterm( 1,N ,xx )
      U%Noff = 1
      Urslt = amp_1( U_mhv ,U_tags ,U_h ,U )
      rslt = rslt + Urslt/( xx*T%Q(cyc(N))%kstr )
      !write(*,'(a6,2e16.8,99i3)') 'amp2 D',Urslt,U_h(1:U%Ntot) !DEBUG
    endif
!
  endif
  end associate
  end function


  function amp_mostly_plus( tag ,T ) result(rslt)
! Evaluate a mostly-plus helicity MHV amplitude
  type(qomentum_list_type),intent(in) :: T
  integer ,intent(in ) :: tag(2)
  complex(fltknd) :: rslt ,xx(3,3)
  integer :: ii
  associate( N=>T%Ntot )
  rslt = 0
  if (N.eq.3) then
    xx(1,2)=T%Q(1)%angL*T%Q(2)%Rang ;if(xx(1,2).eq.0)return ;xx(2,1)=-xx(1,2)
    xx(2,3)=T%Q(2)%angL*T%Q(3)%Rang ;if(xx(2,3).eq.0)return ;xx(3,2)=-xx(2,3)
    xx(3,1)=T%Q(3)%angL*T%Q(1)%Rang ;if(xx(3,1).eq.0)return ;xx(1,3)=-xx(3,1)
    rslt = xx(tag(1),tag(2))**4 / (xx(1,2)*xx(2,3)*xx(3,1))
  else
    xx(1,1)=T%Q(tag(1))%angL*T%Q(tag(2))%Rang ;if(xx(1,1).eq.0)return
    xx(2,2)=T%Q(N)%angL*T%Q(1)%Rang           ;if(xx(2,2).eq.0)return
    xx(1,1)= xx(1,1)**4/xx(2,2)
    do ii=1,N-1
      xx(2,2)=T%Q(ii)%angL*T%Q(ii+1)%Rang ;if(xx(2,2).eq.0)return 
      xx(1,1)= xx(1,1)/xx(2,2)
    enddo
    rslt = xx(1,1)
  endif
  end associate
  do ii=1,T%Noff
    rslt = rslt/T%Q(tag(ii))%kstr
  enddo
  end function


  function amp_mostly_minus( tag ,T ) result(rslt)
! Evaluate a mostly-minus helicity MHV amplitude
  type(qomentum_list_type),intent(in) :: T
  integer ,intent(in ) :: tag(2)
  complex(fltknd) :: rslt ,xx(3,3)
  integer :: ii
  associate( N=>T%Ntot )
  rslt = 0
  if (N.eq.3) then
    xx(1,2)=T%Q(1)%sqrL*T%Q(2)%Rsqr ;if(xx(1,2).eq.0)return ;xx(2,1)=-xx(1,2)
    xx(2,3)=T%Q(2)%sqrL*T%Q(3)%Rsqr ;if(xx(2,3).eq.0)return ;xx(3,2)=-xx(2,3)
    xx(3,1)=T%Q(3)%sqrL*T%Q(1)%Rsqr ;if(xx(3,1).eq.0)return ;xx(1,3)=-xx(3,1)
    rslt = xx(tag(1),tag(2))**4 / (xx(1,3)*xx(3,2)*xx(2,1))
  else
    xx(1,1)=T%Q(tag(1))%sqrL*T%Q(tag(2))%Rsqr ;if(xx(1,1).eq.0)return 
    xx(2,2)=T%Q(1)%sqrL*T%Q(N)%Rsqr           ;if(xx(2,2).eq.0)return
    xx(1,1)=xx(1,1)**4/xx(2,2)
    do ii=1,N-1
      xx(2,2)=T%Q(ii+1)%sqrL*T%Q(ii)%Rsqr ;if(xx(2,2).eq.0)return
      xx(1,1)=xx(1,1)/xx(2,2)
    enddo
    rslt = xx(1,1)
  endif
  end associate
  do ii=1,T%Noff
    rslt = rslt/T%Q(tag(ii))%kapp
  enddo
  end function
  

  function amp_4_2( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(4)
  complex(fltknd) :: rslt
  integer :: o1,o2,i3,i4,hh
  o1=1 ;o2=2 ;i3=3 ;i4=4
  do while (helicity(o1).ne.0.or.helicity(o2).ne.0)
    hh=i4 ;i4=i3 ;i3=o2 ;o2=o1 ;o1=hh
  enddo
  if (helicity(i3).eq.1) then
    if (helicity(i4).eq.1) then !**++
      rslt = 1/(T%Q(o1)%kstr*T%Q(o2)%kstr) * T%ang(o1,o2)**3 / (T%ang(o2,i3)*T%ang(i3,i4)*T%ang(i4,o1))
    else !**+-
      rslt = 1/(T%Q(o1)%kstr*T%Q(o2)%kapp) * T%sqr(o2,i3)**3 * T%ang(i4,o1)**3 /(T%ang(o1,o2,i3)*T%ang(i4,o1,o2)*T%sinv(i4,o1)) &
           + 1/T%Q(o2)%kstr  * T%sqr(o1,i3)**4 * T%ang(o1,o2)**3 &
             /(T%ang(o1,o2,i3)*T%ang(o2,[i3,i4],o1)*T%sqr(o1,i4)*T%sqr(i4,i3)*T%ang(o1,[i3,i4],o1)) &
           + 1/T%Q(o1)%kapp * T%ang(o2,i4)**4 * T%sqr(o1,o2)**3 &
             /(T%ang(o2,[o1,o2],o1)*T%ang(i4,o1,o2)*T%ang(o2,o1,o2)*T%ang(o2,i3)*T%ang(i3,i4))
    endif
  else
    if (helicity(i4).eq.1) then !**-+
      rslt = 1/(T%Q(o1)%kapp*T%Q(o2)%kstr) * T%ang(o2,i3)**3 * T%sqr(i4,o1)**3 /(T%sqr(o1,o2,i3)*T%sqr(i4,o1,o2)*T%sinv(i4,o1)) &
           - 1/T%Q(o2)%kapp  * T%ang(o1,i3)**4 * T%sqr(o1,o2)**3 &
             /(T%sqr(o1,o2,i3)*T%sqr(o2,[i3,i4],o1)*T%ang(o1,i4)*T%ang(i4,i3)*T%sqr(o1,[i3,i4],o1)) &
           - 1/T%Q(o1)%kstr * T%sqr(o2,i4)**4 * T%ang(o1,o2)**3 &
             /(T%sqr(o2,[o1,o2],o1)*T%sqr(i4,o1,o2)*T%sqr(o2,o1,o2)*T%sqr(o2,i3)*T%sqr(i3,i4))
    else !**--
      rslt = 1/(T%Q(o1)%kapp*T%Q(o2)%kapp) * T%sqr(o2,o1)**3 / (T%sqr(i4,i3)*T%sqr(i3,o2)*T%sqr(o1,i4))
    endif
  endif
  end function


  function amp_5_1( helicity ,T) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(5)
  complex(fltknd) :: rslt
  integer :: o1,i2,i3,i4,i5,hh
  o1=1 ;i2=2 ;i3=3 ;i4=4 ;i5=5
  do while (helicity(o1).ne.0)
    hh=i5 ;i5=i4 ;i4=i3 ;i3=i2 ;i2=o1 ;o1=hh
  enddo
  if (helicity(i2).eq.1) then
    if (helicity(i3).eq.1) then
      if (helicity(i4).eq.1) then
        if (helicity(i5).eq.1) then !*++++
          rslt = 0
        else !*+++-
          rslt = T%ang(i5,o1)**3/(T%ang(o1,i2)*T%ang(i2,i3)*T%ang(i3,i4)*T%ang(i4,i5))
          rslt = rslt/T%Q(o1)%kstr
        endif
      else
        if (helicity(i5).eq.1) then !*++-+
          rslt = T%ang(o1,i4)**4/(T%ang(o1,i2)*T%ang(i2,i3)*T%ang(i3,i4)*T%ang(i4,i5)*T%ang(i5,o1))
          rslt = rslt/T%Q(o1)%kstr
        else !*++--
          rslt = 1/T%Q(o1)%kstr * T%ang(o1,i5)**3 * T%sqr(i2,i3)**3 &
                 /(T%ang(i5,o1,i2)*T%ang(o1,[o1,i5],i4)*T%sqr(i3,i4)*T%sinv(o1,i5)) &
               + T%ang(o1,[i3,i2],o1)**3 &
                 /(T%ang(i3,[o1,i2],o1)*T%ang(o1,[o1,i5],i4)*T%sqr(i4,i5)*T%sqr(i5,o1)*T%ang(i3,i2)*T%ang(i2,o1)) &
               + 1/T%Q(o1)%kapp * T%ang(i5,i4)**3 *T%sqr(o1,i2)**3 &
                 /(T%ang(i5,o1,i2)*T%ang(i4,i3)*T%ang(i3,[i2,o1],o1)*T%sinv(o1,i2))
          rslt =-rslt
        endif
      endif
    else
      if (helicity(i4).eq.1) then
        if (helicity(i5).eq.1) then !*+-++
          rslt = T%ang(o1,i3)**4/(T%ang(o1,i2)*T%ang(i2,i3)*T%ang(i3,i4)*T%ang(i4,i5)*T%ang(i5,o1))
          rslt = rslt/T%Q(o1)%kstr
        else !*+-+-
          rslt = 1/T%Q(o1)%kstr * T%ang(o1,i5)**3 * T%sqr(i2,i4)**4 &
                 /(T%ang(i5,o1,i2)*T%ang(o1,[o1,i5],i4)*T%sqr(i2,i3)*T%sqr(i3,i4)*T%sinv(o1,i5)) &
               + T%ang(o1,i3)**4 *T%sqr(i4,o1)**4 &
                 /(T%ang(i3,[o1,i2],o1)*T%ang(o1,[o1,i5],i4)*T%ang(o1,[i3,i2],o1) &
                  *T%sqr(i4,i5)*T%sqr(i5,o1)*T%ang(i3,i2)*T%ang(i2,o1)) &
               + 1/T%Q(o1)%kapp * T%ang(i5,i3)**4 *T%sqr(o1,i2)**3 &
                 /(T%ang(i5,o1,i2)*T%ang(i5,i4)*T%ang(i4,i3)*T%ang(i3,[i2,o1],o1)*T%sinv(o1,i2))
          rslt =-rslt
        endif
      else
        if (helicity(i5).eq.1) then !*+--+
          rslt = 1/T%Q(o1)%kstr * T%ang(o1,[o1,i2],i5)**4 &
                 /(T%ang(o1,i2)*T%sqr(i3,i4)*T%sqr(i4,i5)*T%ang(i2,o1,i5)*T%ang(o1,[o1,i2],i3)*T%sinv(o1,i2)) &
               + T%ang(i4,o1)**4 * T%sqr(o1,i2)**3 &
                 /(T%ang(i4,[o1,i5],o1)*T%ang(o1,[i4,i5],o1)*T%ang(o1,[i4,i5],i3)*T%sqr(i3,i2)*T%ang(i4,i5)*T%ang(i5,o1) ) &
               + 1/T%Q(o1)%kapp * T%ang(i3,i4)**3 *T%sqr(o1,i5)**3 &
                 /(T%ang(i2,o1,i5)*T%ang(i2,i3)*T%ang(i4,[i5,o1],o1)*T%sinv(o1,i5))
        else !*+---
          rslt = T%sqr(i2,o1)**3/(T%sqr(i5,i4)*T%sqr(i4,i3)*T%sqr(i3,i2)*T%sqr(o1,i5))
          rslt = rslt/T%Q(o1)%kapp
        endif
      endif
    endif
  else
    if (helicity(i3).eq.1) then
      if (helicity(i4).eq.1) then
        if (helicity(i5).eq.1) then !*-+++
          rslt = T%ang(o1,i2)**3/(T%ang(i2,i3)*T%ang(i3,i4)*T%ang(i4,i5)*T%ang(i5,o1))
          rslt = rslt/T%Q(o1)%kstr
        else !*-++-
          rslt =-1/T%Q(o1)%kapp * T%sqr(o1,[o1,i2],i5)**4 &
                 /(T%sqr(o1,i2)*T%ang(i3,i4)*T%ang(i4,i5)*T%sqr(i2,o1,i5)*T%sqr(o1,[o1,i2],i3)*T%sinv(o1,i2)) &
               + T%sqr(i4,o1)**4 * T%ang(o1,i2)**3 &
                 /(T%sqr(i4,[o1,i5],o1)*T%sqr(o1,[i4,i5],o1)*T%sqr(o1,[i4,i5],i3)*T%ang(i3,i2)*T%sqr(i4,i5)*T%sqr(i5,o1)) &
               - 1/T%Q(o1)%kstr * T%sqr(i3,i4)**3 *T%ang(o1,i5)**3 &
                 /(T%sqr(i2,o1,i5)*T%sqr(i2,i3)*T%sqr(i4,[i5,o1],o1)*T%sinv(o1,i5))
        endif
      else
        if (helicity(i5).eq.1) then !*-+-+
          rslt = 1/T%Q(o1)%kstr * T%ang(o1,i2)**3 * T%sqr(i5,i3)**4 &
                 /(T%ang(i2,o1,i5)*T%ang(o1,[o1,i2],i3)*T%sqr(i5,i4)*T%sqr(i4,i3)*T%sinv(o1,i2)) &
               + T%ang(o1,i4)**4 *T%sqr(i3,o1)**4 &
                 /(T%ang(i4,[o1,i5],o1)*T%ang(o1,[o1,i2],i3)*T%ang(o1,[i4,i5],o1) &
                  *T%sqr(i3,i2)*T%sqr(i2,o1)*T%ang(i4,i5)*T%ang(i5,o1)) &
               + 1/T%Q(o1)%kapp * T%ang(i2,i4)**4 *T%sqr(o1,i5)**3 &
                 /(T%ang(i2,o1,i5)*T%ang(i2,i3)*T%ang(i3,i4)*T%ang(i4,[i5,o1],o1)*T%sinv(o1,i5))
        else !*-+--
          rslt = T%sqr(i3,o1)**4/(T%sqr(i5,i4)*T%sqr(i4,i3)*T%sqr(i3,i2)*T%sqr(i2,o1)*T%sqr(o1,i5))
          rslt = rslt/T%Q(o1)%kapp
        endif
      endif
    else
      if (helicity(i4).eq.1) then
        if (helicity(i5).eq.1) then !*--++
          rslt = 1/T%Q(o1)%kstr * T%ang(o1,i2)**3 * T%sqr(i5,i4)**3 &
                 /(T%ang(i2,o1,i5)*T%ang(o1,[o1,i2],i3)*T%sqr(i4,i3)*T%sinv(o1,i2)) &
               + T%ang(o1,[i4,i5],o1)**3 &
                 /(T%ang(i4,[o1,i5],o1)*T%ang(o1,[o1,i2],i3)*T%sqr(i3,i2)*T%sqr(i2,o1)*T%ang(i4,i5)*T%ang(i5,o1)) &
               + 1/T%Q(o1)%kapp * T%ang(i2,i3)**3 *T%sqr(o1,i5)**3 &
                 /(T%ang(i2,o1,i5)*T%ang(i3,i4)*T%ang(i4,[i5,o1],o1)*T%sinv(o1,i5))
        else !*--+-
          rslt = T%sqr(i4,o1)**4/(T%sqr(i5,i4)*T%sqr(i4,i3)*T%sqr(i3,i2)*T%sqr(i2,o1)*T%sqr(o1,i5))
          rslt = rslt/T%Q(o1)%kapp
        endif
      else
        if (helicity(i5).eq.1) then !*---+
          rslt = T%sqr(o1,i5)**4/(T%sqr(i5,i4)*T%sqr(i4,i3)*T%sqr(i3,i2)*T%sqr(i2,o1))
          rslt = rslt/T%Q(o1)%kapp
        else !*----
          rslt = 0
        endif
      endif
    endif
  endif
  end function


  function amp_5_2( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(5)
  complex(fltknd) :: rslt ,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10
  complex(fltknd) :: h11,h12,h13,h14,h15,h16,h17,h18,h19
  integer :: o1,o2,i3,i4,i5,hh
  o1=1 ;o2=2 ;i3=3 ;i4=4 ;i5=5
  do while (helicity(o1).ne.0.or.helicity(o2).ne.0)
    hh=i5 ;i5=i4 ;i4=i3 ;i3=o2 ;o2=o1 ;o1=hh
  enddo
  if (helicity(i3).eq.1) then
    if (helicity(i4).eq.1) then
      if (helicity(i5).eq.1) then !**+++
        rslt = 1/(T%Q(o1)%kstr*T%Q(o2)%kstr) * T%ang(o1,o2)**3 /(T%ang(o2,i3)*T%ang(i3,i4)*T%ang(i4,i5)*T%ang(i5,o1))
      else !**++-
        h1=T%ang(o1,o2)**3
        h2=T%ang(i5,o1)**3
        h3=T%ang(o2,[i4,i3],[o1,i5],o1)
        h4=T%ang(o2,i3)*T%ang(i3,i4)
        h5=T%ang(i4,i5)
        h6=T%sqr(i3,[o2,o1],o1)
        h7=T%sqr(o2,[o2,i3],i4)
        h8=T%sqr(o2,o1,i5)
        h9=T%sqr(o1,[o2,o1],o2)
        h10=T%sqr(i3,o2,o1)
        h11=T%sqr(o1,o2,o1)*T%ang(o2,i4)+T%sqr(o1,i3)*T%ang(i3,i4)*T%ang(o2,o1)
        rslt = 1/(T%Q(o2)%kstr*T%Q(o1)%kstr) * h1 * T%sqr(i4,i3)**3 / (T%sqr(i5,[o2,o1],o2)*h6*T%sqr(i5,i4)*T%sinv(o2,o1)) &
             - 1/T%Q(o1)%kstr * h2 * T%sqr(o2,[i4,i3],o2)**3 / (h7*h8*h3*h4*T%sinv(o1,i5)) &
             - 1/T%Q(o1)%kapp * T%ang(o2,i5)**4 *T%sqr(o2,o1)**3 / (h8*h9*T%sqr(o2,o1,o2)*h4*h5) &
             + 1/T%Q(o2)%kstr * T%sqr(o1,[i4,i3],o2)**4 * h1  / (h9*h3*h11*T%sqr(i5,[i4,i3],o2)*T%sqr(o1,i5)*h4) &
             + 1/T%Q(o2)%kstr * h1 *T%sqr(o1,i3)**4 * h2 / (h10*T%sqr(o1,o2,o1)*h6*T%sqr(o1,[o2,i3],o1)*h5*h11)&
             + 1/(T%Q(o2)%kapp*T%Q(o1)%kstr) * T%sqr(i3,o2)**3 * h2  / (h7*h10*h5*T%sinv(o2,i3))
      endif
    else
      if (helicity(i5).eq.1) then !**+-+
        h1 = T%ang(o2,o1)
        h2 = h1**3
        h3 = T%ang(i4,o1)
        h4 = h3**4
        h5 = T%ang(i3,o2)
        h6 = T%ang(i4,o2)**4
        h7 = T%ang(i4,i3)
        h8 = T%ang(i5,o1)
        h9 = T%ang(i4,i5)
        h10 = T%sqr(i5,o1)
        h11 = T%sqr(i5,o2)
        h12 = T%ang(o2,o1,o2)
        h13 = T%ang(o1,o2,i3)
        h14 = T%ang(o2,o1,i5)
        h15 = T%ang(o2,[o1,o2],o1)
        h16 = T%ang(i4,[o1,i5],o1)
        h17 = T%ang(i4,[o1,i5],o2)
        h18 = T%ang(o2,[o1,o2],i5)
        h19 = T%ang(o1,[i4,i5],i3)
        rslt =-1/(T%Q(o1)%kstr*T%Q(o2)%kstr) * T%sqr(i5,i3)**4 * h2 / (h19*h18*T%sqr(i4,i3)*T%sqr(i5,i4)*T%sinv(o1,o2)) &
             - 1/(T%Q(o1)%kstr*T%Q(o2)%kapp) * h4 * T%sqr(o2,i3)**3 / (h17*h13*h8*h9*T%sinv(o2,i3))&
             - h6 *T%sqr(o1,o2)**3 *h4 / (h16*h17*h5*h7*h8*h9*( h15*h3 + h10*h9*h1 )*( h12*h3 + h11*h9*h1 )) &
             - 1/T%Q(o2)%kstr * T%sqr(i3,o1)**4 * h2 * h4 &
                / (h13*T%sqr(o1,o2,o1)*T%sqr(o1,[i4,i5],o1)*h19*h8*h9*( h15*h3 + h10*h9*h1 )) &
             - 1/T%Q(o1)%kstr * h6 * h11**4 * h2  / (h18*h14*h12*T%sqr(o2,[o1,i5],o2)*h7 *h5*( h11*h1*h9+h12*h3 )) &
             + 1/(T%Q(o1)%kapp*T%Q(o2)%kstr) * h6 * h10**3 / (h16*h14*h7*h5*T%sinv(o1,i5))
        rslt =-rslt
      else !**+--
        h1=T%sqr(o2,o1)**3
        h2=T%sqr(i3,o2)**3
        h3=T%sqr(o1,[i4,i5],[o2,i3],o2)
        h4=T%sqr(o1,i5)*T%sqr(i5,i4)
        h5=T%sqr(i4,i3)
        h6=T%ang(i5,[o1,o2],o2)
        h7=T%ang(o1,[o1,i5],i4)
        h8=T%ang(o1,o2,i3)
        h9=T%ang(o2,[o1,o2],o1)
        h10=T%ang(i5,o1,o2)
        h11=T%ang(o2,o1,o2)*T%sqr(o1,i4)+T%ang(o2,i5)*T%sqr(i5,i4)*T%sqr(o1,o2)
        rslt =-1/(T%Q(o1)%kapp*T%Q(o2)%kapp) * h1 * T%ang(i4,i5)**3 / (T%ang(i3,[o1,o2],o1)*h6*T%ang(i3,i4)*T%sinv(o1,o2)) &
             - 1/T%Q(o2)%kapp * h2 * T%ang(o1,[i4,i5],o1)**3 / (h7*h8*h3*h4*T%sinv(o2,i3)) &
             - 1/T%Q(o2)%kstr * T%sqr(o1,i3)**4 *T%ang(o1,o2)**3 / (h8*h9*T%ang(o1,o2,o1)*h4*h5) &
             + 1/T%Q(o1)%kapp * T%ang(o2,[i4,i5],o1)**4 * h1 / (h9*h3*h11*T%ang(i3,[i4,i5],o1)*T%ang(o2,i3)*h4) &
             + 1/T%Q(o1)%kapp * h1 *T%ang(o2,i5)**4 * h2 / (h10*T%ang(o2,o1,o2)*h6*T%ang(o2,[o1,i5],o2)*h5*h11)&
             - 1/(T%Q(o1)%kstr*T%Q(o2)%kapp) * T%ang(i5,o1)**3 * h2 / (h7*h10*h5*T%sinv(o1,i5))
        rslt =-rslt
      endif
    endif
  else
    if (helicity(i4).eq.1) then
      if (helicity(i5).eq.1) then !**-++
        h1=T%ang(o2,o1)**3
        h2=T%ang(i3,o2)**3
        h3=T%ang(o1,[i4,i5],[o2,i3],o2)
        h4=T%ang(o1,i5)*T%ang(i5,i4)
        h5=T%ang(i4,i3)
        h6=T%sqr(i5,[o1,o2],o2)
        h7=T%sqr(o1,[o1,i5],i4)
        h8=T%sqr(o1,o2,i3)
        h9=T%sqr(o2,[o1,o2],o1)
        h10=T%sqr(i5,o1,o2)
        h11=T%sqr(o2,o1,o2)*T%ang(o1,i4)+T%sqr(o2,i5)*T%ang(i5,i4)*T%ang(o1,o2)
        rslt = 1/(T%Q(o1)%kstr*T%Q(o2)%kstr) * h1 * T%sqr(i4,i5)**3 / (T%sqr(i3,[o1,o2],o1)*h6*T%sqr(i3,i4)*T%sinv(o1,o2)) &
             - 1/T%Q(o2)%kstr * h2 * T%sqr(o1,[i4,i5],o1)**3 / (h7*h8*h3*h4*T%sinv(o2,i3)) &
             - 1/T%Q(o2)%kapp * T%ang(o1,i3)**4 *T%sqr(o1,o2)**3 / (h8*h9*T%sqr(o1,o2,o1)*h4*h5) &
             + 1/T%Q(o1)%kstr * T%sqr(o2,[i4,i5],o1)**4 * h1  / (h9*h3*h11*T%sqr(i3,[i4,i5],o1)*T%sqr(o2,i3)*h4) &
             + 1/T%Q(o1)%kstr * h1 *T%sqr(o2,i5)**4 * h2 / (h10*T%sqr(o2,o1,o2)*h6*T%sqr(o2,[o1,i5],o2)*h5*h11)&
             + 1/(T%Q(o1)%kapp*T%Q(o2)%kstr) * T%sqr(i5,o1)**3 * h2  / (h7*h10*h5*T%sinv(o1,i5))
        rslt=-rslt
      else !**-+-
        h1 = T%sqr(o1,o2)
        h2 = h1**3
        h3 = T%sqr(o1,i4)
        h4 = h3**4
        h5 = T%sqr(o2,i3)
        h6 = T%sqr(o2,i4)**4
        h7 = T%sqr(i3,i4)
        h8 = T%sqr(o1,i5)
        h9 = T%sqr(i5,i4)
        h10 = T%ang(o1,i5)
        h11 = T%ang(o2,i5)
        h12 = T%ang(o2,o1,o2)
        h13 = T%ang(i3,o2,o1)
        h14 = T%ang(i5,o1,o2)
        h15 = T%ang(o1,[o1,o2],o2)
        h16 = T%ang(o1,[o1,i5],i4)
        h17 = T%ang(o2,[o1,i5],i4)
        h18 = T%ang(i5,[o1,o2],o2)
        h19 = T%ang(i3,[i4,i5],o1)
        rslt =-1/(T%Q(o1)%kapp*T%Q(o2)%kapp) * T%ang(i3,i5)**4 * h2 / (h19*h18*T%ang(i3,i4)*T%ang(i4,i5)*T%sinv(o1,o2)) &
             - 1/(T%Q(o1)%kapp*T%Q(o2)%kstr) * h4 * T%ang(i3,o2)**3 / (h17*h13*h8*h9*T%sinv(o2,i3))&
             - h6 *T%ang(o2,o1)**3 *h4 / (h16*h17*h5*h7*h8*h9*( h15*h3 + h10*h9*h1 )*( h12*h3 + h11*h9*h1 )) &
             - 1/T%Q(o2)%kapp * T%ang(o1,i3)**4 * h2 * h4 &
                / (h13*T%ang(o1,o2,o1)*T%ang(o1,[i4,i5],o1)*h19*h8*h9*( h15*h3 + h10*h9*h1 )) &
             - 1/T%Q(o1)%kapp * h6 * h11**4 * h2  / (h18*h14*h12*T%ang(o2,[o1,i5],o2)*h7 *h5*( h11*h1*h9+h12*h3 )) &
             + 1/(T%Q(o1)%kstr*T%Q(o2)%kapp) * h6 * h10**3 / (h16*h14*h7*h5*T%sinv(o1,i5))
        rslt =-rslt
      endif
    else
      if (helicity(i5).eq.1) then !**--+
        h1=T%sqr(o1,o2)**3
        h2=T%sqr(i5,o1)**3
        h3=T%sqr(o2,[i4,i3],[o1,i5],o1)
        h4=T%sqr(o2,i3)*T%sqr(i3,i4)
        h5=T%sqr(i4,i5)
        h6=T%ang(i3,[o2,o1],o1)
        h7=T%ang(o2,[o2,i3],i4)
        h8=T%ang(o2,o1,i5)
        h9=T%ang(o1,[o2,o1],o2)
        h10=T%ang(i3,o2,o1)
        h11=T%ang(o1,o2,o1)*T%sqr(o2,i4)+T%ang(o1,i3)*T%sqr(i3,i4)*T%sqr(o2,o1)
        rslt = 1/(T%Q(o2)%kapp*T%Q(o1)%kapp) * h1 * T%ang(i4,i3)**3 / (T%ang(i5,[o2,o1],o2)*h6*T%ang(i5,i4)*T%sinv(o2,o1)) &
             + 1/T%Q(o1)%kapp * h2 * T%ang(o2,[i4,i3],o2)**3 / ( h7*h8*h3*h4*T%sinv(o1,i5)) &
             + 1/T%Q(o1)%kstr * T%sqr(o2,i5)**4 *T%ang(o2,o1)**3 / (h8*h9*T%ang(o2,o1,o2)*h4*h5) &
             - 1/T%Q(o2)%kapp * T%ang(o1,[i4,i3],o2)**4 * h1 / (h9*h3*h11*T%ang(i5,[i4,i3],o2)*T%ang(o1,i5)*h4) &
             - 1/T%Q(o2)%kapp * h1 *T%ang(o1,i3)**4 * h2 / (h10*T%ang(o1,o2,o1)*h6*T%ang(o1,[o2,i3],o1)*h5*h11)&
             + 1/(T%Q(o2)%kstr*T%Q(o1)%kapp) * T%ang(i3,o2)**3 * h2 / (h7*h10*h5*T%sinv(o2,i3))
        rslt =-rslt
      else !**---
        rslt = 1/(T%Q(o1)%kapp*T%Q(o2)%kapp) * T%sqr(o2,o1)**3 /(T%sqr(i5,i4)*T%sqr(i4,i3)*T%sqr(i3,o2)*T%sqr(o1,i5))
      endif
    endif
  endif
  end function


  function colorSum( Nminus2 ,a ) result(rslt)
  integer,intent(in) :: Nminus2
  complex(fltknd),intent(in) :: a(factorial(Nminus2))
  complex(fltknd) :: rslt ,uu,vv
  integer :: i
  select case (Nminus2)
  case (1)
    rslt = conjg(a(1))*a(1)
  case (2)
    rslt = conjg(a(1))*( 2*a(1) + a(2) ) &
         + conjg(a(2))*( 2*a(2) + a(1) )
  case (3)
    associate(m=>mtx_5_sqr)
    rslt = 0
    do i=1,6
      rslt = rslt + conjg(a(i))*( 4*a(m(1,i)) + 2*sum(a(m(2:3,i))) + sum(a(m(4:5,i))) )
    enddo
    end associate
  case (4)
    associate(m=>mtx_6_sqr)
    uu = 0
    do i=1,24
      uu = uu + conjg(a(i))*( 8*a(m(1,i)) + 4*sum(a(m(2:4,i))) &
                             + 2*sum(a(m(5:9,i))) + sum(a(m(10:14,i))) )
    enddo
    vv = 0
    do i=1,24
      vv = vv + conjg(a(i))*sum(a(m(15:19,i)))
    enddo
    rslt = uu + 12*vv/Ncolor(2)
    end associate
  case (5)
    associate(m=>mtx_7_sqr)
    uu = 0
    do i=1,120
      uu = uu + conjg(a(i))*( 16*a(m(1,i)) + 8*sum(a(m(2:5,i))) &
                   + 4*sum(a(m(6:14,i))) + 2*sum(a(m(15:28,i))) + sum(a(m(29:42,i))) )
    enddo
    vv = 0
    do i=1,120
      vv = vv + conjg(a(i))*( 24*sum(a(m(43:52,i))) + 16*sum(a(m(53:57,i))) &
                 + 12*sum(a(m(58:77,i))) + 8*sum(a(m(78:84,i))) + 4*sum(a(m(85:94,i))) &
                 - 4*sum(a(m(95:100,i))) - 16*sum(a(m(101:101,i))) )
    enddo
    rslt = uu + vv/Ncolor(2)
    end associate
  end select
  rslt = rslt * 2 * Ncolor(Nminus2) * (Ncolor(2)-1)
  end function  


  subroutine fill_matrices_ng
  integer :: rUnit
  call check_file('amp4hef.tbl')
  open(newunit=rUnit,file=trim(get_path('amp4hef.tbl')),status='old')
  call read_matrix( 4 ,mtx_4_sqr ,'BEGIN 4-gluon square' )
  call read_matrix( 5 ,mtx_5_sqr ,'BEGIN 5-gluon square' )
  call read_matrix( 6 ,mtx_6_sqr ,'BEGIN 6-gluon square' )
  call read_matrix( 7 ,mtx_7_sqr ,'BEGIN 7-gluon square' )
  close(rUnit)
!
  contains
!
    subroutine read_matrix( nn ,matrix ,tag )
    integer,intent(in) :: nn
    integer,allocatable,intent(out) :: matrix(:,:)
    character(*),intent(in) :: tag
    character(144) :: line
    integer :: ii
    allocate(matrix(factorial(nn-2),factorial(nn-2)))
    do
      read(rUnit,'(A)') line
      if (line(1:len(tag)).eq.tag) exit
    enddo
    do ii=1,factorial(nn-2)
      read(rUnit,*) matrix(1:factorial(nn-2),ii)
    enddo
    end subroutine
!
  end subroutine


  subroutine print_color_matrix( Nminus2 ,writeUnit )
  integer,intent(in) :: Nminus2,writeUnit
  optional :: writeUnit
  integer :: i,wUnit,overall
  integer,allocatable :: mat(:,:)
  wUnit = 6
  if (present(writeUnit)) wUnit=writeUnit
  allocate(mat(factorial(Nminus2),factorial(Nminus2)))
  mat = 0
  overall = 2 * Ncolor(Nminus2) * (Ncolor(2)-1)
  select case (Nminus2)
  case (2)
    mat(1,1:2) = overall*[2,1]
    mat(2,1:2) = overall*[1,2]
  case (3)
    associate(m=>mtx_5_sqr)
    do i=1,6
      mat(i,m(1  ,i)) = overall*4
      mat(i,m(2:3,i)) = overall*2
      mat(i,m(4:5,i)) = overall
    enddo
    end associate
  case (4)
    associate(m=>mtx_6_sqr)
    do i=1,24
      mat(i,m( 1   ,i)) = overall*8
      mat(i,m( 2: 4,i)) = overall*4
      mat(i,m( 5: 9,i)) = overall*2
      mat(i,m(10:14,i)) = overall
    enddo
    overall = (overall*12)/Ncolor(2)
    do i=1,24
      mat(i,m(15:19,i)) = mat(i,m(15:19,i)) + overall
    enddo
    end associate
  case (5)
    associate(m=>mtx_7_sqr)
    do i=1,120
      mat(i,m( 1   ,i)) = overall*16
      mat(i,m( 2: 5,i)) = overall*8
      mat(i,m( 6:14,i)) = overall*4
      mat(i,m(15:28,i)) = overall*2
      mat(i,m(29:42,i)) = overall
    enddo
    overall = overall/Ncolor(2)
    do i=1,120
      mat(i,m(43: 52,i)) = mat(i,m(43: 52,i)) + overall*24
      mat(i,m(53: 57,i)) = mat(i,m(53: 57,i)) + overall*16
      mat(i,m(58: 77,i)) = mat(i,m(58: 77,i)) + overall*12
      mat(i,m(78: 84,i)) = mat(i,m(78: 84,i)) + overall*8
      mat(i,m(85: 94,i)) = mat(i,m(85: 94,i)) + overall*4
      mat(i,m(95:100,i)) = mat(i,m(95:100,i)) - overall*4
      mat(i,m(101   ,i)) = mat(i,m(101   ,i)) - overall*16
    enddo
    end associate
  end select
  do i=1,factorial(Nminus2)
    write(wUnit,'(999i7)') mat(i,:)
  enddo
  end subroutine


end module



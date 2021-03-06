module amp4hef_qq
  use amp4hef_io
  use amp4hef_aux
  use amp4hef_qomentum

  implicit none
  private
  public :: fill_matrices_qq,matrix_element_qq ,amplitude_qq ,all_amplitudes_qq

  integer,parameter :: gluon=0 ,quark=1 ,antiq=-1

  integer,allocatable,save :: mtx_4_sqr(:,:)

contains


  function matrix_element_qq( Tin ) result(rslt)
  class(qomentum_list_type),intent(in) :: Tin
  real(fltknd) :: rslt
  complex(fltknd) :: amp(120)
  integer :: ii,jj,NhelSum,Nminus2
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )
  NhelSum = Ntot-1
  if (Tin%offshell(1).eq.Tin%flavor(1,gluon)) NhelSum = NhelSum-1
  if (Tin%offshell(2).eq.Tin%flavor(1,gluon)) NhelSum = NhelSum-1
  if (Tin%offshell(2).eq.Tin%flavor(2,gluon)) NhelSum = NhelSum-1
!
  Nminus2 = Ntot-2
  rslt = 0
  do ii=1,twoPower(NhelSum-1)
    do jj=1,factorial(Nminus2)
      amp(jj) = amplitude_qq( Tin ,helTable(1:NhelSum,ii) ,perTable(1:Nminus2,jj) )
!      write(*,'(2e16.8,99i3)') amp(jj),helTable(1:NhelSum,ii),perTable(1:Nminus2,jj) !DEBUG
    enddo
!    write(*,*) !DEBUG
    rslt = rslt + colorSum( Nminus2 ,amp(1:factorial(Nminus2)) )
  enddo
  rslt = rslt*2
  do ii=1,Noff
    rslt = rslt * (Tin%Q(Tin%offshell(ii))%kapp*Tin%Q(Tin%offshell(ii))%kstr)
  enddo
  end associate
  end function 


  subroutine all_amplitudes_qq( Tin ,NhelConf ,Nperm ,amplitude ,factor )
  class(qomentum_list_type),intent(in) :: Tin
  integer,intent(out) :: NhelConf,Nperm
  complex(fltknd),intent(out) :: amplitude(:,:),factor
  integer :: ii,jj,NhelSum,Nminus2
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff )
! The following is valid only if there is at most one off-shell parton.
  NhelSum = Ntot-1
  if (Tin%offshell(1).eq.Tin%flavor(1,gluon)) NhelSum = NhelSum-1
!
  Nminus2 = Ntot-2
  NhelConf = twoPower(NhelSum-1)
  Nperm = factorial(Nminus2)
  do ii=1,NhelConf
  do jj=1,Nperm
    amplitude(ii,jj) = amplitude_qq( Tin ,helTable(1:NhelSum,ii) ,perTable(1:Nminus2,jj) )
  enddo
  enddo
  factor = 2 ! only half of the helicity configurations is returned
  do ii=1,Noff
    factor = factor * (Tin%Q(Tin%offshell(ii))%kapp*Tin%Q(Tin%offshell(ii))%kstr)
  enddo
  end associate
  end subroutine


  function amplitude_qq( Tin ,helicity ,perm ) result(rslt)
! The first entries of helicity refer to the on-shell gluons. The
! next helicity to the anti-quark. Off-shell (anti)-quarks have helicity,
! and the quark always gets the opposite helicity of the anti-quark.
! If 1(2) gluon(s) is(are) off-shell, it must be the first (2) gluon(s).
! perm(1:N-2) permutes the gluons.
  class(qomentum_list_type),intent(in) :: Tin
  integer,intent(in) :: helicity(:),perm(:)
  complex(fltknd) :: rslt
  type(qomentum_list_type) :: T
  integer :: hel(-1:NsizeProc),per(NsizeProc) 
  associate( Ntot=>Tin%Ntot ,Noff=>Tin%Noff ,offshell=>Tin%offshell )

  per(1:Ntot-2) = perm(1:Ntot-2)
  per(Ntot-1:Ntot) = [Ntot-1,Ntot]
  hel(-1:0) = 0
  hel(1:size(helicity)) = helicity

  T%Ntot = Ntot
  T%Q(1:Ntot-2) = Tin%Q(Tin%flavor(per(1:Ntot-2),gluon))
  T%Q(Ntot-1)   = Tin%Q(Tin%flavor(    1        ,antiq))             
  T%Q(Ntot)     = Tin%Q(Tin%flavor(    1        ,quark))             

  if     (Ntot.eq.3) then ! 3-POINT AMPLITUDES
!
  if     (Noff.eq.0) then
    rslt = amp_1g0_qb0_q0( hel(per(1:2)) ,T )
  elseif (Noff.eq.1) then
    if     (offshell(1).eq.Tin%flavor(1,quark)) then ;rslt = amp_1g0_qb0_q1(hel(per(1:2)),T)
    elseif (offshell(1).eq.Tin%flavor(1,antiq)) then ;rslt = amp_1g0_qb1_q0(hel(per(1:2)),T)
                                                else ;rslt = amp_1g1_qb0_q0(hel(per(1:2)-1),T)
    endif
  else!if(Noff.eq.2) then
    if     (offshell(1).eq.Tin%flavor(1,quark)) then 
      if     (offshell(2).eq.Tin%flavor(1,antiq)) then ;rslt = amp_1g0_qb1_q1(hel(per(1:2)),T)
                                                  else ;rslt = amp_1g1_qb0_q1(hel(per(1:2)-1),T)
      endif
    elseif (offshell(1).eq.Tin%flavor(1,antiq)) then
      if     (offshell(2).eq.Tin%flavor(1,quark)) then ;rslt = amp_1g0_qb1_q1(hel(per(1:2)),T)
                                                  else ;rslt = amp_1g1_qb1_q0(hel(per(1:2)-1),T)
      endif
    else!if(offshell(1).eq.Tin%flavor(1,gluon)) then
      if     (offshell(2).eq.Tin%flavor(1,quark)) then ;rslt = amp_1g1_qb0_q1(hel(per(1:2)-1),T)
                                                  else ;rslt = amp_1g1_qb1_q0(hel(per(1:2)-1),T)
      endif
    endif
  endif
!
!
  elseif (Ntot.eq.4) then ! 4-POINT AMPLITUDES
!
  if     (Noff.eq.0) then
    rslt = amp_onshell_mhv( hel(per(1:3)) ,T )
  elseif (Noff.eq.1) then
    if     (offshell(1).eq.Tin%flavor(1,quark)) then ;rslt = amp_2g0_qb0_q1(hel(per(1:3)),T)
    elseif (offshell(1).eq.Tin%flavor(1,antiq)) then ;rslt = amp_2g0_qb1_q0(hel(per(1:3)),T)
                                                else ;rslt = amp_2g1_qb0_q0(hel(per(1:3)-1),T)
    endif
  else!if(Noff.eq.2) then
    if     (offshell(1).eq.Tin%flavor(1,quark)) then 
      if     (offshell(2).eq.Tin%flavor(1,antiq)) then ;rslt = amp_2g0_qb1_q1(hel(per(1:3)),T)
                                                  else ;rslt = amp_2g1_qb0_q1(hel(per(1:3)-1),T)
      endif
    elseif (offshell(1).eq.Tin%flavor(1,antiq)) then
      if     (offshell(2).eq.Tin%flavor(1,quark)) then ;rslt = amp_2g0_qb1_q1(hel(per(1:3)),T)
                                                  else ;rslt = amp_2g1_qb1_q0(hel(per(1:3)-1),T)
      endif
    else!if(offshell(1).eq.Tin%flavor(1,gluon)) then
      if     (offshell(2).eq.Tin%flavor(1,quark)) then ;rslt = amp_2g1_qb0_q1(hel(per(1:3)-1),T)
      elseif (offshell(2).eq.Tin%flavor(1,antiq)) then ;rslt = amp_2g1_qb1_q0(hel(per(1:3)-1),T)
                                                  else ;rslt = amp_2g2_qb0_q0(hel(per(1:3)-2),T)
      endif
    endif
  endif
!
!
  elseif (Ntot.eq.5) then ! 5-POINT AMPLITUDES
!
  if     (Noff.eq.0) then
    rslt = amp_onshell_mhv( hel(per(1:4)) ,T )
  elseif (Noff.eq.1) then
    if     (offshell(1).eq.Tin%flavor(1,quark)) then ;rslt=amp_3g0_qb0_q1(hel(per(1:4)),T)
    elseif (offshell(1).eq.Tin%flavor(1,antiq)) then ;rslt=amp_3g0_qb1_q0(hel(per(1:4)),T)
                                                else ;rslt=amp_3g1_qb0_q0(hel(per(1:4)-1),T)
    endif
  else
  endif
!
!
  else
  endif
  end associate
  end function 
  

  function amp_1g0_qb0_q0( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(2)
  complex(fltknd) :: rslt,xx,yy
  integer :: i1,iq,iqb
!
  i1=1 ;iqb=2 ;iq=3
!
  rslt = 0
  if (helicity(iqb).eq.-1) then
    if (helicity(i1).eq.1) then ! g+ qb- q+
      xx = T%sqr(i1,iq)  ;yy = T%sqr(iqb,iq)
    else ! g- qb- q+
      xx = T%ang(i1,iqb)  ;yy = T%ang(iq,iqb)
    endif
  else
    if (helicity(i1).eq.1) then ! g+ qb+ q-
      xx = T%sqr(i1,iqb)  ;yy = T%sqr(iqb,iq)
    else ! g- qb+ q-
      xx = T%ang(i1,iq)  ;yy = T%ang(iq,iqb)
    endif
  endif
  if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy
  end function


  function amp_1g0_qb0_q1( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(2)
  complex(fltknd) :: rslt,xx,yy
  integer :: i1,oq,iqb
!
  i1=1 ;iqb=2 ;oq=3
!
  rslt = 0
  if (helicity(i1).eq.1) then
    if (helicity(iqb).eq.-1) then ! g+ qb- q*
      xx = T%sqr(i1,oq)  ;yy = T%sqr(iqb,oq)
      if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy / T%Q(oq)%kapp
    else ! g- qb- q*
    endif
  else
    if (helicity(i1).eq.1) then ! g+ qb+ q*
    else ! g- qb+ q*
      xx = T%ang(i1,oq)  ;yy = T%ang(oq,iqb)
      if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy / T%Q(oq)%kstr
    endif
  endif
  end function


  function amp_1g0_qb1_q0( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(2)
  complex(fltknd) :: rslt,xx,yy
  integer :: i1,iq,oqb
!
  i1=1 ;oqb=2 ;iq=3
!
  rslt = 0
  if (helicity(oqb).eq.-1) then
    if (helicity(i1).eq.1) then ! g+ qb* q+
    else ! g- qb* q+
      xx = T%ang(i1,oqb)  ;yy = T%ang(iq,oqb)
      if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy / T%Q(oqb)%kstr
    endif
  else
    if (helicity(i1).eq.1) then ! g+ qb* q-
      xx = T%sqr(i1,oqb)  ;yy = T%sqr(oqb,iq)
      if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy / T%Q(oqb)%kapp
    else ! g- qb* q-
    endif
  endif
  end function


  function amp_1g1_qb0_q0( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(2)
  complex(fltknd) :: rslt,xx,yy
  integer :: o1,iq,iqb
!
  o1=1 ;iqb=2 ;iq=3
!
  rslt = 0
  if (helicity(iqb).eq.-1) then ! g* qb- q+
    xx = T%sqr(o1,iq)  ;yy = T%sqr(iqb,iq)
    if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy / T%Q(o1)%kapp
  else ! g* qb+ q-
    xx = T%ang(o1,iq)  ;yy = T%ang(iq,iqb)
    if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy / T%Q(o1)%kstr
  endif
  end function


 
  function amp_1g0_qb1_q1( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(2)
  complex(fltknd) :: rslt
  integer :: i1,oq,oqb
!
  i1=1 ;oqb=2 ;oq=3
!
  if (helicity(i1).eq.-1) then
    if (helicity(oqb).eq.-1) then ! g- qb*(-) q*(+)
      rslt = T%sqr(oq,oqb)**2 /( T%sqr(i1,oq)*T%sqr(oqb,i1)*T%Q(oq)%kapp )
    else ! g- qb*(+) q*(-)
      rslt = T%sqr(oq,oqb)**2 /( T%sqr(i1,oq)*T%sqr(oqb,i1)*T%Q(oqb)%kapp )
    endif
  else
    if (helicity(oqb).eq.-1) then ! g+ qb*(-) q*(+)
      rslt = T%ang(oq,oqb)**2 /( T%ang(i1,oq)*T%ang(oqb,i1)*T%Q(oqb)%kstr )
    else ! g+ qb*(+) q*(-)
      rslt = T%ang(oq,oqb)**2 /( T%ang(i1,oq)*T%ang(oqb,i1)*T%Q(oq)%kstr )
    endif
  endif
  end function

  function amp_1g1_qb1_q0( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(2)
  complex(fltknd) :: rslt,xx,yy
  integer :: o1,iq,oqb
!
  o1=1 ;oqb=2 ;iq=3
!
  rslt = 0
  if (helicity(oqb).eq.-1) then ! g* qb* q+
    xx = T%ang(o1,oqb)  ;yy = T%ang(iq,oqb)
    if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy/( T%Q(o1)%kstr*T%Q(oqb)%kstr )
  else ! g* qb* q-
    xx = T%sqr(o1,oqb)  ;yy = T%sqr(oqb,iq)
    if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy/( T%Q(o1)%kapp*T%Q(oqb)%kapp )
  endif
  end function

  function amp_1g1_qb0_q1( helicity ,T ) result(rslt) 
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(2)
  complex(fltknd) :: rslt,xx,yy
  integer :: o1,iqb,oq
!
  o1=1 ;iqb=2 ;oq=3
!
  rslt = 0
  if (helicity(iqb).eq.-1) then ! g* qb- q*
    xx = T%sqr(o1,oq)  ;yy = T%sqr(iqb,oq)
    if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy/( T%Q(o1)%kapp*T%Q(oq)%kapp )
  else ! g* qb+ q*
    xx = T%ang(o1,oq)  ;yy = T%ang(oq,iqb)
    if (xx.ne.0.and.yy.ne.0) rslt = xx*xx/yy/( T%Q(o1)%kstr*T%Q(oq)%kstr )
  endif
  end function



  function amp_onshell_mhv( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(T%Ntot-1)
  complex(fltknd) :: rslt
  integer :: ii,Nmin,Npls,plus(NsizeProc),mins(NsizeProc)
  associate( N=>T%Ntot )
  Nmin = 0
  Npls = 0
  do ii=1,N-1
    if (helicity(ii).eq. 1) then
      Npls = Npls+1
      plus(Npls) = ii
    else!if(helicity(ii).eq.-1) then
      Nmin = Nmin+1
      mins(Nmin) = ii
    endif
  enddo
  if (helicity(N-1).eq.-1) then
    Npls = Npls+1
    plus(Npls) = N
  else
    Nmin = Nmin+1
    mins(Nmin) = N
  endif
  if (Npls.le.1.or.Nmin.le.1) then
    rslt = 0
  elseif (Nmin.eq.2) then
    if     (mins(2).eq.N) then ! ... qb+ q-
      rslt = T%ang(N,mins(1))**3 * T%ang(N-1,mins(1))    / T%ang(N,1)
    else!if(mins(2).eq.N-1) then ! ... qb- q+
      rslt = T%ang(N,mins(1))    * T%ang(N-1,mins(1))**3 / T%ang(N,1)
    endif
    do ii=1,N-1
      rslt = rslt / T%ang(ii,ii+1)
    enddo 
  elseif (Npls.eq.2) then
    if     (plus(2).eq.N) then ! ... qb- q+
      rslt = T%sqr(N,plus(1))**3 * T%sqr(N-1,plus(1))    / T%sqr(1,N)
    else!if(plus(2).eq.N-1) then ! ... qb+ q-
      rslt = T%sqr(N,plus(1))    * T%sqr(N-1,plus(1))**3 / T%sqr(1,N)
    endif
    do ii=2,N
      rslt = rslt / T%sqr(ii,ii-1)
    enddo 
  else
  endif
  end associate
  end function


  function amp_2g0_qb0_q1( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(3)
  complex(fltknd) :: rslt
  integer :: i1,i2,oq,iqb
!
  i1=1 ;i2=2 ;iqb=3 ;oq=4
!
  if (helicity(iqb).eq.1) then
    if (helicity(i1).eq.1) then
      if (helicity(i2).eq.1) then ! g+ g+ qb+ q*
        rslt = 0
      else ! g+ g- qb+ q*
        rslt = T%ang(oq,i2)**3 / T%Q(oq)%kstr &
               / ( T%ang(oq,iqb) * T%ang(i1,oq) * T%ang(i2,i1) )
      endif
    else
      if (helicity(i2).eq.1) then ! g- g+ qb+ q*
        rslt =-T%ang(oq,i1)**2 /T%Q(oq)%kstr * T%ang(iqb,i1) &
               / ( T%ang(oq,iqb) * T%ang(i2,i1) * T%ang(iqb,i2) )
      else ! g- g- qb+ q*
        rslt = T%sqr(iqb,oq)**2 / ( T%sqr(i1,oq) * T%sqr(i1,i2) * T%sqr(iqb,i2) )
      endif
    endif
  else
    if (helicity(i1).eq.1) then
      if (helicity(i2).eq.1) then ! g+ g+ qb- q*
        rslt = T%ang(oq,iqb)**2 / ( T%ang(oq,i1) * T%ang(i2,i1) * T%ang(i2,iqb) )
      else ! g+ g- qb- q*
        rslt =-T%sqr(i1,oq)**2 /T%Q(oq)%kapp * T%sqr(i1,iqb) &
               / ( T%sqr(iqb,oq) * T%sqr(i1,i2) * T%sqr(i2,iqb) )
      endif
    else
      if (helicity(i2).eq.1) then ! g- g+ qb- q*
        rslt = T%sqr(i2,oq)**3 / T%Q(oq)%kapp &
               / ( T%sqr(iqb,oq) * T%sqr(oq,i1) * T%sqr(i1,i2) )
      else ! g- g- qb- q*
        rslt = 0
      endif
    endif
  endif
  end function


  function amp_2g0_qb1_q0( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(3)
  complex(fltknd) :: rslt
  integer :: i1,i2,iq,oqb
!
  i1=1 ;i2=2 ;oqb=3 ;iq=4
!
  if (helicity(oqb).eq.-1) then
    if (helicity(i1).eq.1) then
      if (helicity(i2).eq.1) then ! g+ g+ qb* q+
        rslt = 0
      else ! g+ g- qb* q+
        rslt =-T%ang(oqb,i2)**2 /T%Q(oqb)%kstr * T%ang(iq,i2) &
               / ( T%ang(oqb,iq) * T%ang(i1,i2) * T%ang(iq,i1) )
      endif
    else
      if (helicity(i2).eq.1) then ! g- g+ qb* q+
        rslt = T%ang(oqb,i1)**3 / T%Q(oqb)%kstr &
               / ( T%ang(oqb,iq) * T%ang(i2,oqb) * T%ang(i1,i2) )
      else ! g- g- qb* q+
        rslt = T%sqr(iq,oqb)**2 / ( T%sqr(i2,oqb) * T%sqr(i2,i1) * T%sqr(iq,i1) )
      endif
    endif
  else
    if (helicity(i1).eq.1) then
      if (helicity(i2).eq.1) then ! g+ g+ qb* q-
        rslt = T%ang(oqb,iq)**2 / ( T%ang(oqb,i2) * T%ang(i1,i2) * T%ang(i1,iq) )
      else ! g+ g- qb* q-
        rslt = T%sqr(i1,oqb)**3 / T%Q(oqb)%kapp &
               / ( T%sqr(iq,oqb) * T%sqr(oqb,i2) * T%sqr(i2,i1) )
      endif
    else
      if (helicity(i2).eq.1) then ! g- g+ qb* q-
        rslt =-T%sqr(i2,oqb)**2 /T%Q(oqb)%kapp * T%sqr(i2,iq) &
               / ( T%sqr(iq,oqb) * T%sqr(i2,i1) * T%sqr(i1,iq) )
      else ! g- g- qb* q-
        rslt = 0
      endif
    endif
  endif
  end function


  function amp_2g1_qb0_q0( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(3)
  complex(fltknd) :: rslt
  integer :: o1,o2,i1,i2,iq,iqb
!
  if (helicity(1).eq.0) then
!
  o1=1 ;i2=2 ;iqb=3 ;iq=4
!
  if (helicity(iqb).eq.-1) then
    if (helicity(i2).eq.1) then ! g* g+ qb- q+
      rslt = T%ang(iqb,o1)**3 / T%Q(o1)%kstr &
             / ( T%ang(iqb,iq) * T%ang(i2,iqb) * T%ang(o1,i2) )
    else ! g* g- qb- q+
      rslt =-T%sqr(o1,iq)**2 /T%Q(o1)%kapp * T%sqr(o1,iqb) &
             / ( T%sqr(iqb,iq) * T%sqr(o1,i2) * T%sqr(i2,iqb) )
    endif
  else
    if (helicity(i2).eq.1) then ! g* g+ qb+ q-
      rslt =-T%ang(iq,o1)**2 /T%Q(o1)%kstr * T%ang(iqb,o1) &
             / ( T%ang(iq,iqb) * T%ang(i2,o1) * T%ang(iqb,i2) )
    else ! g* g- qb+ q-
      rslt = T%sqr(o1,iqb)**3 / T%Q(o1)%kapp &
             / ( T%sqr(iq,iqb) * T%sqr(iqb,i2) * T%sqr(i2,o1) )
    endif
  endif
!
  else!if (helicity(2).eq.0) then
!
  i1=1 ;o2=2 ;iqb=3 ;iq=4
!
  if (helicity(iqb).eq.-1) then
    if (helicity(i1).eq.1) then ! g+ g* qb- q+
      rslt =-T%ang(iqb,o2)**2 /T%Q(o2)%kstr * T%ang(iq,o2) &
             / ( T%ang(iqb,iq) * T%ang(i1,o2) * T%ang(iq,i1) )
    else ! g- g* qb- q+
      rslt = T%sqr(o2,iq)**3 / T%Q(o2)%kapp &
             / ( T%sqr(iqb,iq) * T%sqr(iq,i1) * T%sqr(i1,o2) )
    endif
  else
    if (helicity(i1).eq.1) then ! g+ g* qb+ q-
      rslt = T%ang(iq,o2)**3 / T%Q(o2)%kstr &
             / ( T%ang(iq,iqb) * T%ang(i1,iq) * T%ang(o2,i1) )
    else ! g- g* qb+ q-
      rslt =-T%sqr(o2,iqb)**2 /T%Q(o2)%kapp * T%sqr(o2,iq) &
             / ( T%sqr(iq,iqb) * T%sqr(o2,i1) * T%sqr(i1,iq) )
    endif
  endif
!
  endif
!
  end function



  function amp_2g0_qb1_q1( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(3)
  complex(fltknd) :: rslt
  integer :: i1,i2,oqb,oq
  i1=1 ;i2=2 ;oqb=3 ;oq=4
  if (helicity(i1).eq.-1) then
    if (helicity(i2).eq.-1) then
      if (helicity(oqb).eq.-1) then ! g- g- qb*(-) q*(+)
        !write(*,*) '---+' !DEBUG
        rslt = - T%sqr(oqb,oq)**2 &
               /( T%sqr(oqb,i2)*T%sqr(i2,i1)*T%sqr(i1,oq)*T%Q(oq)%kapp )
      else ! g- g- qb*(+) q*(-)
        !write(*,*) '--+-' !DEBUG
        rslt = - T%sqr(oqb,oq)**2 &
               /( T%sqr(oqb,i2)*T%sqr(i2,i1)*T%sqr(i1,oq)*T%Q(oqb)%kapp )
      endif
    else
      if (helicity(oqb).eq.-1) then ! g- g+ qb*(-) q*(+)
        !write(*,*) '-+-+' !DEBUG
        rslt = +   T%ang(i1,oqb)**3 * T%sqr(oq,oqb)**2 / T%Q(oq)%kapp &
                /( T%ang(i2,oqb)*T%ang(i1,i2)*T%ang(oqb,[i1,i2],oq) ) &
                /( T%ang(oqb,[i1,i2],oqb) ) &
               +   T%sqr(i2,oq)**3 * T%ang(oq,oqb)**2 / T%Q(oqb)%kstr &
                /( T%sqr(oq,i1)*T%sqr(i1,i2)*T%ang(oqb,[i1,i2],oq))&
                /( T%ang(oq,[i1,i2],oq) )

      else ! g- g+ qb*(+) q*(-)
        !write(*,*) '-++-' !DEBUG
        rslt =    T%ang(i1,oqb)**4 * T%sqr(oq,oqb)**2 &
                /( T%ang(i1,i2)*T%ang(i2,oqb)*T%ang(i1,[oqb,i2],oqb)) &
                /( T%ang(oqb,[i1,i2],oq) * T%ang(oqb,[i1,i2],oqb))&
              +    T%sqr(i2,oq)**4 * T%ang(oq,oqb)**2 &
                /( T%sqr(i1,i2)*T%sqr(i1,oq)*T%ang(oq,[oq,i1],i2)) &
                /( T%ang(oqb,[i1,i2],oq) * T%ang(oq,[i1,i2],oq) ) &
              +   T%sqr(i2,oqb)**2 * T%ang(i1,oq)**2 * T%ang(i1,oqb,i2) &
               /( T%Q(oqb)%kapp*T%Q(oq)%kstr*T%sinv(i2,oqb) ) &
               /( T%ang(i1,oqb,i2)*T%ang(oq,[i2,oqb],oqb) &
                 + T%sinv(i2,oqb)*T%ang(i1,oq)*T%sqr(i2,oqb) )
      endif
    endif
  else
    if (helicity(i2).eq.-1) then
      if (helicity(oqb).eq.-1) then ! g+ g- qb*(-) q*(+)
        !write(*,*) '+--+' !DEBUG
        rslt =     T%sqr(oqb,i1)**4 * T%ang(oq,oqb)**2 &
                /( T%sqr(i1,i2)*T%sqr(i2,oqb)*T%ang(oqb,[oqb,i2],i1)*T%ang(oq,[i1,i2],oqb) ) &
                /( T%ang(oqb,[i1,i2],oqb)) &
              +    T%ang(i2,oq)**4 * T%sqr(oq,oqb)**2 &
                /( T%ang(i1,i2)*T%ang(i1,oq)*T%ang(i2,[oq,i1],oq)*T%ang(oq,[i1,i2],oqb) ) &
                /( T%ang(oq,[i1,i2],oq)) &
             +    T%ang(i2,oqb)**2 * T%sqr(i1,oq)**2 * T%ang(i2,oqb,i1) &
               /( T%Q(oq)%kapp*T%Q(oqb)%kstr*T%sinv(i2,oqb) ) &
               /( T%ang(i2,oqb,i1)*T%ang(oqb,[i2,oqb],oq) &
                 +T%sinv(i2,oqb)*T%sqr(i1,oq)*T%ang(i2,oqb) )
      else ! g+ g- qb*(+) q*(-)
        !write(*,*) '+-+-' !DEBUG
        rslt =     T%sqr(oqb,i1)**3 * T%ang(oq,oqb)**2 / T%Q(oq)%kstr &
                /( T%sqr(i1,i2)*T%sqr(i2,oqb)*T%ang(oq,[i1,i2],oqb) ) &
                /( T%ang(oqb,[i1,i2],oqb) )&
                +  T%ang(oq,i2)**3 * T%sqr(oq,oqb)**2 / T%Q(oqb)%kapp &
                /( T%ang(i1,i2)*T%ang(oq,i1)*T%ang(oq,[i1,i2],oqb) )&
                /( T%ang(oq,[i1,i2],oq) )
     endif
    else
      if (helicity(oqb).eq.-1) then ! g+ g+ qb*(-) q*(+)
        !write(*,*) '++-+' !DEBUG
        rslt = - T%ang(oq,oqb)**2 &
               /( T%ang(i1,i2)*T%ang(i2,oqb)*T%ang(oq,i1)*T%Q(oqb)%kstr )
      else ! g+ g+ qb*(+) q*(-)
        !write(*,*) '+++-' !DEBUG
        rslt = - T%ang(oq,oqb)**2 &
               /( T%ang(i1,i2)*T%ang(i2,oqb)*T%ang(oq,i1)*T%Q(oq)%kstr )
      endif
    endif
  endif
  end function


  function amp_2g1_qb0_q1( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(3)
  complex(fltknd) :: rslt
  integer :: o1,o2,i1,i2,iqb,oq
  if (helicity(1).eq.0) then
!
  o1=1 ;i2=2 ;iqb=3 ;oq=4
!
  if (helicity(i2).eq.-1) then
    if (helicity(iqb).eq.-1) then ! g* g- qb- q*
      rslt = T%sqr(o1,oq)**3 * T%sqr(o1,iqb) /( T%Q(o1)%kapp*T%Q(oq)%kapp ) &
             /( T%sqr(o1,oq)*T%sqr(oq,iqb)*T%sqr(iqb,i2)*T%sqr(i2,o1) )
    else ! g* g- qb+ q*
      rslt = T%sqr(oq,iqb)**2 * T%ang(o1,oq)**2 / T%Q(o1)%kstr &
             /( T%sqr(i2,iqb)*T%ang(oq,o1,i2)*T%ang(oq,o1,oq) ) &
           + T%ang(oq,[o1,i2],o1)**3 /( T%Q(o1)%kapp*T%Q(oq)%kstr*T%sinv(o1,i2) ) &
             /( T%sqr(o1,i2)*T%ang(oq,iqb)*T%ang(oq,o1,i2) )
    endif
  else
    if (helicity(iqb).eq.-1) then ! g* g+ qb- q*
      rslt = T%ang(oq,iqb)**2 * T%sqr(o1,oq)**2 / T%Q(o1)%kapp &
             /( T%ang(iqb,i2)*T%ang(i2,o1,oq)*T%ang(oq,o1,oq) ) &
           + T%ang(o1,[o1,i2],oq)**3 /( T%Q(o1)%kstr*T%Q(oq)%kapp*T%sinv(o1,i2) ) &
             /( T%ang(o1,i2)*T%sqr(oq,iqb)*T%ang(i2,o1,oq) )
    else ! g* g+ qb+ q*
      rslt = T%ang(o1,oq)**3 * T%ang(o1,iqb) /( T%Q(o1)%kstr*T%Q(oq)%kstr ) &
             /( T%ang(o1,i2)*T%ang(i2,iqb)*T%ang(iqb,oq)*T%ang(oq,o1) )
    endif
  endif
!
  else
!
  i1=1 ;o2=2 ;iqb=3 ;oq=4
!
  if (helicity(i1).eq.-1) then
    if (helicity(iqb).eq.-1) then ! g- g* qb- q*
      rslt = T%sqr(o2,oq)**3 * T%sqr(o2,iqb) /( T%Q(o2)%kapp*T%Q(oq)%kapp ) &
             /( T%sqr(i1,oq)*T%sqr(oq,iqb)*T%sqr(iqb,o2)*T%sqr(o2,i1) )
    else ! g- g* qb+ q*
      rslt = T%sqr(oq,iqb)**2 * T%ang(o2,oq)**2 / T%Q(o2)%kstr &
             /( T%ang(oq,o2,iqb)*T%sqr(oq,i1)*T%ang(oq,o2,i1) ) &
           - ( T%sqr(o2,iqb)**2 * T%ang(i1,oq)**2 /(T%sinv(o2,iqb)*T%ang(oq,o2,iqb)) &
              -T%ang(oq,[o2,i1],o2)**3 &
               /( T%sinv(i1,o2)*T%sqr(i1,o2)*T%ang(oq,iqb)*T%ang(oq,o2,i1) ) &
             )/( T%Q(o2)%kapp*T%Q(oq)%kstr )
    endif
  else
    if (helicity(iqb).eq.-1) then ! g+ g* qb- q*
      rslt =  T%ang(iqb,oq)**2 * T%sqr(oq,o2)**2 / T%Q(o2)%kapp &
             /( T%ang(iqb,o2,oq)*T%ang(i1,oq)*T%ang(i1,o2,oq) ) &
           - ( T%ang(iqb,o2)**2 * T%sqr(oq,i1)**2 /(T%sinv(o2,iqb)*T%ang(iqb,o2,oq)) &
              -T%ang(o2,[o2,i1],oq)**3 &
               /( T%sinv(i1,o2)*T%ang(i1,o2)*T%sqr(oq,iqb)*T%ang(i1,o2,oq) ) &
             )/( T%Q(o2)%kstr*T%Q(oq)%kapp )
    else ! g+ g* qb+ q*
      rslt = T%ang(o2,oq)**3 * T%ang(o2,iqb) /( T%Q(o2)%kstr*T%Q(oq)%kstr ) &
             /( T%ang(o2,iqb)*T%ang(iqb,oq)*T%ang(oq,i1)*T%ang(i1,o2) )
    endif
  endif
!
  endif
  end function


  function amp_2g1_qb1_q0( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(3)
  complex(fltknd) :: rslt
  integer :: o1,o2,i1,i2,oqb,iq
  if (helicity(1).eq.0) then
!
  o1=1 ;i2=2 ;oqb=3 ;iq=4
!
  if (helicity(i2).eq.-1) then
    if (helicity(oqb).eq.-1) then ! g* g- qb* q+
      rslt = T%sqr(oqb,iq)**2 * T%ang(o1,oqb)**2 / T%Q(o1)%kstr &
             /( T%ang(oqb,o1,iq)*T%sqr(oqb,i2)*T%ang(oqb,o1,i2) ) &
           - ( T%sqr(o1,iq)**2 * T%ang(i2,oqb)**2 /(T%sinv(o1,iq)*T%ang(oqb,o1,iq)) &
              -T%ang(oqb,[o1,i2],o1)**3 &
               /( T%sinv(i2,o1)*T%sqr(i2,o1)*T%ang(oqb,iq)*T%ang(oqb,o1,i2) ) &
             )/( T%Q(o1)%kapp*T%Q(oqb)%kstr )
    else ! g* g- qb* q-
      rslt = T%sqr(o1,oqb)**3 * T%sqr(o1,iq) /( T%Q(o1)%kapp*T%Q(oqb)%kapp ) &
             /( T%sqr(o1,i2)*T%sqr(i2,oqb)*T%sqr(oqb,iq)*T%sqr(iq,o1) )
    endif
  else
    if (helicity(oqb).eq.-1) then ! g* g+ qb* q+
      rslt = T%ang(o1,oqb)**3 * T%ang(o1,iq) /( T%Q(o1)%kstr*T%Q(oqb)%kstr ) &
             /( T%ang(o1,i2)*T%ang(i2,oqb)*T%ang(oqb,iq)*T%ang(iq,o1) )
    else ! g* g+ qb* q-
     rslt =  T%ang(iq,oqb)**2 * T%sqr(oqb,o1)**2 / T%Q(o1)%kapp &
             /( T%ang(iq,o1,oqb)*T%ang(i2,oqb)*T%ang(i2,o1,oqb) ) &
           - ( T%ang(iq,o1)**2 * T%sqr(oqb,i2)**2 /(T%sinv(o1,iq)*T%ang(iq,o1,oqb)) &
              -T%ang(o1,[o1,i2],oqb)**3 &
               /( T%sinv(i2,o1)*T%ang(i2,o1)*T%sqr(oqb,iq)*T%ang(i2,o1,oqb) ) &
             )/( T%Q(o1)%kstr*T%Q(oqb)%kapp )
    endif
  endif
!
  else
!
  i1=1 ;o2=2 ;oqb=3 ;iq=4
!
  if (helicity(i1).eq.-1) then
    if (helicity(oqb).eq.-1) then ! g- g* qb* q+
      rslt = T%sqr(oqb,iq)**2 * T%ang(o2,oqb)**2 / T%Q(o2)%kstr &
             /( T%sqr(i1,iq)*T%ang(oqb,o2,i1)*T%ang(oqb,o2,oqb) ) &
           + T%ang(oqb,[o2,i1],o2)**3 /( T%Q(o2)%kapp*T%Q(oqb)%kstr*T%sinv(o2,i1) ) &
             /( T%sqr(o2,i1)*T%ang(oqb,iq)*T%ang(oqb,o2,i1) )
    else ! g- g* qb* q-
      rslt = T%sqr(o2,oqb)**3 * T%sqr(o2,iq) /( T%Q(o2)%kapp*T%Q(oqb)%kapp ) &
             /( T%sqr(o2,oqb)*T%sqr(oqb,iq)*T%sqr(iq,i1)*T%sqr(i1,o2) )
    endif
  else
    if (helicity(oqb).eq.-1) then ! g+ g* qb* q+
      rslt = T%ang(o2,oqb)**3 * T%ang(o2,iq) /( T%Q(o2)%kstr*T%Q(oqb)%kstr ) &
             /( T%ang(o2,i1)*T%ang(i1,iq)*T%ang(iq,oqb)*T%ang(oqb,o2) )
    else ! g+ g* qb* q-
      rslt = T%ang(oqb,iq)**2 * T%sqr(o2,oqb)**2 / T%Q(o2)%kapp &
             /( T%ang(iq,i1)*T%ang(i1,o2,oqb)*T%ang(oqb,o2,oqb) ) &
           + T%ang(o2,[o2,i1],oqb)**3 /( T%Q(o2)%kstr*T%Q(oqb)%kapp*T%sinv(o2,i1) ) &
             /( T%ang(o2,i1)*T%sqr(oqb,iq)*T%ang(i1,o2,oqb) )
    endif
  endif
!
  endif
  end function


  function amp_2g2_qb0_q0( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(3)
  complex(fltknd) :: rslt
  integer :: o1,o2,iqb,iq
  o1=1 ;o2=2 ;iqb=3 ;iq=4
  if (helicity(iqb).eq.-1) then ! g* g* qb- q+
    !write(*,*) 'yo -+' !DEBUG
    rslt = T%sqr(o2,o1)**3 * T%ang(o1,iqb)**3 / T%Q(o2)%kapp &
           /( T%ang(iqb,iq)*T%ang(o1,[o1,o2],o2)*T%ang(iqb,o2,o1)*T%ang(o1,o2,o1) ) &
         + T%ang(o1,o2)**3 * T%sqr(o2,iq)**3 / T%Q(o1)%kstr &
           /( T%sqr(iq,iqb)*T%ang(o1,[o1,o2],o2)*T%ang(o2,o1,o2)*T%ang(o2,o1,iq) ) &
         - T%ang(o2,iqb)**2 * T%sqr(o1,iq)**2 * T%ang(o2,[o2,iqb],o1) &
           /( T%Q(o1)%kapp*T%Q(o2)%kstr*T%sinv(o2,iqb) ) & 
           /( T%ang(iqb,o2,iq)*T%ang(o2,[o2,iqb],o1) &
             +T%sinv(o2,iqb)*T%ang(iqb,o2)*T%sqr(iq,o1) )
  else ! g* g* qb+ q-
    !write(*,*) 'yo +-' !DEBUG
    rslt = T%ang(o1,o2)**3 * T%sqr(iqb,o1)**3 / T%Q(o2)%kstr &
           /( T%sqr(iq,iqb)*T%ang(o2,[o1,o2],o1)*T%ang(o1,o2,iqb)*T%ang(o1,o2,o1) ) &
         + T%sqr(o2,o1)**3 * T%ang(iq,o2)**3 / T%Q(o1)%kapp &
           /( T%ang(iqb,iq)*T%ang(o2,[o1,o2],o1)*T%ang(o2,o1,o2)*T%ang(iq,o1,o2) ) &
         - T%sqr(iqb,o2)**2 * T%ang(iq,o1)**2 * T%ang(o1,[o2,iqb],o2) &
           /( T%Q(o1)%kstr*T%Q(o2)%kapp*T%sinv(o2,iqb) ) & 
           /( T%ang(iq,o2,iqb)*T%ang(o1,[o2,iqb],o2) &
             +T%sinv(o2,iqb)*T%sqr(o2,iqb)*T%ang(o1,iq) )
!    rslt = T%ang(o2,o1)**3 * T%sqr(o2,iqb)**2 * T%sqr(o2,iq) / T%Q(o1)%kstr &
!           /( T%sqr(iq,iqb)*T%ang(o1,[o1,o2],o2)*T%ang(o2,o1,iq)*T%ang(o2,o1,o2) ) &
!         + T%sqr(o1,o2)**3 * T%ang(o1,iq)**2 * T%ang(o1,iqb) / T%Q(o2)%kapp &
!           /( T%ang(iqb,iq)*T%ang(o1,[o1,o2],o2)*T%ang(iqb,o2,o1)*T%ang(o1,o2,o1) ) &
!         + T%ang(o2,[o2,iqb],o1)**3 /( T%Q(o1)%kapp*T%Q(o2)%kstr*T%sinv(o2,iqb) ) & 
!           /( T%ang(iq,o2,iqb)*T%ang(o2,[o2,iqb],o1) &
!             +T%sinv(o2,iqb)*T%ang(iqb,o2)*T%sqr(iq,o1) )
  endif
  end function



  function amp_3g1_qb0_q0( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(4)
  complex(fltknd) :: rslt
  integer :: o1,o2,o3,i1,i2,i3,iq,iqb
!
  if (helicity(3).eq.0) then
!
  i1=1 ;i2=2 ;o3=3 ;iqb=4 ;iq=5
!
  if (helicity(iqb).eq.1) then
    if (helicity(i1).eq.1) then
      if (helicity(i2).eq.1) then ! g+ g+ g* qb+ q-
        rslt =-T%ang(o3,iq)**3 &
               /(T%ang(o3,i2)*T%ang(i1,i2)*T%ang(i1,iq)*T%ang(iq,iqb)*T%Q(o3)%kstr)
      else ! g+ g- g* qb+ q-
        rslt =-T%sqr(o3,i1)**4 * T%ang(o3,iq)**3 &
               /(T%sqr(o3,i2)*T%sqr(i1,i2)*T%ang(iq,iqb)*T%ang(o3,[o3,i2],i1)&
                 *T%ang(o3,[i1,i2],o3)*T%ang(iq,[i1,i2],o3)) & 
             + T%sqr(o3,iqb)**2  * T%ang(i2,iq)**3 * T%ang(i2,[o3,iqb],o3) &
               /( T%ang(i1,i2) * T%ang(i1,iq) * T%sinv(o3,iqb) * T%Q(o3)%kapp &
                 * (-T%ang(i2,[o3,iqb],o3) * T%ang(iq,o3,iqb) &
                     -( T%sqr(o3,iqb) * T%ang(i2,iq) * T%sinv(o3,iqb) ) ) ) &
             + ((T%sqr(i1,iqb)**3)*(T%ang(o3,i2)**4)) &
               /(T%sqr(iq,iqb)*T%ang(o3,[o3,i2],i1)*T%ang(i2,o3,iqb) &
                 *T%Q(o3)%kstr * T%ang(o3,[o3,i2],o3,i2) )
      endif
    else
      if (helicity(i2).eq.1) then ! g- g+ g* qb+ q-
        rslt = T%sqr(o3,iq) * T%sqr(o3,iqb)**2 * T%ang(o3,i1)**4 &
               /( T%sqr(iq,iqb) * T%ang(o3,i2) * T%ang(i1,i2) * T%ang(o3,[i1,i2],o3) *&
                  T%ang(o3,[i1,i2],iq)*T%ang(i1,[o3,i2],o3)) &
             + ((T%sqr(o3,i2)**4)*(T%ang(i1,iq)**2)*T%ang(i1,iqb))&
               /(T%ang(iq,iqb)*T%ang(i1,[o3,i2],o3)*T%ang(iqb,o3,i2)&
                 *T%sqr(o3,[o3,i2],o3,i2)*T%Q(o3)%kapp)&
             + T%sqr(i2,iq)*(T%ang(o3,[o3,iqb],i2)**3)&
               /(T%sqr(i1,i2)*T%sqr(i1,iq)*T%sinv(o3,iqb) &
                 *(-T%ang(o3,[o3,iqb],i2)*T%ang(iqb,o3,iq)-T%sqr(i2,iq)&
                 *T%ang(o3,iqb)*T%sinv(o3,iqb) )*T%Q(o3)%kstr)
      else ! g- g- g* qb+ q-
        rslt = -(T%sqr(o3,iq)*T%sqr(o3,iqb)**2)&
               /(T%sqr(o3,i2)*T%sqr(i1,i2)*T%sqr(i1,iq)*T%sqr(iq,iqb)*T%Q(o3)%kapp)
      endif
    endif
  else
    if (helicity(i1).eq.1) then
      if (helicity(i2).eq.1) then ! g+ g+ g* qb- q+
        rslt = -(T%ang(iq,o3)*T%ang(iqb,o3)**2)&
               /(T%ang(i2,o3)*T%ang(i2,i1)*T%ang(iq,i1)*T%ang(iqb,iq)*T%Q(o3)%kstr)
      else ! g+ g- g* qb- q+
        rslt = T%ang(iq,o3) * T%ang(iqb,o3)**2 * T%sqr(i1,o3)**4 &
               /( T%ang(iqb,iq) * T%sqr(i2,o3) * T%sqr(i2,i1) * T%ang(o3,[i2,i1],o3) *&
                  T%ang(iq,[i2,i1],o3)*T%ang(o3,[i2,o3],i1)) &
             + ((T%ang(i2,o3)**4)*(T%sqr(iq,i1)**2)*T%sqr(iqb,i1))&
               /(T%sqr(iqb,iq)*T%ang(o3,[i2,o3],i1)*T%ang(i2,o3,iqb)&
                 *T%ang(i2,o3,[i2,o3],o3)*T%Q(o3)%kstr)&
             + T%ang(iq,i2)*(T%ang(i2,[iqb,o3],o3)**3)&
               /(T%ang(i2,i1)*T%ang(iq,i1)*T%sinv(o3,iqb) &
                 *(-T%ang(i2,[iqb,o3],o3)*T%ang(iq,o3,iqb)-T%ang(iq,i2)&
                 *T%sqr(iqb,o3)*T%sinv(o3,iqb) )*T%Q(o3)%kapp)
      endif
    else
      if (helicity(i2).eq.1) then ! g- g+ g* qb- q+
        rslt =-T%ang(i1,o3)**4 * T%sqr(iq,o3)**3 &
               /(T%ang(i2,o3)*T%ang(i2,i1)*T%sqr(iqb,iq)*T%ang(i1,[o3,i2],o3)&
                 *T%ang(o3,[i1,i2],o3)*T%ang(o3,[i1,i2],iq)) & 
             + T%ang(iqb,o3)**2  * T%sqr(iq,i2)**3 * T%ang(o3,[o3,iqb],i2) &
               /( T%sqr(i2,i1) * T%sqr(iq,i1) * T%sinv(o3,iqb) * T%Q(o3)%kstr&
                 * (-T%ang(o3,[o3,iqb],i2) * T%ang(iqb,o3,iq) &
                     -( T%ang(iqb,o3) * T%sqr(iq,i2) * T%sinv(o3,iqb) ) ) ) &
             + ((T%ang(iqb,i1)**3)*(T%sqr(i2,o3)**4)) &
               /(T%ang(iqb,iq)*T%ang(i1,[o3,i2],o3)*T%ang(iqb,o3,i2) &
                 *T%Q(o3)%kapp * T%sqr(i2,o3,[o3,i2],o3))
      else ! g- g- g* qb- q+
        rslt =-T%sqr(iq,o3)**3 &
               /(T%sqr(i2,o3)*T%sqr(i2,i1)*T%sqr(iq,i1)*T%sqr(iqb,iq)*T%Q(o3)%kapp)
      endif
    endif
  endif
!
  elseif (helicity(1).eq.0) then
!
  o1=1 ;i2=2 ;i3=3 ;iqb=4 ;iq=5
!
  if (helicity(iqb).eq.-1) then
    if (helicity(i3).eq.-1) then
      if (helicity(i2).eq.-1) then ! g* g- g- qb- q+
        rslt = -(T%sqr(o1,iqb)*T%sqr(o1,iq)**2)&
               /(T%sqr(o1,i2)*T%sqr(i3,i2)*T%sqr(i3,iqb)*T%sqr(iqb,iq)*T%Q(o1)%kapp)
      else ! g* g- g+ qb- q+
        rslt = T%sqr(o1,iqb) * T%sqr(o1,iq)**2 * T%ang(o1,i3)**4 &
               /( T%sqr(iqb,iq) * T%ang(o1,i2) * T%ang(i3,i2) * T%ang(o1,[i3,i2],o1) *&
                  T%ang(o1,[i3,i2],iqb)*T%ang(i3,[o1,i2],o1)) &
             + ((T%sqr(o1,i2)**4)*(T%ang(i3,iqb)**2)*T%ang(i3,iq))&
               /(T%ang(iqb,iq)*T%ang(i3,[o1,i2],o1)*T%ang(iq,o1,i2)&
                 *T%sqr(o1,[o1,i2],o1,i2)*T%Q(o1)%kapp)&
             + T%sqr(i2,iqb)*(T%ang(o1,[o1,iq],i2)**3)&
               /(T%sqr(i3,i2)*T%sqr(i3,iqb)*T%sinv(o1,iq) &
                 *(-T%ang(o1,[o1,iq],i2)*T%ang(iq,o1,iqb)-T%sqr(i2,iqb)&
                 *T%ang(o1,iq)*T%sinv(o1,iq) )*T%Q(o1)%kstr)
      endif
    else
      if (helicity(i2).eq.-1) then ! g* g+ g- qb- q+
        rslt =-T%sqr(o1,i3)**4 * T%ang(o1,iqb)**3 &
               /(T%sqr(o1,i2)*T%sqr(i3,i2)*T%ang(iqb,iq)*T%ang(o1,[o1,i2],i3)&
                 *T%ang(o1,[i3,i2],o1)*T%ang(iqb,[i3,i2],o1)) & 
             + T%sqr(o1,iq)**2  * T%ang(i2,iqb)**3 * T%ang(i2,[o1,iq],o1) &
               /( T%ang(i3,i2) * T%ang(i3,iqb) * T%sinv(o1,iq) * T%Q(o1)%kapp &
                 * (-T%ang(i2,[o1,iq],o1) * T%ang(iqb,o1,iq) &
                     -( T%sqr(o1,iq) * T%ang(i2,iqb) * T%sinv(o1,iq) ) ) ) &
             + ((T%sqr(i3,iq)**3)*(T%ang(o1,i2)**4)) &
               /(T%sqr(iqb,iq)*T%ang(o1,[o1,i2],i3)*T%ang(i2,o1,iq) &
                 *T%Q(o1)%kstr * T%ang(o1,[o1,i2],o1,i2) )
      else ! g* g+ g+ qb- q+
        rslt =-T%ang(o1,iqb)**3 &
               /(T%ang(o1,i2)*T%ang(i3,i2)*T%ang(i3,iqb)*T%ang(iqb,iq)*T%Q(o1)%kstr)
      endif
    endif
  else
    if (helicity(i3).eq.-1) then
      if (helicity(i2).eq.-1) then ! g* g- g- qb+ q-
        rslt =-T%sqr(iqb,o1)**3 &
               /(T%sqr(i2,o1)*T%sqr(i2,i3)*T%sqr(iqb,i3)*T%sqr(iq,iqb)*T%Q(o1)%kapp)
      else ! g* g- g+ qb+ q-
        rslt =-T%ang(i3,o1)**4 * T%sqr(iqb,o1)**3 &
               /(T%ang(i2,o1)*T%ang(i2,i3)*T%sqr(iq,iqb)*T%ang(i3,[o1,i2],o1)&
                 *T%ang(o1,[i3,i2],o1)*T%ang(o1,[i3,i2],iqb)) & 
             + T%ang(iq,o1)**2  * T%sqr(iqb,i2)**3 * T%ang(o1,[o1,iq],i2) &
               /( T%sqr(i2,i3) * T%sqr(iqb,i3) * T%sinv(o1,iq) * T%Q(o1)%kstr&
                 * (-T%ang(o1,[o1,iq],i2) * T%ang(iq,o1,iqb) &
                     -( T%ang(iq,o1) * T%sqr(iqb,i2) * T%sinv(o1,iq) ) ) ) &
             + ((T%ang(iq,i3)**3)*(T%sqr(i2,o1)**4)) &
               /(T%ang(iq,iqb)*T%ang(i3,[o1,i2],o1)*T%ang(iq,o1,i2) &
                 *T%Q(o1)%kapp * T%sqr(i2,o1,[o1,i2],o1))
      endif
    else
      if (helicity(i2).eq.-1) then ! g* g+ g- qb+ q-
        rslt = T%ang(iqb,o1) * T%ang(iq,o1)**2 * T%sqr(i3,o1)**4 &
               /( T%ang(iq,iqb) * T%sqr(i2,o1) * T%sqr(i2,i3) * T%ang(o1,[i2,i3],o1) *&
                  T%ang(iqb,[i2,i3],o1)*T%ang(o1,[i2,o1],i3)) &
             + ((T%ang(i2,o1)**4)*(T%sqr(iqb,i3)**2)*T%sqr(iq,i3))&
               /(T%sqr(iq,iqb)*T%ang(o1,[i2,o1],i3)*T%ang(i2,o1,iq)&
                 *T%ang(i2,o1,[i2,o1],o1)*T%Q(o1)%kstr)&
             + T%ang(iqb,i2)*(T%ang(i2,[iq,o1],o1)**3)&
               /(T%ang(i2,i3)*T%ang(iqb,i3)*T%sinv(o1,iq) &
                 *(-T%ang(i2,[iq,o1],o1)*T%ang(iqb,o1,iq)-T%ang(iqb,i2)&
                 *T%sqr(iq,o1)*T%sinv(o1,iq) )*T%Q(o1)%kapp)
      else ! g* g+ g+ qb+ q-
        rslt = -(T%ang(iqb,o1)*T%ang(iq,o1)**2)&
               /(T%ang(i2,o1)*T%ang(i2,i3)*T%ang(iqb,i3)*T%ang(iq,iqb)*T%Q(o1)%kstr)
      endif
    endif
  endif
!
  else!if (helicity(2).eq.0) then
!
  i1=1 ;o2=2 ;i3=3 ;iqb=4 ;iq=5
!
  if (helicity(iqb).eq.1) then
    if (helicity(i3).eq.1) then
      if (helicity(i1).eq.1) then ! g+ g* g+ qb+ q-
        rslt = -((T%ang(o2,iq)**3)*T%ang(o2,iqb))&
               /(T%ang(o2,i3)*T%ang(o2,i1)*T%ang(i3,iqb)&
               *T%ang(i1,iq)*T%ang(iq,iqb)*T%Q(o2)%kstr)
      else ! g- g* g+ qb+ q-
        rslt = ( T%sqr(o2,iq) * T%ang(o2,iqb) *(T%ang(o2,[i1,iq],o2)**2))&
            /(T%sqr(o2,i1)*T%sqr(i1,iq)*T%ang(o2,i3)*T%ang(i3,iqb)&
            *T%ang(o2,[o2,i1],iq)*T%ang(iqb,[i1,iq],o2))&
            -((T%sqr(o2,i3)**3)*(T%ang(i1,iq)**2)*T%ang(i1,iqb))&
            /(T%ang(iq,iqb)*T%ang(i1,o2,i3)*T%ang(iqb,[o2,i3],o2)&
            *T%sinv(o2,i3)*T%Q(o2)%kapp )&
            +(T%sqr(i3,iq)*(T%sqr(i3,iqb)**2)*(T%ang(o2,i1)**4))&
            /(T%sqr(iq,iqb)*T%ang(o2,[o2,i1],iq)*T%ang(i1,o2,i3)&
            *T%ang(i1,o2,[o2,i1],o2)*T%Q(o2)%kstr)
      endif
    else
      if (helicity(i1).eq.1) then ! g+ g* g- qb+ q-
        rslt =-((T%sqr(o2,iqb)**3)*(T%ang(o2,iq)**3))&
            /(T%sqr(o2,i3)*T%sqr(i3,iqb)*T%ang(o2,i1)*T%ang(i1,iq)&
            *T%ang(o2,[i1,iq],iqb)*T%ang(iq,[o2,i1],o2))&
            -((T%sqr(o2,i1)**4)*(T%ang(i3,iq)**3))&
            /(T%ang(iq,iqb)*T%ang(i3,o2,i1)*T%ang(iq,[o2,i1],o2)&
            *T%sqr(o2,[o2,i1],o2,i1)*T%Q(o2)%kapp)&
            -((T%sqr(i1,iqb)**3)*(T%ang(o2,i3)**3))&
            /(T%sqr(iq,iqb)*T%ang(o2,[o2,i3],iqb)&
            *T%ang(i3,o2,i1)*T%sinv(o2,i3)*T%Q(o2)%kstr)
      else ! g- g* g- qb+ q-
        rslt = -(T%sqr(o2,iq)*(T%sqr(o2,iqb)**3))&
            /(T%sqr(o2,i3)*T%sqr(o2,i1)*T%sqr(i3,iqb)&
            *T%sqr(i1,iq)*T%sqr(iq,iqb)*T%Q(o2)%kapp)
      endif
    endif
  else
    if (helicity(i3).eq.1) then
      if (helicity(i1).eq.1) then ! g+ g* g+ qb- q+
        rslt = -(T%ang(iq,o2)*(T%ang(iqb,o2)**3))&
            /(T%ang(i3,o2)*T%ang(i1,o2)*T%ang(iqb,i3)&
            *T%ang(iq,i1)*T%ang(iqb,iq)*T%Q(o2)%kstr)
      else ! g- g* g+ qb- q+
        rslt =-((T%ang(iqb,o2)**3)*(T%sqr(iq,o2)**3))&
            /(T%ang(i3,o2)*T%ang(iqb,i3)*T%sqr(i1,o2)*T%sqr(iq,i1)&
            *T%ang(iqb,[iq,i1],o2)*T%ang(o2,[i1,o2],iq))&
            -((T%ang(i1,o2)**4)*(T%sqr(iq,i3)**3))&
            /(T%sqr(iqb,iq)*T%ang(i1,o2,i3)*T%ang(o2,[i1,o2],iq)&
            *T%ang(i1,o2,[i1,o2],o2)*T%Q(o2)%kstr)&
            -((T%ang(iqb,i1)**3)*(T%sqr(i3,o2)**3))&
            /(T%ang(iqb,iq)*T%ang(iqb,[i3,o2],o2)&
            *T%ang(i1,o2,i3)*T%sinv(o2,i3)*T%Q(o2)%kapp)
      endif
    else
      if (helicity(i1).eq.1) then ! g+ g* g- qb- q+
        rslt = ( T%ang(iq,o2) * T%sqr(iqb,o2) *(T%ang(o2,[iq,i1],o2)**2))&
            /(T%ang(i1,o2)*T%ang(iq,i1)*T%sqr(i3,o2)*T%sqr(iqb,i3)&
            *T%ang(iq,[i1,o2],o2)*T%ang(o2,[iq,i1],iqb))&
            -((T%ang(i3,o2)**3)*(T%sqr(iq,i1)**2)*T%sqr(iqb,i1))&
            /(T%sqr(iqb,iq)*T%ang(i3,o2,i1)*T%ang(o2,[i3,o2],iqb)&
            *T%sinv(o2,i3)*T%Q(o2)%kstr )&
            +(T%ang(iq,i3)*(T%ang(iqb,i3)**2)*(T%sqr(i1,o2)**4))&
            /(T%ang(iqb,iq)*T%ang(iq,[i1,o2],o2)*T%ang(i3,o2,i1)&
            *T%sqr(o2,[i1,o2],o2,i1)*T%Q(o2)%kapp)
      else ! g- g* g- qb- q+
        rslt = -((T%sqr(iq,o2)**3)*T%sqr(iqb,o2))&
               /(T%sqr(i3,o2)*T%sqr(i1,o2)*T%sqr(iqb,i3)&
               *T%sqr(iq,i1)*T%sqr(iqb,iq)*T%Q(o2)%kapp)
      endif
    endif
  endif
!
  endif
!
  end function


  function amp_3g0_qb0_q1( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(4)
  complex(fltknd) :: rslt
  integer :: i1,i2,i3,oq,iqb
!
  i1=1 ;i2=2 ;i3=3 ;iqb=4 ;oq=5
!
  if (helicity(iqb).eq.1) then
    if (helicity(i3).eq.1) then
      if (helicity(i2).eq.1) then
        if (helicity(i1).eq.1) then ! g+ g+ g+ qb+ q*
          rslt = 0
        else ! g- g+ g+ qb+ q*
           rslt =-(T%ang(oq,i1)**2)*T%ang(i1,iqb)/(T%ang(oq,iqb)*T%ang(i1,i2)*T%ang(i3,i2) &
            *T%ang(i3,iqb)*T%Q(oq)%kstr)
        endif
      else
        if (helicity(i1).eq.1) then ! g+ g- g+ qb+ q*
           rslt =(T%ang(i2,oq)**3)*T%ang(iqb,i2)/(T%ang(i1,oq)*T%ang(iqb,oq)*T%ang(i2,i1)&
            *T%ang(i2,i3)*T%ang(iqb,i3)*T%Q(oq)%kstr)
        else ! g- g- g+ qb+ q*
          rslt = - T%ang(iqb,i2) * T%ang(i2,[iqb,i3],oq)**2 &
            /( T%ang(i3,i2) * T%ang(iqb,i3) * T%sqr(i1,oq) * T%ang(iqb,oq,i1) * T%sinv(oq,i1) ) &
            - T%ang(iqb,oq,i3) * T%ang(oq,[iqb,oq],i3)**2 &
            /( T%ang(iqb,oq) * T%sqr(i1,i2) * T%sqr(i3,i2) * T%ang(iqb,oq,i1) * T%sinv(oq,iqb) * T%Q(oq)%kstr  )
        endif
      endif
    else
      if (helicity(i2).eq.1) then
        if (helicity(i1).eq.1) then ! g+ g+ g- qb+ q*
          rslt =T%ang(i3,oq)**3/(T%ang(i1,oq)*T%ang(iqb,oq)*T%ang(i2,i1)*T%ang(i2,i3)*T%Q(oq)%kstr)
        else ! g- g+ g- qb+ q*
          rslt = -((T%ang(i3,oq)**3)*(T%sqr(i2,oq)**4))&
            /(T%ang(iqb,oq)*T%sqr(i1,oq)*T%sqr(i1,i2)&
            *T%ang(oq,[i1,i2],oq)*T%ang(i3,[i1,i2],oq)*T%ang(oq,[i1,oq],i2))&
            -((T%ang(i3,i1)**4)*(T%sqr(iqb,oq)**2))&
            /(T%ang(i1,i2)*T%ang(i3,i2)*T%ang(i3,[iqb,oq],oq)&
            *T%ang(i1,oq,iqb)*T%sinv(oq,iqb))&
            +   T%ang(i1,oq)**2 * T%sqr(iqb,i2)**3 * T%ang(i1,oq,i2) &
            /( T%sqr(i3,i2) * T%sqr(iqb,i3) * T%ang(oq,[i1,oq],i2) &
            * T%ang(i1,oq,iqb) * T%sinv(oq,i1) * T%Q(oq)%kstr )
        endif
      else
        if (helicity(i1).eq.1) then ! g+ g- g- qb+ q*
          rslt = -((T%ang(i3,oq)**3)*(T%sqr(i1,oq)**3))&
            /(T%ang(iqb,oq)*T%sqr(i1,i2)*T%ang(oq,[i1,i2],oq)&
            *T%ang(i3,[i1,i2],oq)*T%ang(oq,[i1,oq],i2))&
            -((T%ang(i3,i2)**3)*T%sqr(iqb,oq)**2)&
            /(T%ang(i1,i2)*T%ang(i3,[iqb,oq],oq)*T%ang(i1,oq,iqb)*T%sinv(oq,iqb)) &
            +(T%ang(oq,[i1,oq],iqb)**3)&
            /(T%ang(i1,oq)*T%sqr(i3,i2)*T%sqr(iqb,i3)&
            *T%ang(oq,[i1,oq],i2)*T%ang(i1,oq,iqb)*T%Q(oq)%kstr)
        else ! g- g- g- qb+ q*
          rslt =T%sqr(iqb,oq)**2/(T%sqr(oq,i1)*T%sqr(i2,i1)*T%sqr(i2,i3)*T%sqr(iqb,i3))
        endif
      endif
    endif
  else
    if (helicity(i3).eq.1) then
      if (helicity(i2).eq.1) then
        if (helicity(i1).eq.1) then ! g+ g+ g+ qb- q*
          rslt =T%ang(oq,iqb)**2/(T%ang(i1,oq)*T%ang(i1,i2)*T%ang(i3,i2)*T%ang(i3,iqb))
        else ! g- g+ g+ qb- q*
          rslt = -((T%sqr(oq,i3)**3)*(T%ang(oq,i1)**3))&
            /(T%sqr(oq,iqb)*T%ang(i2,i1)*T%ang(oq,[i2,i1],oq)&
            *T%ang(oq,[i2,i1],i3)*T%ang(i2,[oq,i1],oq))&
            -((T%sqr(i2,i3)**3)*T%ang(oq,iqb)**2)&
            /(T%sqr(i2,i1)*T%ang(oq,[oq,iqb],i3)*T%ang(iqb,oq,i1)*T%sinv(oq,iqb)) &
            +(T%ang(iqb,[oq,i1],oq)**3)&
            /(T%sqr(oq,i1)*T%ang(i2,i3)*T%ang(i3,iqb)&
            *T%ang(i2,[oq,i1],oq)*T%ang(iqb,oq,i1)*T%Q(oq)%kapp)
        endif
      else
        if (helicity(i1).eq.1) then ! g+ g- g+ qb- q*
          rslt = -((T%sqr(oq,i3)**3)*(T%ang(oq,i2)**4))&
            /(T%sqr(oq,iqb)*T%ang(oq,i1)*T%ang(i2,i1)&
            *T%ang(oq,[i2,i1],oq)*T%ang(oq,[i2,i1],i3)*T%ang(i2,[oq,i1],oq))&
            -((T%sqr(i1,i3)**4)*(T%ang(oq,iqb)**2))&
            /(T%sqr(i2,i1)*T%sqr(i2,i3)*T%ang(oq,[oq,iqb],i3)&
            *T%ang(iqb,oq,i1)*T%sinv(oq,iqb))&
            +   T%sqr(oq,i1)**2 * T%ang(i2,iqb)**3 * T%ang(i2,oq,i1) &
            /( T%ang(i2,i3) * T%ang(i3,iqb) * T%ang(i2,[oq,i1],oq) &
            * T%ang(iqb,oq,i1) * T%sinv(oq,i1) * T%Q(oq)%kapp )
        else ! g- g- g+ qb- q*
          rslt =T%sqr(oq,i3)**3/(T%sqr(oq,i1)*T%sqr(oq,iqb)*T%sqr(i1,i2)*T%sqr(i3,i2)*T%Q(oq)%kapp)
        endif
      endif
    else
      if (helicity(i2).eq.1) then
        if (helicity(i1).eq.1) then ! g+ g+ g- qb- q*
          rslt = - T%sqr(i2,iqb) * T%ang(oq,[i3,iqb],i2)**2 &
            /( T%sqr(i2,i3) * T%sqr(i3,iqb) * T%ang(oq,i1) * T%ang(i1,oq,iqb) * T%sinv(oq,i1) ) &
            - T%ang(i3,oq,iqb) * T%ang(i3,[oq,iqb],oq)**2 &
            /( T%sqr(oq,iqb) * T%ang(i2,i1) * T%ang(i2,i3) * T%ang(i1,oq,iqb) * T%sinv(oq,iqb) * T%Q(oq)%kapp  )
        else ! g- g+ g- qb- q*
          rslt =(T%sqr(oq,i2)**3)*T%sqr(i2,iqb)/(T%sqr(oq,i1)*T%sqr(oq,iqb)*T%sqr(i1,i2)&
            *T%sqr(i3,i2)*T%sqr(i3,iqb)*T%Q(oq)%kapp)
        endif
      else
        if (helicity(i1).eq.1) then ! g+ g- g- qb- q*
          rslt =((T%sqr(oq,i1)**2)*T%sqr(i1,iqb))/(T%sqr(oq,iqb)*T%sqr(i1,i2)*T%sqr(i3,i2) &
            *T%sqr(i3,iqb)*T%Q(oq)%kapp)
        else ! g- g- g- qb- q*
          rslt = 0
        endif
      endif
    endif
  endif
  end function


  function amp_3g0_qb1_q0( helicity ,T ) result(rslt)
  type(qomentum_list_type),intent(in) :: T
  integer,intent(in) :: helicity(4)
  complex(fltknd) :: rslt
  integer :: i1,i2,i3,oqb,iq
!
  i1=1 ;i2=2 ;i3=3 ;oqb=4 ;iq=5
!
  if (helicity(oqb).eq.-1) then
    if (helicity(i1).eq.1) then
      if (helicity(i2).eq.1) then
        if (helicity(i3).eq.1) then ! g+ g+ g+ qb* q+
          rslt = 0
        else ! g+ g+ g- qb* q+
          rslt =-(T%ang(oqb,i3)**2)*T%ang(i3,iq)/(T%ang(oqb,iq)*T%ang(i3,i2)*T%ang(i1,i2) &
              *T%ang(i1,iq)*T%Q(oqb)%kstr)
        endif
      else
        if (helicity(i3).eq.1) then ! g+ g- g+ qb* q+
          rslt =(T%ang(i2,oqb)**3)*T%ang(iq,i2)/(T%ang(i3,oqb)*T%ang(iq,oqb)*T%ang(i2,i3)&
              *T%ang(i2,i1)*T%ang(iq,i1)*T%Q(oqb)%kstr)
        else ! g+ g- g- qb* q+
          rslt = - T%ang(iq,i2) * T%ang(i2,[iq,i1],oqb)**2 &
              /( T%ang(i1,i2) * T%ang(iq,i1) * T%sqr(i3,oqb) * T%ang(iq,oqb,i3) * T%sinv(oqb,i3) ) &
              - T%ang(iq,oqb,i1) * T%ang(oqb,[iq,oqb],i1)**2 &
              /( T%ang(iq,oqb) * T%sqr(i3,i2) * T%sqr(i1,i2) * T%ang(iq,oqb,i3) * T%sinv(oqb,iq) * T%Q(oqb)%kstr  )
        endif
      endif
    else
      if (helicity(i2).eq.1) then
        if (helicity(i3).eq.1) then ! g- g+ g+ qb* q+
          rslt =T%ang(i1,oqb)**3/(T%ang(i3,oqb)*T%ang(iq,oqb)*T%ang(i2,i3)*T%ang(i2,i1)*T%Q(oqb)%kstr)
        else ! g- g+ g- qb* q+
          rslt = -((T%ang(i1,oqb)**3)*(T%sqr(i2,oqb)**4))&
              /(T%ang(iq,oqb)*T%sqr(i3,oqb)*T%sqr(i3,i2)&
              *T%ang(oqb,[i3,i2],oqb)*T%ang(i1,[i3,i2],oqb)*T%ang(oqb,[i3,oqb],i2))&
              -((T%ang(i1,i3)**4)*(T%sqr(iq,oqb)**2))&
              /(T%ang(i3,i2)*T%ang(i1,i2)*T%ang(i1,[iq,oqb],oqb)&
              *T%ang(i3,oqb,iq)*T%sinv(oqb,iq))&
              +   T%ang(i3,oqb)**2 * T%sqr(iq,i2)**3 * T%ang(i3,oqb,i2) &
              /( T%sqr(i1,i2) * T%sqr(iq,i1) * T%ang(oqb,[i3,oqb],i2) &
              * T%ang(i3,oqb,iq) * T%sinv(oqb,i3) * T%Q(oqb)%kstr )
        endif
      else
        if (helicity(i3).eq.1) then ! g- g- g+ qb* q+
          rslt = -((T%ang(i1,oqb)**3)*(T%sqr(i3,oqb)**3))&
              /(T%ang(iq,oqb)*T%sqr(i3,i2)*T%ang(oqb,[i3,i2],oqb)&
              *T%ang(i1,[i3,i2],oqb)*T%ang(oqb,[i3,oqb],i2))&
              -((T%ang(i1,i2)**3)*T%sqr(iq,oqb)**2)&
              /(T%ang(i3,i2)*T%ang(i1,[iq,oqb],oqb)*T%ang(i3,oqb,iq)*T%sinv(oqb,iq)) &
              +(T%ang(oqb,[i3,oqb],iq)**3)&
              /(T%ang(i3,oqb)*T%sqr(i1,i2)*T%sqr(iq,i1)&
              *T%ang(oqb,[i3,oqb],i2)*T%ang(i3,oqb,iq)*T%Q(oqb)%kstr)
        else ! g- g- g- qb* q+
          rslt =T%sqr(iq,oqb)**2/(T%sqr(oqb,i3)*T%sqr(i2,i3)*T%sqr(i2,i1)*T%sqr(iq,i1))
        endif
      endif
    endif
  else
    if (helicity(i1).eq.1) then
      if (helicity(i2).eq.1) then
        if (helicity(i3).eq.1) then ! g+ g+ g+ qb* q-
          rslt =T%ang(oqb,iq)**2/(T%ang(i3,oqb)*T%ang(i3,i2)*T%ang(i1,i2)*T%ang(i1,iq))
        else ! g+ g+ g- qb* q-
          rslt = -((T%sqr(oqb,i1)**3)*(T%ang(oqb,i3)**3))&
              /(T%sqr(oqb,iq)*T%ang(i2,i3)*T%ang(oqb,[i2,i3],oqb)&
              *T%ang(oqb,[i2,i3],i1)*T%ang(i2,[oqb,i3],oqb))&
              -((T%sqr(i2,i1)**3)*T%ang(oqb,iq)**2)&
              /(T%sqr(i2,i3)*T%ang(oqb,[oqb,iq],i1)*T%ang(iq,oqb,i3)*T%sinv(oqb,iq)) &
              +(T%ang(iq,[oqb,i3],oqb)**3)&
              /(T%sqr(oqb,i3)*T%ang(i2,i1)*T%ang(i1,iq)&
              *T%ang(i2,[oqb,i3],oqb)*T%ang(iq,oqb,i3)*T%Q(oqb)%kapp)
        endif
      else
        if (helicity(i3).eq.1) then ! g+ g- g+ qb* q-
          rslt = -((T%sqr(oqb,i1)**3)*(T%ang(oqb,i2)**4))&
              /(T%sqr(oqb,iq)*T%ang(oqb,i3)*T%ang(i2,i3)&
              *T%ang(oqb,[i2,i3],oqb)*T%ang(oqb,[i2,i3],i1)*T%ang(i2,[oqb,i3],oqb))&
              -((T%sqr(i3,i1)**4)*(T%ang(oqb,iq)**2))&
              /(T%sqr(i2,i3)*T%sqr(i2,i1)*T%ang(oqb,[oqb,iq],i1)&
              *T%ang(iq,oqb,i3)*T%sinv(oqb,iq))&
              +   T%sqr(oqb,i3)**2 * T%ang(i2,iq)**3 * T%ang(i2,oqb,i3) &
              /( T%ang(i2,i1) * T%ang(i1,iq) * T%ang(i2,[oqb,i3],oqb) &
              * T%ang(iq,oqb,i3) * T%sinv(oqb,i3) * T%Q(oqb)%kapp ) 
        else ! g+ g- g- qb* q-
          rslt =T%sqr(oqb,i1)**3/(T%sqr(oqb,i3)*T%sqr(oqb,iq)*T%sqr(i3,i2)*T%sqr(i1,i2)*T%Q(oqb)%kapp)
        endif
      endif
    else
      if (helicity(i2).eq.1) then
        if (helicity(i3).eq.1) then ! g- g+ g+ qb* q-
          rslt = - T%sqr(i2,iq) * T%ang(oqb,[i1,iq],i2)**2 &
              /( T%sqr(i2,i1) * T%sqr(i1,iq) * T%ang(oqb,i3) * T%ang(i3,oqb,iq) * T%sinv(oqb,i3) ) &
              - T%ang(i1,oqb,iq) * T%ang(i1,[oqb,iq],oqb)**2 &
              /( T%sqr(oqb,iq) * T%ang(i2,i3) * T%ang(i2,i1) * T%ang(i3,oqb,iq) * T%sinv(oqb,iq) * T%Q(oqb)%kapp  )
        else ! g- g+ g- qb* q-
          rslt =(T%sqr(oqb,i2)**3)*T%sqr(i2,iq)/(T%sqr(oqb,i3)*T%sqr(oqb,iq)*T%sqr(i3,i2)&
              *T%sqr(i1,i2)*T%sqr(i1,iq)*T%Q(oqb)%kapp)
        endif
      else
        if (helicity(i3).eq.1) then ! g- g- g+ qb* q-
          rslt =((T%sqr(oqb,i3)**2)*T%sqr(i3,iq))/(T%sqr(oqb,iq)*T%sqr(i3,i2)*T%sqr(i1,i2) &
              *T%sqr(i1,iq)*T%Q(oqb)%kapp)
        else ! g- g- g- qb* q-
          rslt = 0
        endif
      endif
    endif
  endif
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


  subroutine fill_matrices_qq
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




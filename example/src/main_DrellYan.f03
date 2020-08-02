module main_DrellYan
  use amp4hef
  implicit none
  private :: build_momenta_2x3, build_momenta_2x2, build_proton_momenta
  public  :: matrix_element_2x3, matrix_element_2x2

contains


  subroutine matrix_element_2x2(S,  xF, qT, M, xq, kT2, fi_k, ampSquared)
    real(fltknd),intent(in)  :: qT, kT2, fi_k, S, xF, M, xq
    real(fltknd), intent(out) :: ampSquared
    real(fltknd) :: momenta(0:3,8),directions(0:3,5), P1(0:3), P2(0:3), qT_tab(1:2), sqrt_S, amp_algeb, z
    integer :: id1, Ntotal, Noff, NZ

    id1 = 1; Ntotal = 4; Noff = 1; NZ = 1
    call put_process( id1 ,Ntotal ,Noff , NZ, [0, -1, 1, 2] )
    sqrt_S = sqrt(S)
    call build_proton_momenta (sqrt_S, P1, P2)
    qT_tab = [qT*cos(3.14/6.), qT*sin(3.14/6.)]

    call build_momenta_2x2(P1, P2, xF, qT_tab, M, xq, sqrt(kT2)*[cos(fi_k), sin(fi_k)], momenta(0:3,1), &
                   momenta(0:3,2), momenta(0:3,3), momenta(0:3,4))

    directions(0:3,1) = - P1(0:3)
    directions(0:3,2) = - P2(0:3)

    directions(0,Ntotal) = 1
    directions(1,Ntotal) = 0
    directions(2,Ntotal) = 0
    directions(3,Ntotal) = 1
    call put_momenta( id1 ,momenta ,directions )
    call matrix_element_b( id1 ,ampSquared )
    z = xF/xq
  end subroutine


  subroutine matrix_element_2x3(S, x1, x2, k1T2, k2T2, fi_k1, fi_k2,  qT, xF ,M ,z &
                               ,fi_kappa, ampSquared)
    real(fltknd),intent(in)  :: qT, k1T2, k2T2, fi_k1, fi_k2, S, xF, M, x1, x2, z, fi_kappa
    real(fltknd), intent(out) :: ampSquared
    real(fltknd) :: momenta(0:3,8),directions(0:3,5), P1(0:3), P2(0:3), k1T(1:2), k2T(1:2), qT_tab(1:2), sqrt_S
    integer :: id1, Ntotal, Noff, NZ
    k1T(1) = sqrt(k1T2) * cos(fi_k1)
    k1T(2) = sqrt(k1T2) * sin(fi_k1)

    k2T(1) = sqrt(k2T2) * cos(fi_k2)
    k2T(2) = sqrt(k2T2) * sin(fi_k2)
    qT_tab = [0._fltknd, qT]
    sqrt_S = sqrt(S)
    id1 = 1; Ntotal = 5; Noff = 2; NZ = 1
    call put_process( id1 ,Ntotal ,Noff , NZ, [0, 0, 1, -1, 2] )
    call build_proton_momenta (sqrt_S, P1, P2)
    call build_momenta_2x3(P1, P2, xF, qT_tab, M, x1, x2, k1T, k2T, z, fi_kappa, momenta(0:3,1), &
                   momenta(0:3,2), momenta(0:3,3), momenta(0:3,4), momenta(0:3,5))

    directions(0:3,1) = -P1(0:3)
    directions(0:3,2) = -P2(0:3)

    directions(0,Ntotal) = 1
    directions(1,Ntotal) = 0
    directions(2,Ntotal) = 0
    directions(3,Ntotal) = 1

    call put_momenta( id1 ,momenta ,directions )
    call matrix_element_b( id1 ,ampSquared )
  end subroutine


  subroutine build_proton_momenta (sqrt_S, P1, P2)
    real(fltknd), intent(in)  :: sqrt_S
    real(fltknd), intent(out) :: P1(0:3), P2(0:3)

    P1(0:3) = [sqrt_S/2, 0._fltknd, 0._fltknd,  sqrt_S/2]
    P2(0:3) = [sqrt_S/2, 0._fltknd, 0._fltknd, -sqrt_S/2]
  end subroutine



  subroutine build_momenta_2x2(P1, P2, xF, qT, M, xq, kT, k1, pp2, p3, q)
    real(fltknd),intent(in)  :: qT(1:2), kT(1:2), P1(0:3), P2(0:3), xF, M, xq
    real(fltknd),intent(out) :: k1(0:3), pp2(0:3), q(0:3), p3(0:3)
    real(fltknd) :: MT, S, xg, z, qTabs, D1, D2, kTqT, kTY, ab, gd, kTabs


    MT = sqrt( M*M + qT(1)**2 + qT(2)**2 )
    qTabs= sqrt(qT(1)**2 + qT(2)**2 )
    S = (P1(0)+P2(0))**2 - (P1(1)+P2(1))**2 - (P1(2)+P2(2))**2 - (P1(3)+P2(3))**2
    z = xF/xq
    xg = ((1-z)*(M**2)+qT(1)**2+qT(2)**2) &
       /(S*xF*(1-z)) &
       + z*(kT(1)**2+kT(2)**2-2*kT(1)*qT(1)-2*kT(2)*qT(2)) &
       /(S*xF*(1-z))


    k1(0:3) = xg*P1(0:3)
    k1(1:2) = k1(1:2) + kT(1:2)
    pp2(0:3) = xq*P2(0:3)

    q(0:3) = xF*P2(0:3) + (MT**2)*P1(0:3)/(xF*S)
    q(1:2) = q(1:2) + qT(1:2)

    p3(0:3) = k1(0:3) + pp2(0:3) - q(0:3)
    k1(0:3) = -k1(0:3)
    pp2(0:3) = -pp2(0:3)
    kTqT = kT(1)*qT(1)+kT(2)*qT(2)
    kTY = (kT(1)*qT(2)-kT(2)*qT(1))/qTabs
    kTabs= sqrt(kT(1)**2 + kT(2)**2 )
    D2 = MT**2-z*M**2
    D1 = M**2*(1-z) + qTabs**2 + z**2*kTabs**2-2*z*kTqT

    write(*,*) "algebr", (xF-xQ)*S
  end subroutine



  subroutine build_momenta_2x3(P1, P2, xF, qT, M, x1, x2, k1T, k2T, z, fi_kappa, k1, k2, p3, p4, q)
    real(fltknd),intent(in)  :: qT(1:2), k1T(1:2), k2T(1:2), P1(0:3), P2(0:3), xF, M, x1, x2, z, fi_kappa
    real(fltknd),intent(out) :: k1(0:3), k2(0:3), q(0:3), p3(0:3), p4(0:3)
    real(fltknd) :: MT, S, xqq, p3T(1:2), p4T(1:2), delt(1:2), kapp(1:2)

    k1(0:3) = x1*P1(0:3)
    k1(1:2) = k1(1:2) + k1T(1:2)

    k2(0:3) = x2*P2(0:3)
    k2(1:2) = k2(1:2) + k2T(1:2)

    MT = sqrt( M*M + qT(1)**2 + qT(2)**2 )
    S = (P1(0)+P2(0))**2 - (P1(1)+P2(1))**2 - (P1(2)+P2(2))**2 - (P1(3)+P2(3))**2

    q(0:3) = xF*P1(0:3) + MT**2/(xF*S)*P2(0:3)
    q(1:2) = q(1:2) +  qT(1:2)

    xqq = x1 - xF
    delt(1:2) = k1T(1:2) + k2T(1:2) -qT(1:2)
    kapp(1) = sqrt(abs(z*(1-z)*(xqq*x2*S-xqq*MT**2/xF-delt(1)**2 - delt(2)**2)))*cos(fi_kappa)
    kapp(2) = sqrt(abs(z*(1-z)*(xqq*x2*S-xqq*MT**2/xF-delt(1)**2 - delt(2)**2)))*sin(fi_kappa)
!    write(*,*) "kapp:", kapp
!    write(*,*) "un k:", z*(1-z)*(xqq*x2*S-xqq*MT**2/xF-delt(1)**2 - delt(2)**2)
    p3T(1:2) = z*delt(1:2) + kapp(1:2)
    p4T(1:2) = (1-z)*delt(1:2) - kapp(1:2)

    p3(0:3) = z*xqq*P1(0:3) + (p3T(1)**2 + p3T(2)**2)/(z*xqq*S)*P2(0:3)

    p3(1:2) = p3(1:2) + p3T(1:2)

    p4(0:3) = (1-z)*xqq*P1(0:3) + (p4T(1)**2 + p4T(2)**2)/((1-z)*xqq*S)*P2(0:3)
    p4(1:2) = p4(1:2) + p4T(1:2)

    k1(0:3) = -k1(0:3)
    k2(0:3) = -k2(0:3)

  end subroutine




end module




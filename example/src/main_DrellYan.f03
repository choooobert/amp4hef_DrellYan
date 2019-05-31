module main_DrellYan
  use amp4hef
  implicit none
  private :: build_momenta_2x3, build_momenta_2x2, build_proton_momenta
  public  :: matrix_element_2x3, matrix_element_2x2

contains

  subroutine matrix_element_2x3(sqrt_S, Y, qT, M, x1, x2, k1T, k2T, z, fi, ampSquared)
    real(fltknd),intent(in)  :: qT(1:2), k1T(1:2), k2T(1:2), sqrt_S, Y, M, x1, x2, z, fi
    real(fltknd), intent(out) :: ampSquared
    real(fltknd) :: momenta(0:3,8),directions(0:3,5), P1(0:3), P2(0:3)
    integer :: id1, Ntotal, Noff, NZ

    id1 = 1; Ntotal = 5; Noff = 2; NZ = 1
    call put_process( id1 ,Ntotal ,Noff , NZ, [0, 0, 1, -1, 2] )
    call build_proton_momenta (sqrt_S, P1, P2)
    call build_momenta_2x3(P1, P2, Y, qT, M, x1, x2, k1T, k2T, z, fi, momenta(0:3,1), &
                   momenta(0:3,2), momenta(0:3,3), momenta(0:3,4), momenta(0:3,5))

    directions(0,1:2) = momenta(0,1:2)
    directions(1,1:2) = 0
    directions(2,1:2) = 0
    directions(3,1:2) = momenta(3,1:2)

    directions(0,Ntotal) = 1
    directions(1,Ntotal) = 0
    directions(2,Ntotal) = 0
    directions(3,Ntotal) = 1

    call put_momenta( id1 ,momenta ,directions )
    call matrix_element_b( id1 ,ampSquared )
  end subroutine


  subroutine matrix_element_2x2(sqrt_S,  Y, qT, M, xq, kT, ampSquared)
    real(fltknd),intent(in)  :: qT(1:2), kT(1:2), sqrt_S, Y, M, xq
    real(fltknd), intent(out) :: ampSquared
    real(fltknd) :: momenta(0:3,8),directions(0:3,5), P1(0:3), P2(0:3)
    integer :: id1, Ntotal, Noff, NZ

    id1 = 2; Ntotal = 4; Noff = 1; NZ = 1
    call put_process( id1 ,Ntotal ,Noff , NZ, [0, -1, 1, 2] )
    call build_proton_momenta (sqrt_S, P1, P2)
    call build_momenta_2x2(P1, P2, Y, qT, M, xq, kT, momenta(0:3,1), &
                   momenta(0:3,2), momenta(0:3,3), momenta(0:3,4))

    directions(0,1) = momenta(0,1)
    directions(1,1) = 0
    directions(2,1) = 0
    directions(3,1) = momenta(3,1)

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


  subroutine build_momenta_2x3(P1, P2, Y, qT, M, x1, x2, k1T, k2T, z, fi, k1, k2, p3, p4, q)
    real(fltknd),intent(in)  :: qT(1:2), k1T(1:2), k2T(1:2), P1(0:3), P2(0:3), Y, M, x1, x2, z, fi
    real(fltknd),intent(out) :: k1(0:3), k2(0:3), q(0:3), p3(0:3), p4(0:3)
    real(fltknd) :: MT, S, xF, xqq, p3T(1:2), p4T(1:2), delt(1:2), kapp(1:2)

    k1(0:3) = x1*P1(0:3)
    k1(1:2) = k1(1:2) + k1T(1:2)
    k1(0:3) = -k1(0:3)

    k2(0:3) = x2*P2(0:3)
    k2(1:2) = k2(1:2) + k2T(1:2)
    k2(0:3) = -k2(0:3)

    MT = sqrt( M*M + qT(1)**2 + qT(2)**2 )
    S = (P1(0)+P2(0))**2 - (P1(1)+P2(1))**2 - (P1(2)+P2(2))**2 - (P1(3)+P2(3))**2
    xF = MT/sqrt(S) * exp(Y)

    q(0:3) = xF*P1(0:3) + MT**2/(xF*S)* P2(0:3)
    q(1:2) = q(1:2) +  qT(1:2)

    write(*,*) sqrt(q(0)**2 - q(1)**2 - q(2)**2 -q(3)**2)

    xqq = x1 - xF
    delt(1:2) = k1T(1:2) + k2T(1:2) -qT(1:2)
    kapp(1) = sqrt(z*(1-z)*(xqq*x2*S-xqq*MT**2/xF-delt(1)**2 - delt(2)**2))*cos(fi)
    kapp(2) = sqrt(z*(1-z)*(xqq*x2*S-xqq*MT**2/xF-delt(1)**2 - delt(2)**2))*sin(fi)
    p3T(1:2) = z*delt(1:2) + kapp(1:2)
    p4T(1:2) = (1-z)*delt(1:2) - kapp(1:2)

    p3(0:3) = z*xqq*P1(0:3) + (p3T(1)**2 + p3T(2)**2)/(z*xqq*S)*P2(0:3)
    p3(1:2) = p3(1:2) + p3T(1:2)

    p4(0:3) = (1-z)*xqq*P1(0:3) + (p4T(1)**2 + p4T(2)**2)/((1-z)*xqq*S)*P2(0:3)
    p4(1:2) = p4(1:2) + p4T(1:2)
  end subroutine


  subroutine build_momenta_2x2(P1, P2, Y, qT, M, xq, kT, k1, pp2, p3, q)
    real(fltknd),intent(in)  :: qT(1:2), kT(1:2), P1(0:3), P2(0:3), Y, M, xq
    real(fltknd),intent(out) :: k1(0:3), pp2(0:3), q(0:3), p3(0:3)
    real(fltknd) :: MT, S, xF, xg, z


    MT = sqrt( M*M + qT(1)**2 + qT(2)**2 )
    S = (P1(0)+P2(0))**2 - (P1(1)+P2(1))**2 - (P1(2)+P2(2))**2 - (P1(3)+P2(3))**2
    xF = MT/sqrt(S) * exp(Y)
    z = xF/xq
    xg = ((1-z)*(M**2)+qT(1)**2+qT(2)**2)/(S*xF*(1-z))

    k1(0:3) = xg*P1(0:3)
    k1(1:2) = k1(1:2) + kT(1:2)
    k1(0:3) = -k1(0:3)

    pp2(0:3) = xq*P2(0:3)
    pp2(0:3) = -pp2(0:3)

    q(0:3) = xF*P1(0:3) + MT**2/(xF*S)* P2(0:3)
    q(1:2) = q(1:2) +  qT(1:2)

    p3(0:3) = -k1(0:3) -pp2(0:3) -q(0:3)
!    write(*,*) "p3  ", p3(0:3)
!    write(*,*) "p3^2", p3(0)**2 - p3(1)**2 - p3(2)**2 -p3(3)**2
!    write(*,*) " q^2", sqrt(q(0)**2 - q(1)**2 - q(2)**2 -q(3)**2)

  end subroutine

end module




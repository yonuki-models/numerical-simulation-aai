module Floquet_Governing_Equation
  implicit none

  integer, parameter :: NZ = 2

contains

  subroutine equation(N, Omega, D, M, K, a, b, term)
    implicit none
    double precision, intent(in) :: N
    double precision, intent(in) :: Omega(1:3)
    double precision, intent(in) :: D(1:3,1:3)
    double precision, intent(in) :: M(1:3)
    double precision, intent(in) :: K(1:3)
    double precision, intent(in) :: a(1:3)
    double precision, intent(in) :: b
    double precision, intent(out) :: term(1:4)

    double precision :: kappa_square
    double precision :: tmp

  ! Levi-Civita symbol
    integer :: j_plus(1:3) = 0
    integer :: k_plus(1:3) = 0
    integer :: j_minus(1:3) = 0
    integer :: k_minus(1:3) = 0

    integer :: i
    integer :: j

    j_plus(1) = 2
    j_plus(2) = 3
    j_plus(3) = 1

    k_plus(1) = 3
    k_plus(2) = 1
    k_plus(3) = 2

    j_minus(1) = 3
    j_minus(2) = 1
    j_minus(3) = 2

    k_minus(1) = 2
    k_minus(2) = 3
    k_minus(3) = 1

    term(:) = 0.d0

    kappa_square = 0.d0
    do i = 1, 3
      kappa_square = kappa_square + K(i)**2
    enddo

    do i = 1, 3
      do j = 1, 3
        term(i) = term(i) - a(j) * D(j,i)
      enddo
    enddo

    do i = 1, 3
      term(i) = term(i) - 2 * (  &
        &  Omega(j_plus(i)) * a(k_plus(i))  &
        &  - Omega(j_minus(i)) * a(k_minus(i)) )
    enddo

    term(3) = term(3) + N * b

! Pressure terms
    tmp = 0.d0
    do i = 1, 3
      do j = 1, 3
        tmp = tmp + 2 * a(i) * D(i,j) * K(j)
      enddo
    enddo

    do i = 1, 3
      tmp = tmp + 2 * (  &
        &  Omega(j_plus(i)) * a(k_plus(i))  &
        &  - Omega(j_minus(i)) * a(k_minus(i)) )  &
        &  * K(i)
    enddo

    tmp = tmp - N * b * K(3)

    do i = 1, 3
      term(i) = term(i) + K(i) / kappa_square * tmp
    enddo

! Buoyancy
    do i = 1, 3
      term(4) = term(4) - a(i) * M(i)
    enddo
    term(4) = term(4) - N * a(3)

    return
  end subroutine equation

end module Floquet_Governing_Equation

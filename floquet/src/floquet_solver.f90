module Floquet_Solver
  implicit none
  private

  integer, parameter :: NZ = 4
  double precision, parameter :: K0 = 1.d0

  public elliptic_instability_solver

contains

  subroutine elliptic_instability_solver(  &
    &  alpha_plus, alpha_minus, K3,  &
    &  growth_rate, unstable_eigen_vector, neutral_eigen_vector)
    use Floquet_Control, only : N => buoyancy_frequency
    use Floquet_Control, only : Omega => angular_velocity
    use Floquet_Control, only : N_time
    use Floquet_Eigen_Problem, only : eigen_solver
    use Floquet_Governing_Equation, only : equation
    implicit none
    double precision, intent(in) :: alpha_plus
    double precision, intent(in) :: alpha_minus
    double precision, intent(in) :: K3
    double precision, intent(out) :: growth_rate
    double precision, intent(out) :: unstable_eigen_vector(1:NZ)
    double precision, intent(out) :: neutral_eigen_vector(1:NZ)

    complex(kind(0d0)) :: lambda(1:NZ)
    complex(kind(0d0)) :: eigen_zeta(1:NZ, 1:NZ)
    complex(kind(0d0)) :: unstable_C(1:NZ)
    complex(kind(0d0)) :: neutral_C(1:NZ)

    complex(kind(0d0)) :: tmp(1:NZ)

    double precision :: lambda_abs(1:NZ)

    double precision, parameter :: PI = acos(-1.d0)
    complex(kind(0d0)) :: UI = (0.d0, 1.d0)
    double precision, parameter :: incompressible_tolerance = 1.d-5

    double precision :: current_time = 0.d0
    double precision :: wavenumber_frequency
    double precision :: end_time
    double precision :: tmp_time = 0.d0
    double precision :: delta_time

    double precision :: angle

    double precision :: D(1:3, 1:3)
    double precision :: M(1:3)
    double precision :: K(1:3)
    double precision :: K_tmp(1:3)

    double precision :: e

    complex(kind(0d0)) :: d2
    complex(kind(0d0)) :: d3

    integer :: sorted_index(1:NZ)

    integer :: it
    integer :: i
    integer :: j

    double precision :: zeta(1:NZ, 1:NZ) = 0.d0
    double precision :: zeta_tmp(1:NZ) = 0.d0
    double precision :: term(1:NZ, 4) = 0.d0

    logical :: exist_NAN = .false.

    D(1:3, 1:3) = 0.d0
    D(2,1) = alpha_plus
    D(1,2) = - alpha_minus
    M(1:3) = 0.d0
    K(3) = K3

    wavenumber_frequency = sqrt(alpha_minus * alpha_plus)

    end_time = 2 * PI / wavenumber_frequency
    delta_time = end_time / N_time
    current_time = 0.d0

    e = sign(sqrt(alpha_plus / alpha_minus), alpha_plus)

    zeta(:,:) = (0.d0, 0.d0)
    do i = 1, NZ
      zeta(i,i) = (1.d0, 0.d0)
    enddo
    K(1) = K0
    K(2) = 0.d0

    do it = 0, N_time - 1
      do i = 1, NZ
        zeta_tmp(:) = zeta(:,i)

        tmp_time = current_time
        K_tmp(1) = K0 * cos(wavenumber_frequency * tmp_time)
        K_tmp(2) = - e * K0 * sin(wavenumber_frequency * tmp_time)
        K_tmp(3) = K3
        call equation(N, Omega, D, M, K_tmp, zeta_tmp(1:3), zeta_tmp(4),  &
          &  term(:,1))

        zeta_tmp(:) = zeta(:,i) + term(:,1) * delta_time * 0.5
        tmp_time = current_time + delta_time * 0.5
        K_tmp(1) = K0 * cos(wavenumber_frequency * tmp_time)
        K_tmp(2) = - e * K0 * sin(wavenumber_frequency * tmp_time)
        K_tmp(3) = K3
        call equation(N, Omega, D, M, K_tmp, zeta_tmp(1:3), zeta_tmp(4),  &
          &  term(:,2))

        zeta_tmp(:) = zeta(:,i) + term(:,2) * delta_time * 0.5
        tmp_time = current_time + delta_time * 0.5
        K_tmp(1) = K0 * cos(wavenumber_frequency * tmp_time)
        K_tmp(2) = - e * K0 * sin(wavenumber_frequency * tmp_time)
        K_tmp(3) = K3
        call equation(N, Omega, D, M, K_tmp, zeta_tmp(1:3), zeta_tmp(4),  &
          &  term(:,3))

        zeta_tmp(:) = zeta(:,i) + term(:,3) * delta_time
        tmp_time = current_time + delta_time
        K_tmp(1) = K0 * cos(wavenumber_frequency * tmp_time)
        K_tmp(2) = - e * K0 * sin(wavenumber_frequency * tmp_time)
        K_tmp(3) = K3
        call equation(N, Omega, D, M, K_tmp, zeta_tmp(1:3), zeta_tmp(4),  &
          &  term(:,4))

        zeta(:,i) = zeta(:,i)  &
          &  + (term(:,1) + 2 * term(:,2) + 2 * term(:,3) + term(:,4))  &
          &  * delta_time / 6

      end do
      current_time = current_time + delta_time
      K(:) = K_tmp(:)
    end do

    exist_NAN = .false.
    do i = 1, NZ
      do j = 1, NZ
        if (zeta(i,j) .ne. zeta(i,j)) exist_NAN = .true.
      end do
    end do

    if (exist_NAN) then
      lambda(:) = (1.d0, 0.d0)
      eigen_zeta(:,:) = (0.d0, 0.d0)
    else
      call eigen_solver(4, zeta, lambda, eigen_zeta)
    end if

    call sort_eigenvalues(lambda, sorted_index)

    do i = 1, NZ
      lambda_abs(i) = abs(lambda(sorted_index(i)))
      growth_rate = log(abs(lambda(sorted_index(i)))) / end_time
    enddo
    growth_rate = log(lambda_abs(1)) / end_time

    K(1) = K0
    K(2) = 0.d0
    d2 = K(1) * eigen_zeta(1, sorted_index(2))  &
      &  + K(2) * eigen_zeta(2, (sorted_index(2))) &
      &  + K(3) * eigen_zeta(3, (sorted_index(2)))
    d3 = K(1) * eigen_zeta(1, sorted_index(3))  &
      &  + K(2) * eigen_zeta(2, (sorted_index(3))) &
      &  + K(3) * eigen_zeta(3, (sorted_index(3)))

    unstable_C(1:NZ) = eigen_zeta(:, sorted_index(1))
    neutral_C(1:NZ) = d3 * eigen_zeta(:, sorted_index(2))  &
      &  - d2 * eigen_zeta(:, sorted_index(3))

    unstable_C(1:NZ) = unstable_C(1:NZ)  &
      &  / sum(abs(unstable_C(1:NZ))**2)**0.5  &
      &  * (dble(unstable_C(4))  &
      &  - UI * dimag(unstable_C(4)))  &
      &  / abs(unstable_C(4))
    unstable_eigen_vector(1:NZ) = dble(unstable_C(1:NZ))

    neutral_C(1:NZ) = neutral_C(1:NZ)  &
      &  / sum(abs(neutral_C(1:NZ))**2)**0.5  &
      &  * (dble(neutral_C(4))  &
      &  - UI * dimag(neutral_C(4)))  &
      &  / abs(neutral_C(4))
    neutral_eigen_vector(1:NZ) = dble(neutral_C(1:NZ))

    return
  end subroutine elliptic_instability_solver

! *****************************************************************************

  subroutine calculate_wavenumber_tmp(K_in, D, delta_time, time_factor, K_out)
    implicit none
    double precision, intent(in) :: K_in(1:3)
    double precision, intent(in) :: D(1:3, 1:3)
    double precision, intent(in) :: delta_time
    double precision, intent(in) :: time_factor
    double precision, intent(out) :: K_out(1:3)

    double precision :: K_term(1:3)

    integer :: i
    integer :: j

    K_term(:) = 0.d0
    do i = 1, 3
      do j = 1, 3
        K_term(i) = K_term(i) - D(i,j) * K_in(j)
      enddo
    enddo

    K_out(:) = K_in(:) - K_term(:) * delta_time * time_factor

    return
  end subroutine calculate_wavenumber_tmp

! *****************************************************************************

  subroutine increment_wavenumber(K_in, D, delta_time, K_out)
    implicit none
    double precision, intent(in) :: K_in(1:3)
    double precision, intent(in) :: D(1:3, 1:3)
    double precision, intent(in) :: delta_time
    double precision, intent(out) :: K_out(1:3)

    double precision :: K_term(1:3, 1:4)

    integer :: i
    integer :: j

    K_term(:, :) = 0.d0

    K_term(i, 1) = - D(i,j) * K_in(j)

    do i = 1, 3
      do j = 1, 3
        K_term(i, 2) = K_term(i, 2) - D(i,j) * K_in(j) * 0.5d0
      enddo
    enddo

    do i = 1, 3
      do j = 1, 3
        K_term(i, 3) = K_term(i, 3) - D(i,j) * K_in(j) * 0.5d0
      enddo
    enddo

    do i = 1, 3
      do j = 1, 3
        K_term(i, 4) = K_term(i, 4) - D(i,j) * K_in(j)
      enddo
    enddo

    K_out(:) = K_out(:)  &
      &  + (K_term(:,1) + 2 * K_term(:,2) + 2 * K_term(:,3) + K_term(:,4))  &
      &  * delta_time / 6

    return
  end subroutine increment_wavenumber

! *****************************************************************************

  subroutine sort_eigenvalues(lambda, sorted_index)
    implicit none
    complex(kind(0d0)), intent(in) :: lambda(1:NZ)
    integer, intent(out) :: sorted_index(1:NZ)

    double precision :: lambda_abs(1:NZ)
    integer :: index_tmp
    integer :: i
    integer :: j

    do i = 1, NZ
      lambda_abs(i) = abs(lambda(i))
      sorted_index(i) = i
    enddo

    do i = 1, NZ-1
      do j = i+1, NZ
        if (lambda_abs(sorted_index(i)) < lambda_abs(sorted_index(j))) then
          index_tmp = sorted_index(j)
          sorted_index(j) = sorted_index(i)
          sorted_index(i) = index_tmp
        endif
      enddo
    enddo

    return
  end subroutine sort_eigenvalues

end module Floquet_Solver

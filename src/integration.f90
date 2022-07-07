module Integration
  implicit none
  private

  public runge_kutta_3
  public CFL_condition

contains

! *****************************************************************************

  subroutine runge_kutta_3(delta_time)
    use parallel, only : my_rank
    use Geometry, only : dealias_truncate
    use Geometry, only : calculate_wavenumbers
    use Geometry, only : incompressible_projection
    use Geometry, only : state
    use Geometry, only : spectrum_HV
    use Geometry, only : n1C_local, n2C_local, n3
    use Time_Module, only : current_time
    use Time_Module, only : evaluate_delta_time
    use External_Field, only : derive_transform_tensor
    use External_Field, only : derive_velocity_gradient
    use External_Field, only : derive_temperature_gradient
    use Dynamics, only : linear_terms
    use Dynamics, only : nonlinear_terms
    use Dynamics, only : wave_vortex_production
    use Dynamics, only : derive_viscous_factor
    use Dynamics, only : viscous_parameterization
    use Analysis, only : initialize_index
    use Analysis, only : calc_flux_convergence
    use Log, only : log_energy
    implicit none
    double precision, parameter :: FACTOR_A = 8.0d0 / 15
    double precision, parameter :: FACTOR_B = 5.0d0 / 12
    double precision, parameter :: FACTOR_C = 17.0d0 / 60
    double precision, parameter :: FACTOR_D = 3.0d0 / 4
    double precision, parameter :: FACTOR_E = 5.0d0 / 12

    double precision, parameter :: FACTOR_A_RCP = 1.0d0 / FACTOR_A
    double precision, parameter :: FACTOR_B_RCP = 1.0d0 / FACTOR_B
    double precision, parameter :: FACTOR_ABC = FACTOR_A + FACTOR_B - FACTOR_C
                        ! = 2.0d0 / 3.0

    double precision, parameter :: FACTOR_ENERGY_1 = 0.25d0
    double precision, parameter :: FACTOR_ENERGY_2 = 0.0d0
    double precision, parameter :: FACTOR_ENERGY_3 = 0.75d0

    double precision, intent(out) :: delta_time

    double precision :: time_a
    double precision :: time_b
    double precision :: delta_time_CFL
    double precision :: delta_time_until_b

    double precision :: A_a(1:3, 1:3)
    double precision :: A_b(1:3, 1:3)
    double precision :: D(1:3, 1:3)
    double precision :: M(1:3)

    complex(kind(0d0)), pointer :: U(:,:,:,:)
    complex(kind(0d0)), pointer :: T(:,:,:)
    double precision, pointer :: HV_KFC(:,:)
    double precision, pointer :: HV_PFC(:,:)
    double precision :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: K_a(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: K_b(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision :: viscous_factor_a(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: viscous_factor_b(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision ::  &
      &  diffusive_factor_a(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision ::  &
      &  diffusive_factor_b(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    complex(kind(0d0)) :: term_U_a(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: term_T_a(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    complex(kind(0d0)) :: term_U_b(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: term_T_b(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    complex(kind(0d0)) :: PK(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    integer :: i

    U => state%velocity(:,:,:,:)
    T => state%temperature(:,:,:)

    viscous_factor_a(:,:,:) = 1.d0
    viscous_factor_b(:,:,:) = 1.d0
    diffusive_factor_a(:,:,:) = 1.d0
    diffusive_factor_b(:,:,:) = 1.d0

    call derive_transform_tensor(current_time, A_a)
    call calculate_wavenumbers(A_a, K)
    delta_time_CFL = CFL_condition(U, K)
    delta_time = evaluate_delta_time(delta_time_CFL)

    call viscous_parameterization(K, U, T)
    call initialize_index(K)

! Runge-Kutta step 1
    time_a = current_time
    delta_time_until_b = delta_time * FACTOR_A
    time_b = current_time + delta_time_until_b

    call derive_transform_tensor(time_a, A_a)
    call calculate_wavenumbers(A_a, K_a)
    call derive_transform_tensor(time_b, A_b)
    call calculate_wavenumbers(A_b, K_b)

    call derive_velocity_gradient(time_a, D)
    call derive_temperature_gradient(time_a, M)

    call linear_terms(K_a, D, M, delta_time, FACTOR_ENERGY_1,  &
    &  term_U_a, term_T_a)

    call nonlinear_terms(K_a, delta_time, FACTOR_ENERGY_1,  &
    &  term_U_a, term_T_a)

    call wave_vortex_production(K_a, D, delta_time, FACTOR_ENERGY_1)

    do i = 1, 3
      U(i,:,:,:) = U(i,:,:,:) + term_U_a(i,:,:,:) * FACTOR_A
    enddo
    T(:,:,:) = T(:,:,:) + term_T_a(:,:,:) * FACTOR_A

    call derive_viscous_factor(delta_time_until_b,  &
      &  viscous_factor_a, diffusive_factor_a)

    call incompressible_projection(K_a, K_b, U, T, PK,  &
      &  viscous_factor_a, diffusive_factor_a)

    term_U_a(:,:,:,:) = term_U_a(:,:,:,:) - PK(:,:,:,:) * FACTOR_A_RCP

    call dealias_truncate(U, T)

! Runge-Kutta step 2
    time_a = time_b
    delta_time_until_b = delta_time * FACTOR_ABC
    time_b = current_time + delta_time_until_b

    call derive_transform_tensor(time_a, A_a)
    call calculate_wavenumbers(A_a, K_a)
    call derive_transform_tensor(time_b, A_b)
    call calculate_wavenumbers(A_b, K_b)

    call derive_velocity_gradient(time_a, D)
    call derive_temperature_gradient(time_a, M)

    call linear_terms(K_a, D, M, delta_time, FACTOR_ENERGY_2,  &
    &  term_U_b, term_T_b)

    call nonlinear_terms(K_a, delta_time, FACTOR_ENERGY_2,  &
      &  term_U_b, term_T_b)

    call wave_vortex_production(K_a, D, delta_time, FACTOR_ENERGY_2)

    do i = 1, 3
      U(i,:,:,:) = viscous_factor_a(:,:,:)  &
        &  * (U(i,:,:,:) + term_U_b(i,:,:,:) * FACTOR_B)  &
        &  - term_U_a(i,:,:,:) * FACTOR_C
    enddo
    T(:,:,:) = diffusive_factor_a(:,:,:)  &
      &  * (T(:,:,:) + term_T_b(:,:,:) * FACTOR_B)  &
      &  - term_T_a(:,:,:) * FACTOR_C

    call derive_viscous_factor(delta_time_until_b,  &
      &  viscous_factor_b, diffusive_factor_b)

    call incompressible_projection(K_a, K_b, U, T, PK,  &
      &  viscous_factor_b, diffusive_factor_b)

    do i = 1, 3
      term_U_b(i,:,:,:) = viscous_factor_a(:,:,:) * term_U_b(i,:,:,:)  &
        &  - PK(i,:,:,:) * FACTOR_B_RCP
    enddo

    term_T_b(:,:,:) = diffusive_factor_a(:,:,:) * term_T_b(:,:,:)

    call dealias_truncate(U, T)

! Runge-Kutta step 3
    time_a = time_b
    delta_time_until_b = delta_time
    time_b = current_time + delta_time_until_b

    call derive_transform_tensor(time_a, A_a)
    call calculate_wavenumbers(A_a, K_a)
    call derive_transform_tensor(time_b, A_b)
    call calculate_wavenumbers(A_b, K_b)

    call derive_velocity_gradient(time_a, D)
    call derive_temperature_gradient(time_a, M)

    call linear_terms(K_a, D, M, delta_time, FACTOR_ENERGY_3,  &
      &  term_U_a, term_T_a)

    call nonlinear_terms(K_a, delta_time, FACTOR_ENERGY_3,  &
      &  term_U_a, term_T_a)

    call wave_vortex_production(K_a, D, delta_time, FACTOR_ENERGY_3)

    do i = 1, 3
      U(i,:,:,:) = viscous_factor_b(:,:,:)  &
        &  * (U(i,:,:,:) + term_U_a(i,:,:,:) * FACTOR_D)  &
        &  - term_U_b(i,:,:,:) * FACTOR_E
    enddo
    T(:,:,:) = diffusive_factor_b(:,:,:)  &
      &  * (T(:,:,:) + term_T_a(:,:,:) * FACTOR_D)  &
      &  - term_T_b(:,:,:) * FACTOR_E

    call derive_viscous_factor(delta_time_until_b,  &
      &  viscous_factor_a, diffusive_factor_a)

    call incompressible_projection(K_a, K_b, U, T, PK,  &
      &  viscous_factor_a, diffusive_factor_a)

    call dealias_truncate(U, T)

    HV_KFC => spectrum_HV%kinetic_flux_convergence
    HV_PFC => spectrum_HV%potential_flux_convergence
    call calc_flux_convergence(U, T, K_b, HV_KFC, HV_PFC)

    return
  end subroutine runge_kutta_3

! *****************************************************************************

  function CFL_condition(U, K)
    use mpi, only : MPI_COMM_WORLD
    use mpi, only : MPI_DOUBLE_PRECISION
    use mpi, only : MPI_MAX
    use mpi, only : MPI_SUM
    use parallel, only : my_rank
    use FFTW_Interface, only : fftw_backward_parallel
    use Geometry, only : get_mask
    use Geometry, only : n1C_local, n2C_local, n3
    use Geometry, only : n1, n2R_local, n3R_local
    use Geometry, only : n2
    use Geometry, only : n1_truncated, n2_truncated, n3_truncated
    implicit none
    complex(kind(0d0)), intent(in) ::  &
      &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) ::  &
      &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: UC(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: UR(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    integer :: mask(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: U_max_local(1:3)
    double precision :: K_max_local(1:3)
    double precision :: U_max(1:3)
    double precision :: K_max(1:3)
    double precision :: UK_max
    double precision :: CFL_factor = 1.83d0
    double precision :: CFL_condition
    integer :: i

    integer :: ierr

    call get_mask(mask)

    do i = 1, 3
      UC(:,:,:) = U(i,:,:,:)
      call fftw_backward_parallel(UC(:,:,:), UR(:,:,:))
      U_max_local(i) = maxval(abs(UR(:,:,:)))
      K_max_local(i) = maxval(abs(K(i,:,:,:) * mask(:,:,:)))
    enddo

    call MPI_Allreduce(U_max_local(1), U_max(1), 3, MPI_DOUBLE_PRECISION,  &
      &  MPI_MAX, MPI_COMM_WORLD, ierr)

    call MPI_Allreduce(K_max_local(1), K_max(1), 3, MPI_DOUBLE_PRECISION,  &
      &  MPI_MAX, MPI_COMM_WORLD, ierr)

    UK_max = maxval(U_max(:) * K_max(:))

    CFL_condition = CFL_factor / UK_max

    return
  end function CFL_condition

end module Integration

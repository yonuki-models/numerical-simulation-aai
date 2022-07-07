module Initial_Condition
  use Parallel, only : rank_1
  use Parallel, only : rank_2
  use Geometry, only : n1C_local, n2C_local, n3
  use Geometry, only : n1_truncated, n2_truncated, n3_truncated
  use Parameters, only : initial_noise_level
  implicit none
  private

  public :: white_noise
  public :: white_noise_low_modes
  public :: interaction_test
  public :: preconditioning
  public :: eliminate_vortical_mode
  public :: adjust_to_power_law

contains

! *****************************************************************************

  subroutine white_noise()
    use Parallel, only : my_rank
    use Geometry, only : domain_length
    use Geometry, only : state
    implicit none
    double precision, parameter :: PI = acos(-1.d0)
    complex(kind(0d0)), parameter :: UI = (0.d0, 1.d0)
    double precision ::  &
      &  homogeneous(1:6, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: gauss(1:6, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    integer, allocatable:: seed(:)
    integer :: is, seedsize

    integer :: i
    integer :: i1
    integer :: i2
    integer :: i3

    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    do is = 1, seedsize
      call system_clock(count=seed(is))
      seed(is) = seed(is) * (my_rank + 1)
    end do
    call random_seed(put=seed(:))

    call random_number(homogeneous)
    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local - 1
        do i1 = 0, n1C_local - 1
          do i = 1, 6
            if (homogeneous(i, i1, i2, i3) == 0.d0) then
              homogeneous(i, i1, i2, i3) = 0.5d0
            endif
          enddo
        enddo
      enddo
    enddo

    gauss(1:3,:,:,:) = sqrt(-2 * log(homogeneous(1:3,:,:,:)))  &
    &  * cos(2 * PI * homogeneous(4:6,:,:,:))
    gauss(4:6,:,:,:) = sqrt(-2 * log(homogeneous(1:3,:,:,:)))  &
    &  * sin(2 * PI * homogeneous(4:6,:,:,:))

    state%velocity(1:3,:,:,:) = initial_noise_level  &
    &  * (gauss(1:3,:,:,:) + UI * gauss(4:6,:,:,:))

    call preconditioning()

    return
  end subroutine white_noise

! *****************************************************************************

  subroutine white_noise_low_modes()
    use Parallel, only : my_rank
    use Geometry, only : domain_length
    use Geometry, only : state
    use Geometry, only : n2
    use Input_Output, only : print_main
    implicit none
    double precision, parameter :: PI = acos(-1.d0)
    complex(kind(0d0)), parameter :: UI = (0.d0, 1.d0)
    double precision ::  &
      &  homogeneous(1:6, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: gauss(1:6, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision :: tmp
    integer, allocatable:: seed(:)
    integer :: mode_truncate
    integer :: is, seedsize
    integer :: i
    integer :: i1_global
    integer :: i2_global

    integer :: i1
    integer :: i2
    integer :: i3

    call random_seed(size=seedsize)

    allocate(seed(seedsize))

    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local - 1
        i2_global = rank_2 * n2C_local + i2
        do i1 = 0, n1C_local - 1
          i1_global = rank_1 * n1C_local + i1
          do is = 1, seedsize
            call system_clock(count=seed(is))
              ! To fix random seed
              ! seed(is) = i1_global * 37**2  &
              ! &  + sign(min(i2_global, n2 - i2_global), n2/2 - i2_global)  &
              ! &  + sign(min(i3, n3 - i3), n3/2 - i3) * 37 + is
          end do
          call random_seed(put=seed(:))
          call random_number(homogeneous(1:6, i1, i2, i3))
          do i = 1, 6
            if (homogeneous(i, i1, i2, i3) == 0.d0) then
              homogeneous(i, i1, i2, i3) = 0.5d0
            endif
          enddo
        enddo
      enddo
    enddo

    gauss(1:3,:,:,:) = sqrt(-2 * log(homogeneous(1:3,:,:,:)))  &
    &  * cos(2 * PI * homogeneous(4:6,:,:,:))
    gauss(4:6,:,:,:) = sqrt(-2 * log(homogeneous(1:3,:,:,:)))  &
    &  * sin(2 * PI * homogeneous(4:6,:,:,:))

    state%velocity(1:3,:,:,:) = (gauss(1:3,:,:,:) + UI * gauss(4:6,:,:,:))

    call preconditioning()
    call eliminate_vortical_mode()

    mode_truncate = 1
    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local-1
        i2_global = rank_2 * n2C_local + i2
        do i1 = 0, n1C_local-1
          i1_global = rank_1 * n1C_local + i1
          if (i1_global <= mode_truncate  &
          &  .and. (i2_global <= mode_truncate  &
          &  .or. i2_global >= n2-mode_truncate)  &
          &  .and. (i3 <= mode_truncate  &
          &  .or. i3 >= n3-mode_truncate)) then
            tmp = sum(abs(state%velocity(1:3, i1, i2, i3))**2)  &
            &  + abs(state%temperature(i1, i2, i3))**2
            state%velocity(1:3, i1, i2, i3)  &
            &  = state%velocity(1:3, i1, i2, i3) / sqrt(tmp)
            state%temperature(i1, i2, i3)  &
            &  = state%temperature(i1, i2, i3) / sqrt(tmp)
          else
            state%velocity(1:3, i1, i2, i3) = 0.d0
            state%temperature(i1, i2, i3) = 0.d0
          endif
        enddo
      enddo
    enddo

    state%velocity(:,:,:,:) = (state%velocity(:,:,:,:)  &
    &  * initial_noise_level)
    state%temperature(:,:,:) = (state%temperature(:,:,:)  &
    &  * initial_noise_level)

    return
  end subroutine white_noise_low_modes

! *****************************************************************************

  subroutine interaction_test()
    use Parallel, only : my_rank
    use Geometry, only : state
    implicit none

    if (my_rank == 0) then
      state%velocity(1,2,2,0) = 0.1d0
    endif

    if (my_rank == 1) then
      state%velocity(2,2,2,0) = 0.2d0
    endif

    call preconditioning()

    return
  end subroutine interaction_test

! *****************************************************************************

  subroutine preconditioning()
    use mpi, only : MPI_COMM_WORLD
    use mpi, only : MPI_DOUBLE_PRECISION
    use mpi, only : MPI_SUM
    use Parallel, only : my_rank
    use FFTW_Interface, only : fftw_backward_parallel
    use FFTW_Interface, only : fftw_forward_parallel
    use Geometry, only : state
    use Geometry, only : dealias_truncate
    use Geometry, only : calculate_wavenumbers
    use Geometry, only : incompressible_projection
    use Geometry, only : n1, n2R_local, n3R_local
    use Time_Module, only: current_time
    use External_Field, only : derive_transform_tensor
    use Log, only : log_energy
    implicit none
    double precision :: A(1:3, 1:3)
    double precision :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: R(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    complex(kind(0d0)), pointer :: U(:,:,:,:)
    complex(kind(0d0)), pointer :: T(:,:,:)
    complex(kind(0d0)) :: PK(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    integer :: i

    U => state%velocity(:,:,:,:)
    T => state%temperature(:,:,:)

    do i = 1, 3
      call fftw_backward_parallel(U(i,:,:,:), R(:,:,:))
      call fftw_forward_parallel(R(:,:,:), U(i,:,:,:))
    enddo

    call fftw_backward_parallel(T(:,:,:), R(:,:,:))
    call fftw_forward_parallel(R(:,:,:), T(:,:,:))

    call derive_transform_tensor(current_time, A)
    call calculate_wavenumbers(A, K)
    call incompressible_projection(K, K, U, T, PK)

    call dealias_truncate(U, T)

    return
  end subroutine preconditioning

! *****************************************************************************

  subroutine eliminate_vortical_mode()
    use Geometry, only : state
    use Geometry, only : calculate_wavenumbers
    use External_Field, only : derive_transform_tensor
    use Analysis, only : wave_vortex_decomposition
    implicit none
    double precision :: time = 0.d0
    double precision :: A(1:3, 1:3)
    double precision :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), pointer :: U(:,:,:,:)
    complex(kind(0d0)), pointer :: T(:,:,:)
    complex(kind(0d0)) :: UW(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: TW(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: UV(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: TV(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    U => state%velocity(:,:,:,:)
    T => state%temperature(:,:,:)

    call derive_transform_tensor(time, A)
    call calculate_wavenumbers(A, K)
    call wave_vortex_decomposition(K, U, T, UW, TW, UV, TV)

    U(:,:,:,:) = UW(:,:,:,:)
    T(:,:,:) = TW(:,:,:)

    return
  end subroutine eliminate_vortical_mode

! *****************************************************************************

  subroutine adjust_to_power_law(power)
    use Parallel, only : my_rank
    use Geometry, only : state
    use Geometry, only : calculate_wavenumbers
    use Geometry, only : delta_K
    use External_Field, only : derive_transform_tensor
    implicit none
    double precision, intent(in) :: power
    double precision :: time = 0.d0
    double precision :: A(1:3, 1:3)
    complex(kind(0d0)), pointer :: U(:,:,:,:)
    complex(kind(0d0)), pointer :: T(:,:,:)
    double precision :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: K_abs(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: factor(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: K_min
    integer :: i

    U => state%velocity(:,:,:,:)
    T => state%temperature(:,:,:)

    call derive_transform_tensor(time, A)
    call calculate_wavenumbers(A, K)
    K_abs(:,:,:) = 0.d0
    do i = 1, 3
      K_abs(:,:,:) = K_abs(:,:,:) + K(i,:,:,:)**2 * 0.5d0
    enddo
    K_abs = sqrt(K_abs)

    K_min = minval(delta_K)

    factor(:,:,:) = (K_abs(:,:,:) / K_min)**(-1 + power/2)
    if (my_rank == 0) then
      factor(0, 0, 0) = 1.d0
    endif
    do i = 1, 3
      U(i,:,:,:) = U(i,:,:,:) * factor
    enddo
    T(:,:,:) = T(:,:,:) * factor

    return
  end subroutine adjust_to_power_law

end module Initial_Condition

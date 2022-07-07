module Analysis
  use mpi, only : MPI_COMM_WORLD
  use mpi, only : MPI_DOUBLE_PRECISION
  use mpi, only : MPI_SUM
  use Parallel, only : rank_1
  use Parallel, only : rank_2
  use Geometry, only : n1C_local, n2C_local, n3
  use Geometry, only : n1, n2R_local, n3R_local
  use Geometry, only : n2
  implicit none

  integer, allocatable :: index_HV(:,:,:,:)

  public initialize_index
  public calc_flux_convergence
  public derive_wave_vortex_energy
  public wave_vortex_decomposition
  public get_primary_energy
  public operator_Q
  public calculate_PV
  public sum_real_data
  public spectrum_integrate_azimuth
  public analysis_Richardson_distribution
  public Thorpe_displacements

contains

! *****************************************************************************

  subroutine initialize_index(K)
    use Geometry, only : NH_max
    use Geometry, only : NV_max
    use Geometry, only : delta_KH
    use Geometry, only : delta_KV
    implicit none
    double precision, intent(in) ::  &
      &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision :: KH
    double precision :: KV
    integer :: i1
    integer :: i2
    integer :: i3

    integer :: jH
    integer :: jV

    if (.not. allocated(index_HV)) then
      allocate(index_HV(1:2, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1), source=0)
    endif
    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local - 1
        do i1 = 0, n1C_local - 1

          KH = sqrt(K(1, i1, i2, i3)**2 + K(2, i1, i2, i3)**2)
          KV = abs(K(3, i1, i2, i3))

          jH = int(KH / delta_KH + 0.5d0)
          jV = int(KV / delta_KV + 0.5d0)

          index_HV(1, i1, i2, i3) = jH
          index_HV(2, i1, i2, i3) = jV

        enddo
      enddo
    enddo

    return
  end subroutine initialize_index

! *****************************************************************************

  subroutine calc_flux_convergence(U, T, K, HV_KFC, HV_PFC)
    use Parallel, only : my_rank
    use Geometry, only : NH_max
    use Geometry, only : NV_max
    use Geometry, only : delta_KH
    use Geometry, only : delta_KV
    implicit none
    complex(kind(0d0)), intent(in) ::  &
      &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(in) ::  &
      &  T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) ::  &
      &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(out) ::  &
      &  HV_KFC(0:NH_max, 0:NV_max)
    double precision, intent(out) ::  &
      &  HV_PFC(0:NH_max, 0:NV_max)
    double precision :: KFC_tmp(0:NH_max, 0:NV_max)
    double precision :: PFC_tmp(0:NH_max, 0:NV_max)
    double precision :: spec_global(0:NH_max, 0:NV_max)

    double precision :: KH
    double precision :: KV
    double precision :: KE
    double precision :: PE
    double precision :: factor

    integer :: i
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: i1_global

    integer :: jH
    integer :: jV
    integer :: jH_past
    integer :: jV_past

    integer :: ierr

    KFC_tmp(:,:) = 0.d0
    PFC_tmp(:,:) = 0.d0
    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local - 1
        do i1 = 0, n1C_local - 1
          i1_global = rank_1 * n1C_local + i1

          KH = sqrt(K(1, i1, i2, i3)**2 + K(2, i1, i2, i3)**2)
          KV = abs(K(3, i1, i2, i3))

          jH = int(KH / delta_KH + 0.5d0)
          jV = int(KV / delta_KV + 0.5d0)
          jH_past = index_HV(1, i1, i2, i3)
          jV_past = index_HV(2, i1, i2, i3)

          if ((jH <= NH_max .and. jV <= NV_max)  &
            &  .and. (jH_past <= NH_max .and. jV_past <= NV_max)  &
            &  .and. (jH /= jH_past .or. jV /= jV_past)) then

            if (i1_global == 0) then
              factor = 1.d0
            else
              factor = 2.d0
            endif

            KE = sum(abs(U(:, i1, i2, i3))**2) * 0.5d0 * factor
            KFC_tmp(jH_past, jV_past) = KFC_tmp(jH_past, jV_past) - KE
            KFC_tmp(jH, jV) = KFC_tmp(jH, jV) + KE

            PE = abs(T(i1, i2, i3))**2 * 0.5d0 * factor
            PFC_tmp(jH_past, jV_past) = PFC_tmp(jH_past, jV_past) - PE
            PFC_tmp(jH, jV) = PFC_tmp(jH, jV) + PE
          endif
          index_HV(1, i1, i2, i3) = jH
          index_HV(2, i1, i2, i3) = jV
        enddo
      enddo
    enddo

    spec_global(:,:) = 0.d0
    do i = 0, NV_max
      call MPI_Reduce(KFC_tmp(0,i), spec_global(0,i), NH_max+1,  &
        &  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    enddo
    HV_KFC(:,:) = HV_KFC(:,:) + spec_global(:,:)

    spec_global(:,:) = 0.d0
    do i = 0, NV_max
      call MPI_Reduce(PFC_tmp(0,i), spec_global(0,i), NH_max+1,  &
        &  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    enddo
    HV_PFC(:,:) = HV_PFC(:,:) + spec_global(:,:)

    return
  end subroutine calc_flux_convergence

! *****************************************************************************

  subroutine derive_wave_vortex_energy(U, T, time,  &
      &  wave_spectrum, vortex_spectrum)
    use Parameters, only : N => buoyancy_frequency
    use Geometry, only : calculate_wavenumbers
    use External_Field, only : derive_transform_tensor

    implicit none
    complex(kind(0d0)), intent(in) ::  &
      &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(in) :: T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) :: time
    double precision, intent(out) ::  &
      &  wave_spectrum(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(out) ::  &
      &  vortex_spectrum(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    complex(kind(0d0)) :: UW(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: TW(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: UV(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: TV(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision :: A(1:3, 1:3)
    double precision :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    integer :: i

    call derive_transform_tensor(time, A)
    call calculate_wavenumbers(A, K)
    call wave_vortex_decomposition(K, U, T, UW, TW, UV, TV)

    wave_spectrum(:,:,:) = 0.d0
    do i = 1, 3
      wave_spectrum(:,:,:) = wave_spectrum(:,:,:)  &
        &  + abs(UW(i,:,:,:))**2 * 0.5d0
    enddo
    wave_spectrum(:,:,:) = wave_spectrum(:,:,:) + abs(TW(:,:,:))**2 * 0.5d0

    vortex_spectrum(:,:,:) = 0.d0
    do i = 1, 3
      vortex_spectrum(:,:,:) = vortex_spectrum(:,:,:)  &
        &  + abs(UV(i,:,:,:))**2 * 0.5d0
    enddo
    vortex_spectrum(:,:,:) = vortex_spectrum(:,:,:) + abs(TV(:,:,:))**2 * 0.5d0

    return
  end subroutine derive_wave_vortex_energy

! *****************************************************************************

  subroutine wave_vortex_decomposition(K, U, T,  &
    &  UW, TW, UV, TV)
  use Parameters, only : N => buoyancy_frequency
  use Parameters, only : Omega => angular_velocity
  use External_Field, only : alpha_plus
  use External_Field, only : alpha_minus
  implicit none
  double precision, intent(in) ::  &
    &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
  complex(kind(0d0)), intent(in) ::  &
    &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
  complex(kind(0d0)), intent(in) ::  &
    &  T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
  complex(kind(0d0)), intent(out) ::  &
    &  UW(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
  complex(kind(0d0)), intent(out) ::  &
    &  TW(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
  complex(kind(0d0)), intent(out) ::  &
    &  UV(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
  complex(kind(0d0)), intent(out) ::  &
    &  TV(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

  complex(kind(0d0)), parameter :: ui = (0.d0, 1.d0)
  complex(kind(0d0)) :: Q(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
  double precision :: C(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
  double precision :: f

  double precision :: epsilon = 1.d-20

  integer :: i1
  integer :: i2
  integer :: i3

  f = 2 * Omega(3)
  C = N**2 * K(1,:,:,:)**2 + N**2 * K(2,:,:,:)**2  &
    &  + (f - alpha_plus - alpha_minus)**2 * K(3,:,:,:)**2

  do i1 = 0, n1C_local - 1
    do i2 = 0, n2C_local - 1
      do i3 = 0, n3 - 1
        if (abs(C(i1, i2, i3)) < epsilon) then
          C(i1, i2, i3) = 1.d0
        else
          C(i1, i2, i3) = 1.d0 / C(i1, i2, i3)
        endif
      enddo
    enddo
  enddo

  call operator_Q(K, U, T, Q)

  UV(1,:,:,:) = C * ui * N * K(2,:,:,:) * Q(:,:,:)
  UV(2,:,:,:) = - C * ui * N * K(1,:,:,:) * Q(:,:,:)
  UV(3,:,:,:) = 0.d0
  TV(:,:,:) = - C * ui * (f - alpha_plus - alpha_minus)  &
    &  * K(3,:,:,:) * Q(:,:,:)

  UW(:,:,:,:) = U(:,:,:,:) - UV(:,:,:,:)
  TW(:,:,:) = T(:,:,:) - TV(:,:,:)

  return
  end subroutine wave_vortex_decomposition

! *****************************************************************************

  function get_primary_energy(wave_spectrum) result(primary_energy)
    implicit none
    double precision, intent(in) ::  &
      &  wave_spectrum(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: primary_energy_local
    double precision :: primary_energy

    integer :: i1
    integer :: i2
    integer :: i3
    integer :: i1g
    integer :: i2g
    integer :: ierr

    logical :: flag

    primary_energy_local = 0.d0
    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local - 1
        i2g = rank_2 * n2C_local + i2
        do i1 = 0, n1C_local - 1
          i1g = rank_1 * n1C_local + i1
          flag = (((i1g == 1) .and. (i2g == 0) .and. (i3 == 1))  &
            &  .or. ((i1g == 1) .and. (i2g == 0) .and. (i3 == n3 - 1)))
          if (flag) then
            primary_energy_local = primary_energy_local  &
              &  + wave_spectrum(i1, i2, i3) * 2
          endif
          flag = (((i1g == 0) .and. (i2g == 1) .and. (i3 == 1))  &
            &  .or. ((i1g == 0) .and. (i2g == 1) .and. (i3 == n3 - 1))  &
            &  .or. ((i1g == 0) .and. (i2g == n2 - 1) .and. (i3 == 1))  &
            &  .or. ((i1g == 0) .and. (i2g == n2 - 1) .and. (i3 == n3 - 1)))
          if (flag) then
            primary_energy_local = primary_energy_local  &
              &  + wave_spectrum(i1, i2, i3)
          endif
        enddo
      enddo
    enddo

    primary_energy = 0.d0
    call MPI_Reduce(primary_energy_local, primary_energy, 1,  &
      &  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    return
  end function get_primary_energy

! *****************************************************************************

  subroutine operator_Q(K, U, T, Q)
    use Parameters, only : N => buoyancy_frequency
    use Parameters, only : Omega => angular_velocity
    use External_Field, only : alpha_plus
    use External_Field, only : alpha_minus
    implicit none
    double precision, intent(in) ::  &
      &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(in) ::  &
      &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(in) ::  &
      &  T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(out) ::  &
      &  Q(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    complex(kind(0d0)), parameter :: ui = (0.d0, 1.d0)
    double precision :: f

    f = 2 * Omega(3)

    Q(:,:,:) = - ui * N * K(2,:,:,:) * U(1,:,:,:)  &
      &  + ui * N * K(1,:,:,:) * U(2,:,:,:)  &
      &  + ui * (f - alpha_plus - alpha_minus) * K(3,:,:,:) * T(:,:,:)

    return
  end subroutine operator_Q

! *****************************************************************************

  subroutine calculate_PV(U, T, time, Q1, Q2)
    use FFTW_Interface, only : fftw_forward_parallel
    use FFTW_Interface, only : fftw_backward_parallel
    use Geometry, only : calculate_wavenumbers
    use Geometry, only : wavenumber_moving_frame_local
    use Geometry, only : K1_truncated
    use Geometry, only : K2_truncated
    use Geometry, only : K3_truncated
    use Geometry, only : get_mask
    use External_Field, only : derive_transform_tensor
    implicit none
    complex(kind(0d0)), intent(in) ::  &
      &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(in) ::  &
      &  T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) :: time
    complex(kind(0d0)), intent(out) ::  &
      &  Q1(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(out) ::  &
      &  Q2(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    complex(kind(0d0)), parameter :: ui = (0.d0, 1.d0)
    double precision :: A(1:3, 1:3)
    double precision :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    complex(kind(0d0)) :: vor1(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: vor2(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: vor3(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: dTdx1(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: dTdx2(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: dTdx3(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: vor1R(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: vor2R(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: vor3R(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: dTdx1R(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: dTdx2R(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: dTdx3R(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: Q2R(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)

    integer :: mask(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    call derive_transform_tensor(time, A)
    call calculate_wavenumbers(A, K)
    call operator_Q(K, U, T, Q1)

    vor1(:,:,:) = ui * (K(2,:,:,:) * U(3,:,:,:) - K(3,:,:,:) * U(2,:,:,:))
    vor2(:,:,:) = ui * (K(3,:,:,:) * U(1,:,:,:) - K(1,:,:,:) * U(3,:,:,:))
    vor3(:,:,:) = ui * (K(1,:,:,:) * U(2,:,:,:) - K(2,:,:,:) * U(1,:,:,:))
    call fftw_backward_parallel(vor1, vor1R)
    call fftw_backward_parallel(vor2, vor2R)
    call fftw_backward_parallel(vor3, vor3R)

    dTdx1(:,:,:) = ui * K(1,:,:,:) * T(:,:,:)
    dTdx2(:,:,:) = ui * K(2,:,:,:) * T(:,:,:)
    dTdx3(:,:,:) = ui * K(3,:,:,:) * T(:,:,:)
    call fftw_backward_parallel(dTdx1, dTdx1R)
    call fftw_backward_parallel(dTdx2, dTdx2R)
    call fftw_backward_parallel(dTdx3, dTdx3R)

    Q2R(:,:,:) = vor1R(:,:,:) * dTdx1R(:,:,:)  &
      &  + vor2R(:,:,:) * dTdx2R(:,:,:)  &
      &  + vor3R(:,:,:) * dTdx3R(:,:,:)
    call fftw_forward_parallel(Q2R, Q2)
    call get_mask(mask)
    Q2(:,:,:) = Q2(:,:,:) * mask(:,:,:)

    return
  end subroutine calculate_PV

! *****************************************************************************

  function sum_real_data(data_) result(global_sum)
    implicit none
    double precision, intent(in) ::  &
      &  data_(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: global_sum
    double precision :: local_sum
    integer :: i1
    integer :: i1_global

    integer :: ierr

    local_sum = 0.d0
    do i1 = 0, n1C_local - 1
      i1_global = rank_1 * n1C_local + i1
      if (i1_global == 0) then
        local_sum = local_sum + sum(data_(i1,:,:))
      else
        local_sum = local_sum + sum(data_(i1,:,:)) * 2
      endif
    enddo

    global_sum = 0.d0
    call MPI_Reduce(local_sum, global_sum, 1, MPI_DOUBLE_PRECISION,  &
      &  MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    return
  end function sum_real_data

! *****************************************************************************

  subroutine spectrum_integrate_azimuth(data_, K, spectrum)
    use Geometry, only : NH_max
    use Geometry, only : NV_max
    use Geometry, only : delta_KH
    use Geometry, only : delta_KV
    implicit none
    double precision, intent(in) ::  &
      &  data_(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) ::  &
      &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(out) :: spectrum(0:NH_max, 0:NV_max)

    double precision, allocatable :: spectrum_local(:,:)
    double precision, allocatable :: spectrum_global(:,:)
    double precision :: spectrum_sum_local
    double precision :: spectrum_sum_global

    double precision :: factor
    double precision :: KH
    double precision :: KV

    integer :: i1
    integer :: i2
    integer :: i3
    integer :: i1_global
    integer :: i
    integer :: jH
    integer :: jV

    integer :: N_data_size
    integer :: ierr

    allocate(  &
      &  spectrum_local(0:NH_max, 0:NV_max),  &
      &  spectrum_global(0:NH_max, 0:NV_max),  &
      &  source=0.d0)

    spectrum(:,:) = 0.d0

    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local-1
        do i1 = 0, n1C_local-1
          i1_global = rank_1 * n1C_local + i1
          if (i1_global == 0) then
            factor = 1.d0
          else
            factor = 2.d0
          endif

          KH = sqrt(K(1, i1, i2, i3)**2 + K(2, i1, i2, i3)**2)
          KV = abs(K(3, i1, i2, i3))

          jH = int(KH / delta_KH + 0.5d0)
          jV = int(KV / delta_KV + 0.5d0)

          if (jH <= NH_max .and. jV <= NV_max) then
            spectrum_local(jH, jV) = spectrum_local(jH, jV)  &
              &  + data_(i1, i2, i3) * factor
          endif

        enddo
      enddo
    enddo

    N_data_size = (NH_max + 1) * (NV_max + 1)

    do i = 0, NV_max
      call MPI_Reduce(spectrum_local(0,i), spectrum_global(0,i), NH_max+1,  &
        &  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    enddo

    spectrum_sum_local = sum(spectrum_local)
    call MPI_Reduce(spectrum_sum_local, spectrum_sum_global, 1,  &
      &  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    spectrum(:,:) = spectrum_global(:,:)

    return
  end subroutine spectrum_integrate_azimuth

! *****************************************************************************

  subroutine spectrum_project_horizontal(data_, K, spectrum)
    use Geometry, only : calculate_wavenumbers
    use Geometry, only : n2_truncated, n1_truncated
    use Geometry, only : NH_max, NH_min
    use Geometry, only : delta_KH
    use External_Field, only : derive_transform_tensor
    implicit none
    double precision, intent(in) ::  &
      &  data_(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) ::  &
      &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(out) :: spectrum(0:NH_max, 0:NH_min*2)

    double precision, allocatable :: spectrum_local(:,:)
    double precision, allocatable :: spectrum_global(:,:)

    double precision :: factor

    integer :: i
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: i1_global

    integer :: j1
    integer :: j2
    integer :: j1_s
    integer :: j2_s

    integer :: ierr

    allocate(  &
      &  spectrum_local(0:NH_max, 0:NH_min*2),  &
      &  spectrum_global(0:NH_max, 0:NH_min*2),  &
      &  source=0.d0)

    spectrum(:,:) = 0.d0

    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local - 1
        do i1 = 0, n1C_local - 1

          i1_global = rank_1 * n1C_local + i1
          if (i1_global == 0) then
            factor = 1.d0
          else
            factor = 2.d0
          endif

          j1 = int(K(1, i1, i2, i3) / delta_KH + 0.5d0)
          j2 = int(K(2, i1, i2, i3) / delta_KH + 0.5d0)

          if (abs(j1) <= NH_min .and. abs(j2) <= NH_max) then
            if (j2 >= 0) then
              j1_s = j1 + NH_min
              j2_s = j2
            else
              j1_s = - j1 + NH_min
              j2_s = - j2
            endif
            spectrum_local(j2_s, j1_s) = spectrum_local(j2_s, j1_s)  &
              &  + data_(i1, i2, i3) * factor
          endif

        enddo
      enddo
    enddo

    do i = 0, NH_min*2
      call MPI_Reduce(spectrum_local(0,i), spectrum_global(0,i), NH_max+1,  &
        &  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    enddo

    spectrum(:,:) = spectrum_global(:,:)

    return
  end subroutine spectrum_project_horizontal

! *****************************************************************************

  subroutine analysis_Richardson_distribution(U, T, time,  &
    &  Richardson_distribution)
    use Geometry, only : Ri_max
    use Geometry, only : Ri_min
    use Geometry, only : N_Ri
    implicit none
    complex(kind(0d0)), intent(in) ::  &
      &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(in) :: T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) :: time
    double precision, intent(out) :: Richardson_distribution(0:N_Ri-1)

    double precision :: Richardson_number(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: distribution_local(0:N_Ri-1)
    double precision :: Del_Ri
    double precision :: tmp
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: j

    integer :: ierr

    Del_Ri = (Ri_max - Ri_min) / (N_Ri - 2)

    call analysis_Richardson_number(U, T, time, Richardson_number)

    distribution_local(:) = 0.0d0

    do i1 = 0, n1 - 1
      do i2 = 0, n2R_local - 1
        do i3 = 0, n3R_local - 1
          tmp = (Richardson_number(i1, i2, i3) - Ri_min) / Del_Ri
          if (tmp < 0) then
            distribution_local(0) = distribution_local(0) + 1.0d0
          elseif (tmp >= N_Ri) then
            distribution_local(N_Ri-1) = distribution_local(N_Ri-1) + 1.0d0
          else
            j = int(tmp)
            distribution_local(j) = distribution_local(j) + 1.0d0
          end if
        end do
      end do
    end do

    call MPI_Allreduce(distribution_local(0), Richardson_distribution(0),  &
      &  N_Ri, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    Richardson_distribution(:) = Richardson_distribution(:)  &
      &  / (dble(n1) * n2 * n3 * Del_Ri)

    return
  end subroutine analysis_Richardson_distribution

! *****************************************************************************

  subroutine analysis_Richardson_number(U, T, time, Richardson_number)
    use Geometry, only : domain_length
    use Geometry, only : calculate_wavenumbers
    use FFTW_Interface, only : fftw_backward_parallel
    use Parameters, only : buoyancy_frequency
    use External_Field, only : derive_transform_tensor
    use Parallel, only : my_rank
    implicit none
    complex(kind(0d0)), intent(in) ::  &
      &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(in) :: T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) :: time
    double precision, intent(out) :: Richardson_number(  &
      &  0:n1-1, 0:n2R_local-1, 0:n3R_local-1)

    complex(kind(0d0)), parameter :: ui = (0.d0, 1.d0)
    double precision :: A(1:3, 1:3)
    double precision :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    complex(kind(0d0)) :: work(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: dUdZ(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: dVdZ(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: dTdZ(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)

    call derive_transform_tensor(time, A)
    call calculate_wavenumbers(A, K)

    work(:,:,:) = ui * K(3,:,:,:) * U(1,:,:,:)
    call fftw_backward_parallel(work, dUdZ)

    work(:,:,:) = ui * K(3,:,:,:) * U(2,:,:,:)
    call fftw_backward_parallel(work, dVdZ)

    work(:,:,:) = ui * K(3,:,:,:) * T(:,:,:)
    call fftw_backward_parallel(work, dTdZ)

    Richardson_number(:,:,:) = buoyancy_frequency * (  &
      &  buoyancy_frequency + dTdZ(:,:,:)) / (dUdZ(:,:,:)**2 + dVdZ(:,:,:)**2)

    return
  end subroutine analysis_Richardson_number

! *****************************************************************************

  subroutine Thorpe_displacements(T, PDF)
    use Parallel, only : comm_2
    use Parallel, only : rank_1
    use Parallel, only : rank_2
    use Geometry, only : redistribute_real_data
    use Geometry, only : n1R_local
    use Geometry, only : division
    use Geometry, only : domain_length
    use Geometry, only : N_Thorpe
    use Parameters, only : N => buoyancy_frequency
    use FFTW_Interface, only : fftw_backward_parallel
    implicit none
    complex(kind(0d0)), intent(in) :: T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(out) :: PDF(0:N_Thorpe-1)

    double precision :: TR(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: T_redist(0:n1R_local-1, 0:n2R_local-1, 0:n3-1)
    double precision :: profile(0:n3-1)
    double precision :: displacement
    double precision :: PDF_local(0:N_Thorpe-1)

    double precision :: delta_z
    double precision :: vertical_axis(0:n3-1)

    double precision :: T_tmp1
    double precision :: T_tmp2
    double precision :: displacement_min
    double precision :: displacement_max

    integer :: index_stable(0:n3-1)

    integer :: i1
    integer :: i2
    integer :: i3
    integer :: i3s
    integer :: i_tmp
    integer :: j
    integer :: ierr

    call fftw_backward_parallel(T, TR)
    call redistribute_real_data(TR, T_redist)

    delta_z = domain_length(3) / n3
    displacement_max = domain_length(3) - delta_z / 2
    displacement_min = - domain_length(3) + delta_z / 2

    PDF_local = 0.d0
    do i3 = 0, n3-1
      vertical_axis(i3) = i3 * delta_z
    enddo

    do i2 = 0, n2R_local - 1
      do i1 = 0, n1R_local - 1

        do i3 = 0, n3 - 1
          profile(i3) = T_redist(i1, i2, i3) + N * vertical_axis(i3)
          index_stable(i3) = i3
        enddo

        do i3 = 0, n3 - 2
          T_tmp1 = profile(i3)
          do i3s = i3 + 1, n3 - 1
            T_tmp2 = profile(i3s)
            if (T_tmp1 > T_tmp2) then
              profile(i3) = T_tmp2
              profile(i3s) = T_tmp1
              T_tmp1 = T_tmp2
              i_tmp = index_stable(i3)
              index_stable(i3) = index_stable(i3s)
              index_stable(i3s) = i_tmp
            endif
          enddo
        enddo

        do i3 = 0, n3 - 1
          displacement = vertical_axis(index_stable(i3))  &
            &  - vertical_axis(i3)
          j = int((displacement - displacement_min) / delta_z)
          PDF_local(j) = PDF_local(j) + 1.d0
        enddo
      enddo
    enddo

    call MPI_Allreduce(PDF_local(0), PDF(0), N_Thorpe,  &
      &  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    PDF(:) = PDF(:) / sum(PDF)

    return
  end subroutine Thorpe_displacements

end module Analysis

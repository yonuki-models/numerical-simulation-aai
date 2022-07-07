module Dynamics
  use Geometry, only : n1C_local, n2C_local, n3
  use Geometry, only : n1, n2R_local, n3R_local
  use Geometry, only : n2
  use Geometry, only : NH_max, NH_min, NV_max
  use Analysis, only : spectrum_project_horizontal
  use Analysis, only : spectrum_integrate_azimuth
  use Geometry, only : state
  use Geometry, only : energy
  use Geometry, only : spectrum_H
  use Geometry, only : spectrum_HV
  implicit none
  private
  double precision, allocatable, public :: decay_rate(:,:,:)
  double precision, public :: decay_coefficient_H
  double precision, public :: decay_coefficient_V

  public linear_terms
  public nonlinear_terms
  public wave_vortex_production
  public derive_viscous_factor
  public viscous_parameterization

contains

! *****************************************************************************

  subroutine linear_terms(K, D, M, delta_time, factor_energy,  &
    &  term_U, term_T)
    use Parameters, only : N => buoyancy_frequency
    use Parameters, only : Omega => angular_velocity
    use Geometry, only : NH_max
    use Geometry, only : NV_max
    use Geometry, only : delta_KH
    use Geometry, only : delta_KV
    implicit none
    double precision, intent(in) ::  &
      &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) :: D(1:3, 1:3)
    double precision, intent(in) :: M(1:3)
    double precision, intent(in) :: delta_time
    double precision, intent(in) :: factor_energy
    complex(kind(0d0)), intent(out) ::  &
      &  term_U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(out) ::  &
      &  term_T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    complex(kind(0d0)), pointer :: U(:,:,:,:)
    complex(kind(0d0)), pointer :: T(:,:,:)
    double precision, pointer :: PKE(:,:,:)
    double precision, pointer :: PPE(:,:,:)
    double precision, pointer :: CE(:,:,:)
    double precision, pointer :: DKE(:,:,:)
    double precision, pointer :: DPE(:,:,:)

    double precision, pointer :: HS_PKE(:,:)
    double precision, pointer :: HS_PPE(:,:)
    double precision, pointer :: HS_CE(:,:)
    double precision, pointer :: HS_DKE(:,:)
    double precision, pointer :: HS_DPE(:,:)

    double precision, pointer :: HV_PKE(:,:)
    double precision, pointer :: HV_PPE(:,:)
    double precision, pointer :: HV_CE(:,:)
    double precision, pointer :: HV_DKE(:,:)
    double precision, pointer :: HV_DPE(:,:)

    double precision :: K_squared(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision :: energy_tmp(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: HS_tmp(0:NH_max, 0:NH_min*2)
    double precision :: HV_tmp(0:NH_max, 0:NV_max)

    ! Levi-Civita symbol
    integer :: j_plus(1:3) = 0
    integer :: k_plus(1:3) = 0
    integer :: j_minus(1:3) = 0
    integer :: k_minus(1:3) = 0

    integer :: i1
    integer :: i2
    integer :: i3
    integer :: i
    integer :: j
    integer :: jH
    integer :: jV

    U => state%velocity(:,:,:,:)
    T => state%temperature(:,:,:)

    PKE => energy%kinetic_energy_production(:,:,:)
    PPE => energy%potential_energy_production(:,:,:)
    CE => energy%energy_conversion(:,:,:)
    DKE => energy%kinetic_energy_dissipation(:,:,:)
    DPE => energy%potential_energy_dissipation(:,:,:)

    HS_PKE => spectrum_H%kinetic_energy_production(:,:)
    HS_PPE => spectrum_H%potential_energy_production(:,:)
    HS_CE => spectrum_H%energy_conversion(:,:)
    HS_DKE => spectrum_H%kinetic_energy_dissipation(:,:)
    HS_DPE => spectrum_H%potential_energy_dissipation(:,:)

    HV_PKE => spectrum_HV%kinetic_energy_production(:,:)
    HV_PPE => spectrum_HV%potential_energy_production(:,:)
    HV_CE => spectrum_HV%energy_conversion(:,:)
    HV_DKE => spectrum_HV%kinetic_energy_dissipation(:,:)
    HV_DPE => spectrum_HV%potential_energy_dissipation(:,:)

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

!! Energy variables
    energy_tmp(:,:,:) = 0.d0
    do i = 1, 3
      do j = 1, 3
        energy_tmp(:,:,:) = energy_tmp(:,:,:)  &
        &  + factor_energy * delta_time * (-D(i,j))  &
        &  * dble(conjg(U(i,:,:,:)) * U(j,:,:,:))
      enddo
    enddo
    PKE(:,:,:) = PKE(:,:,:) + energy_tmp(:,:,:)
    call spectrum_project_horizontal(energy_tmp, K, HS_tmp)
    HS_PKE(:,:) = HS_PKE(:,:) + HS_tmp(:,:)
    call spectrum_integrate_azimuth(energy_tmp, K, HV_tmp)
    HV_PKE(:,:) = HV_PKE(:,:) + HV_tmp(:,:)

    energy_tmp(:,:,:) = 0.d0
    do i = 1, 3
      energy_tmp(:,:,:) = energy_tmp(:,:,:)  &
      &  + factor_energy * delta_time * (-M(i))  &
      &  * dble(conjg(T(:,:,:)) * U(i,:,:,:))
    enddo
    PPE(:,:,:) = PPE(:,:,:) + energy_tmp(:,:,:)
    call spectrum_project_horizontal(energy_tmp, K, HS_tmp)
    HS_PPE(:,:) = HS_PPE(:,:) + HS_tmp(:,:)
    call spectrum_integrate_azimuth(energy_tmp, K, HV_tmp)
    HV_PPE(:,:) = HV_PPE(:,:) + HV_tmp(:,:)

    energy_tmp(:,:,:) = factor_energy * delta_time * (-N)  &
    &  * dble(conjg(T(:,:,:)) * U(3,:,:,:))
    CE(:,:,:) = CE(:,:,:) + energy_tmp(:,:,:)
    call spectrum_project_horizontal(energy_tmp, K, HS_tmp)
    HS_CE(:,:) = HS_CE(:,:) + HS_tmp(:,:)
    call spectrum_integrate_azimuth(energy_tmp, K, HV_tmp)
    HV_CE(:,:) = HV_CE(:,:) + HV_tmp(:,:)

    K_squared(:,:,:) = 0.d0
    do i = 1, 3
      K_squared(:,:,:) = K_squared(:,:,:) + K(i,:,:,:)**2
    enddo

    energy_tmp(:,:,:) = 0.d0
    do i = 1, 3
      energy_tmp(:,:,:) = energy_tmp(:,:,:) + factor_energy * delta_time  &
      &  * decay_rate(:,:,:) * dble(conjg(U(i,:,:,:)) * U(i,:,:,:))
    enddo
    DKE(:,:,:) = DKE(:,:,:) + energy_tmp(:,:,:)
    call spectrum_project_horizontal(energy_tmp, K, HS_tmp)
    HS_DKE(:,:) = HS_DKE(:,:) + HS_tmp(:,:)
    call spectrum_integrate_azimuth(energy_tmp, K, HV_tmp)
    HV_DKE(:,:) = HV_DKE(:,:) + HV_tmp(:,:)

    energy_tmp(:,:,:) = factor_energy * delta_time  &
    &  * decay_rate(:,:,:) * dble(conjg(T(:,:,:)) * T(:,:,:))
    DPE(:,:,:) = DPE(:,:,:) + energy_tmp(:,:,:)
    call spectrum_project_horizontal(energy_tmp, K, HS_tmp)
    HS_DPE(:,:) = HS_DPE(:,:) + HS_tmp(:,:)
    call spectrum_integrate_azimuth(energy_tmp, K, HV_tmp)
    HV_DPE(:,:) = HV_DPE(:,:) + HV_tmp(:,:)

!! Governing equations
    term_U(:,:,:,:) = (0.d0, 0.d0)
    do i = 1, 3
      do j = 1, 3
        term_U(i,:,:,:) = term_U(i,:,:,:) - D(i,j) * U(j,:,:,:)
      enddo
    enddo

    do i = 1, 3
      term_U(i,:,:,:) = term_U(i,:,:,:) - 2 * (  &
        &  Omega(j_plus(i)) * U(k_plus(i),:,:,:)  &
        &  - Omega(j_minus(i)) * U(k_minus(i),:,:,:) )
    enddo
    term_U(3,:,:,:) = term_U(3,:,:,:) + N * T(:,:,:)

    term_T(:,:,:) = (0.d0, 0.d0)
    do i = 1, 3
      term_T(:,:,:) = term_T(:,:,:) - M(i) * U(i,:,:,:)
    enddo
    term_T(:,:,:) = term_T(:,:,:) - N * U(3,:,:,:)

    term_U(:,:,:,:) = term_U(:,:,:,:) * delta_time
    term_T(:,:,:) = term_T(:,:,:) * delta_time

    return
  end subroutine linear_terms

! *****************************************************************************

  subroutine nonlinear_terms(K, delta_time, factor_energy,  &
    &  term_U, term_T)
    use FFTW_Interface, only : fftw_forward_parallel
    use FFTW_Interface, only : fftw_backward_parallel
    implicit none
    double precision, intent(in) ::  &
      &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) :: delta_time
    double precision, intent(in) :: factor_energy
    complex(kind(0d0)), intent(inout) ::  &
      &  term_U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(inout) ::  &
      &  term_T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    complex(kind(0d0)), pointer :: U(:,:,:,:)
    complex(kind(0d0)), pointer :: T(:,:,:)
    double precision, pointer :: TKE(:,:,:)
    double precision, pointer :: TPE(:,:,:)
    double precision, pointer :: HS_TKE(:,:)
    double precision, pointer :: HS_TPE(:,:)
    double precision, pointer :: HV_TKE(:,:)
    double precision, pointer :: HV_TPE(:,:)

    double precision :: energy_tmp(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: HS_tmp(0:NH_max, 0:NH_min*2)
    double precision :: HV_tmp(0:NH_max, 0:NV_max)

    complex(kind(0d0)), parameter :: UI = (0.d0, 1.d0)
    double precision :: UR(1:3, 0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: TR(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision :: work_r(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    complex(kind(0d0)) :: work_a(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: work_b(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    integer :: i
    integer :: j

    U => state%velocity(:,:,:,:)
    T => state%temperature(:,:,:)
    TKE => energy%kinetic_energy_transfer(:,:,:)
    TPE => energy%potential_energy_transfer(:,:,:)
    HS_TKE => spectrum_H%kinetic_energy_transfer(:,:)
    HS_TPE => spectrum_H%potential_energy_transfer(:,:)
    HV_TKE => spectrum_HV%kinetic_energy_transfer(:,:)
    HV_TPE => spectrum_HV%potential_energy_transfer(:,:)

    do i = 1, 3
      call fftw_backward_parallel(U(i,:,:,:), UR(i,:,:,:))
    enddo
    call fftw_backward_parallel(T, TR)

    energy_tmp(:,:,:) = 0.d0
    do i = 1, 3
      work_r(:,:,:) = UR(i,:,:,:) * UR(i,:,:,:)
      call fftw_forward_parallel(work_r, work_a)
      work_b(:,:,:) = - UI * K(i,:,:,:) * work_a(:,:,:) * delta_time
      term_U(i,:,:,:) = term_U(i,:,:,:) + work_b(:,:,:)
      energy_tmp(:,:,:) = energy_tmp(:,:,:) + factor_energy  &
        &  * dble(conjg(U(i,:,:,:)) * work_b(:,:,:))
    enddo

    do i = 1, 3
      do j = 1, 3
        if (j >= i) cycle
        work_r(:,:,:) = UR(i,:,:,:) * UR(j,:,:,:)
        call fftw_forward_parallel(work_r, work_a)

        work_b(:,:,:) = - UI * K(j,:,:,:) * work_a(:,:,:) * delta_time
        term_U(i,:,:,:) = term_U(i,:,:,:) + work_b(:,:,:)
        energy_tmp(:,:,:) = energy_tmp(:,:,:) + factor_energy  &
          &  * dble(conjg(U(i,:,:,:)) * work_b(:,:,:))

        work_b(:,:,:) = - UI * K(i,:,:,:) * work_a(:,:,:) * delta_time
        term_U(j,:,:,:) = term_U(j,:,:,:) + work_b(:,:,:)
        energy_tmp(:,:,:) = energy_tmp(:,:,:) + factor_energy  &
          &  * dble(conjg(U(j,:,:,:)) * work_b(:,:,:))
      enddo
    enddo
    TKE(:,:,:) = TKE(:,:,:) + energy_tmp(:,:,:)
    call spectrum_project_horizontal(energy_tmp, K, HS_tmp)
    HS_TKE(:,:) = HS_TKE(:,:) + HS_tmp(:,:)
    call spectrum_integrate_azimuth(energy_tmp, K, HV_tmp)
    HV_TKE(:,:) = HV_TKE(:,:) + HV_tmp(:,:)

    energy_tmp(:,:,:) = 0.d0
    do i = 1, 3
      work_r(:,:,:) = UR(i,:,:,:) * TR(:,:,:)
      call fftw_forward_parallel(work_r, work_a)
      work_b(:,:,:) = - UI * K(i,:,:,:) * work_a(:,:,:) * delta_time
      term_T(:,:,:) = term_T(:,:,:) + work_b(:,:,:)
      energy_tmp(:,:,:) = energy_tmp(:,:,:) + factor_energy  &
        &  * dble(conjg(T(:,:,:)) * work_b(:,:,:))
    enddo
    TPE(:,:,:) = TPE(:,:,:) + energy_tmp(:,:,:)
    call spectrum_project_horizontal(energy_tmp, K, HS_tmp)
    HS_TPE(:,:) = HS_TPE(:,:) + HS_tmp(:,:)
    call spectrum_integrate_azimuth(energy_tmp, K, HV_tmp)
    HV_TPE(:,:) = HV_TPE(:,:) + HV_tmp(:,:)

    return
  end subroutine nonlinear_terms

! *****************************************************************************

  subroutine wave_vortex_production(K, D, delta_time, factor_energy)
    use Parallel, only : my_rank
    use Parameters, only : N => buoyancy_frequency
    use Parameters, only : Omega => angular_velocity
    use External_Field, only : alpha_plus
    use External_Field, only : alpha_minus
    use Analysis, only : wave_vortex_decomposition
    use Analysis, only : operator_Q
    implicit none
    double precision, intent(in) :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) :: D(1:3, 1:3)
    double precision, intent(in) :: delta_time
    double precision, intent(in) :: factor_energy

    complex(kind(0d0)), parameter :: UI = (0.d0, 1.d0)
    complex(kind(0d0)), pointer :: U(:,:,:,:)
    complex(kind(0d0)), pointer :: T(:,:,:)

    complex(kind(0d0)) :: Q(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: UW(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: TW(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: UV(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: TV(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision :: gamma_sq(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: gamma_U(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: gamma_V(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: gamma_T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision :: dot_gamma_sq(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: dot_gamma_U(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: dot_gamma_V(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: dot_gamma_T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision :: energy_tmp(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision, pointer :: HV_WP1(:,:)
    double precision, pointer :: HV_WP2(:,:)
    double precision, pointer :: HV_VP(:,:)
    double precision :: HV_tmp(0:NH_max, 0:NV_max)

    double precision :: absolute_vorticity
    double precision :: epsilon = 1.d-20

    integer :: i1
    integer :: i2
    integer :: i3

    integer :: i
    integer :: j

    U => state%velocity(:,:,:,:)
    T => state%temperature(:,:,:)
    HV_WP1 => spectrum_HV%wave_production_1(:,:)
    HV_WP2 => spectrum_HV%wave_production_2(:,:)
    HV_VP => spectrum_HV%vortex_production(:,:)

    absolute_vorticity = 2 * Omega(3) - alpha_minus - alpha_plus

    gamma_sq(:,:,:) = N**2 * K(1,:,:,:)**2 + N**2 * K(2,:,:,:)**2  &
      &  + absolute_vorticity**2 * K(3,:,:,:)**2
    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local - 1
        do i1 = 0, n1C_local - 1
          if (abs(gamma_sq(i1, i2, i3)) < epsilon) then
            gamma_sq(i1, i2, i3) = 1.d0
          endif
        enddo
      enddo
    enddo

    dot_gamma_sq = 2 * N**2 * K(1,:,:,:) * (- D(2,1) * K(2,:,:,:))  &
      &  + 2 * N**2 * K(2,:,:,:) * (- D(1,2) * K(1,:,:,:))
    gamma_U(:,:,:) = UI * N * K(2,:,:,:) / gamma_sq
    gamma_V(:,:,:) = - UI * N * K(1,:,:,:) / gamma_sq
    gamma_T(:,:,:) = - UI * absolute_vorticity * K(3,:,:,:) / gamma_sq

    dot_gamma_U(:,:,:) = (- dot_gamma_sq / gamma_sq * gamma_U  &
      &  + UI * N * (- alpha_plus * K(1,:,:,:)) / gamma_sq)
    dot_gamma_V(:,:,:) = (- dot_gamma_sq / gamma_sq * gamma_V  &
      &  - UI * N * (alpha_minus * K(2,:,:,:)) / gamma_sq)
    dot_gamma_T(:,:,:) = - dot_gamma_sq / gamma_sq * gamma_T

    call wave_vortex_decomposition(K, U, T, UW, TW, UV, TV)
    call operator_Q(K, U, T, Q)

    energy_tmp(:,:,:) = 0.d0
    do i = 1, 3
      do j = 1, 3
        energy_tmp(:,:,:) = energy_tmp(:,:,:)  &
          &  + factor_energy * delta_time * (-D(i,j))  &
          &  * dble(conjg(UW(i,:,:,:)) * UW(j,:,:,:))
      enddo
    enddo
    call spectrum_integrate_azimuth(energy_tmp, K, HV_tmp)
    HV_WP1(:,:) = HV_WP1(:,:) + HV_tmp(:,:)

    energy_tmp(:,:,:) = factor_energy * delta_time *  &
      &  (dble(conjg(UW(1,:,:,:))  &
      &  * (alpha_minus * gamma_V - dot_gamma_U) * Q)  &
      &  + dble(conjg(UW(2,:,:,:))  &
      &  * (- alpha_plus * gamma_U - dot_gamma_V) * Q)  &
      &  + dble(conjg(TW) * (- dot_gamma_T) * Q))
    call spectrum_integrate_azimuth(energy_tmp, K, HV_tmp)
    HV_WP2(:,:) = HV_WP2(:,:) + HV_tmp(:,:)

    energy_tmp(:,:,:) = factor_energy * delta_time *  &
      &  (- dot_gamma_sq / gamma_sq**2) * abs(Q)**2 * 0.5d0
    call spectrum_integrate_azimuth(energy_tmp, K, HV_tmp)
    HV_VP(:,:) = HV_VP(:,:) + HV_tmp(:,:)

    return
  end subroutine wave_vortex_production

! *****************************************************************************

  subroutine derive_viscous_factor(delta_time_until_b,  &
    &  viscous_factor, diffusive_factor)
    use Parameters, only : Prandtl_number
    implicit none
    double precision, intent(in) :: delta_time_until_b
    double precision, intent(out) ::  &
      &  viscous_factor(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(out) ::  &
      &  diffusive_factor(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    viscous_factor(:,:,:) = exp(decay_rate(:,:,:) * delta_time_until_b)
    diffusive_factor(:,:,:) = exp(decay_rate(:,:,:) * delta_time_until_b  &
      &  / Prandtl_number)

    return
  end subroutine derive_viscous_factor

! *****************************************************************************

  subroutine viscous_parameterization(K, U, T)
    use mpi, only : MPI_COMM_WORLD
    use mpi, only : MPI_DOUBLE_PRECISION
    use mpi, only : MPI_SUM
    use Parallel, only : my_rank
    use Parallel, only : rank_1
    use Geometry, only : wavenumber_moving_frame_local
    use Geometry, only : K1_truncated
    use Geometry, only : K2_truncated
    use Geometry, only : K3_truncated
    use Geometry, only : get_mask
    use Parameters, only : LES_on
    use Parameters, only : LES_coefficient_H
    use Parameters, only : LES_coefficient_V
    use Parameters, only : viscosity
    implicit none
    double precision, intent(in) ::  &
      &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(in) ::  &
      &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(in) ::  &
      &  T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: E_on_shell_H_local
    double precision :: E_on_shell_H
    double precision :: E_on_shell_V_local
    double precision :: E_on_shell_V
    double precision :: scaled_wavenumber_H_tmp
    double precision :: scaled_wavenumber_V_tmp
    double precision ::  &
      &  scaled_wavenumber_H(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision ::  &
      &  scaled_wavenumber_V(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    integer :: mask(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision :: factor
    double precision :: thickness_factor = 0.1d0
    double precision :: thickness
    double precision :: KH_max
    double precision :: KV_max
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: i1_global

    integer :: ierr

    if (.not. allocated(decay_rate)) then
      allocate(decay_rate(0:n1C_local-1, 0:n2C_local-1, 0:n3-1))
    endif
    if (.not. LES_on) then
      decay_rate(:,:,:) = viscosity  &
      &  * (K(1,:,:,:)**2 + K(2,:,:,:)**2 + K(3,:,:,:)**2)
      return
    endif

    E_on_shell_H_local = 0.d0
    E_on_shell_V_local = 0.d0

    KV_max = K3_truncated
    KH_max = K1_truncated

    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local - 1
        do i1 = 0, n1C_local - 1
          i1_global = rank_1 * n1C_local + i1
          if (i1_global == 0) then
            factor = 1.d0
          else
            factor = 2.d0
          endif
          scaled_wavenumber_H_tmp  &
          &  = (wavenumber_moving_frame_local(1)%val(i1)  &
          &  / K1_truncated)**2  &
          &  + (wavenumber_moving_frame_local(2)%val(i2)  &
          &  / K2_truncated)**2
          scaled_wavenumber_H_tmp = sqrt(scaled_wavenumber_H_tmp)
          scaled_wavenumber_H(i1, i2, i3) = scaled_wavenumber_H_tmp

          scaled_wavenumber_V_tmp  &
          &  = wavenumber_moving_frame_local(3)%val(i3)  &
          &  / K3_truncated
          scaled_wavenumber_V(i1, i2, i3) = abs(scaled_wavenumber_V_tmp)

          if (scaled_wavenumber_H_tmp > (1.d0 - thickness_factor)  &
          &  .and. scaled_wavenumber_H_tmp <= 1.d0) then
            E_on_shell_H_local = E_on_shell_H_local  &
            &  + sum(abs(U(:, i1, i2, i3))**2) * 0.5d0 * factor
          endif

          if (scaled_wavenumber_V_tmp > (1.d0 - thickness_factor)  &
          &  .and. scaled_wavenumber_V_tmp <= 1.d0) then
            E_on_shell_V_local = E_on_shell_V_local  &
            &  + sum(abs(U(:, i1, i2, i3))**2) * 0.5d0 * factor
          endif

        enddo
      enddo
    enddo

    call MPI_Allreduce(E_on_shell_H_local, E_on_shell_H,  &
    &  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    thickness = thickness_factor * KH_max
    E_on_shell_H = E_on_shell_H / thickness

    call MPI_Allreduce(E_on_shell_V_local, E_on_shell_V,  &
    &  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    thickness = thickness_factor * KV_max
    E_on_shell_V = E_on_shell_V / thickness

    decay_coefficient_H = LES_coefficient_H * sqrt(E_on_shell_H / KH_max)
    decay_coefficient_V = LES_coefficient_V * sqrt(E_on_shell_V / KV_max)

    decay_rate(:,:,:) = (viscosity  &
    &  + decay_coefficient_H * scaled_wavenumber_H(:,:,:)**2  &
    &  + decay_coefficient_V * scaled_wavenumber_V(:,:,:)**4)  &
    &  * (K(1,:,:,:)**2 + K(2,:,:,:)**2 + K(3,:,:,:)**2)

    call get_mask(mask)
    decay_rate(:,:,:) = mask(:,:,:) * decay_rate(:,:,:)

    return
  end subroutine viscous_parameterization

end module Dynamics

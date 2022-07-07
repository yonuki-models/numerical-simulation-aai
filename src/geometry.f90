module Geometry
  use Parallel, only : total_process
  use Parallel, only : my_rank
  use Parallel, only : comm_1
  use Parallel, only : comm_2
  use Parallel, only : rank_1
  use Parallel, only : rank_2
  implicit none
  private

  type :: Real_1D_Array
    double precision, pointer, public :: val(:)
  end type Real_1D_Array
  type(Real_1D_Array), public :: wavenumber_moving_frame_local(1:3)

  type :: State_3D_Complex
    complex(kind(0d0)), pointer, public :: velocity(:,:,:,:)
    complex(kind(0d0)), pointer, public :: temperature(:,:,:)
    contains
    procedure, public, pass :: calculate_kinetic_energy
    procedure, public, pass :: calculate_potential_energy
  end type State_3D_Complex
  type(State_3D_Complex), public :: state

  type :: Energy_Variables
    double precision, pointer, public :: kinetic_energy_production(:,:,:)
    double precision, pointer, public :: potential_energy_production(:,:,:)
    double precision, pointer, public :: energy_conversion(:,:,:)
    double precision, pointer, public :: kinetic_energy_dissipation(:,:,:)
    double precision, pointer, public :: potential_energy_dissipation(:,:,:)
    double precision, pointer, public :: kinetic_energy_transfer(:,:,:)
    double precision, pointer, public :: potential_energy_transfer(:,:,:)
  end type Energy_Variables
  type(Energy_Variables), public :: energy

  type :: Horizontal_Spectrum
    double precision, pointer, public :: kinetic_energy_production(:,:)
    double precision, pointer, public :: potential_energy_production(:,:)
    double precision, pointer, public :: energy_conversion(:,:)
    double precision, pointer, public :: kinetic_energy_dissipation(:,:)
    double precision, pointer, public :: potential_energy_dissipation(:,:)
    double precision, pointer, public :: kinetic_energy_transfer(:,:)
    double precision, pointer, public :: potential_energy_transfer(:,:)
  end type Horizontal_Spectrum
  type(Horizontal_Spectrum), public :: spectrum_H

  type :: HV_Spectrum
    double precision, pointer, public :: kinetic_energy_production(:,:)
    double precision, pointer, public :: potential_energy_production(:,:)
    double precision, pointer, public :: energy_conversion(:,:)
    double precision, pointer, public :: kinetic_energy_dissipation(:,:)
    double precision, pointer, public :: potential_energy_dissipation(:,:)
    double precision, pointer, public :: kinetic_energy_transfer(:,:)
    double precision, pointer, public :: potential_energy_transfer(:,:)
    double precision, pointer, public :: wave_production_1(:,:)
    double precision, pointer, public :: wave_production_2(:,:)
    double precision, pointer, public :: vortex_production(:,:)
    double precision, pointer, public :: kinetic_flux_convergence(:,:)
    double precision, pointer, public :: potential_flux_convergence(:,:)
  end type HV_Spectrum
  type(HV_Spectrum), public :: spectrum_HV

  double precision, public :: K1_truncated
  double precision, public :: K2_truncated
  double precision, public :: K3_truncated
  double precision, public :: delta_K(1:3)

  integer, public :: number_of_grid(1:3)
  integer, public :: division(1:2)
  integer, public :: n1
  integer, public :: n2
  integer, public :: n3
  integer, public :: n1_half
  integer, public :: n2_half
  integer, public :: n3_half
  integer, public :: n1C_local
  integer, public :: n2C_local
  integer, public :: n1R_local
  integer, public :: n2R_local
  integer, public :: n3R_local

  integer, public :: j1C_begin
  integer, public :: j1C_end
  integer, public :: j2C_begin
  integer, public :: j2C_end
  integer, public :: j2R_begin
  integer, public :: j2R_end
  integer, public :: j3R_begin
  integer, public :: j3R_end

  integer, public :: n1_truncated
  integer, public :: n2_truncated
  integer, public :: n3_truncated

  logical, public :: restart
  logical, public :: save_initial_data
  integer, public :: file_number_initial

  integer, public :: n1_low = 10
  integer, public :: n2_low = 10
  integer, public :: n3_low = 10

  double precision, public :: domain_length(1:3)

  ! For spectrum analysis
  double precision, allocatable, public :: energy_mask(:,:,:)
  double precision, public :: KH_max
  double precision, public :: KV_max
  double precision, public :: delta_KH
  double precision, public :: delta_KV
  integer, public :: NH_max
  integer, public :: NH_min
  integer, public :: NV_max

  ! For Richardson number analysis
  double precision, public :: Ri_max = 10.d0
  double precision, public :: Ri_min = - 5.d0
  integer, public :: N_Ri = 1500

  ! For Thorpe displacements analysis
  integer, public :: N_Thorpe

  public State_3D_Complex
  public Energy_Variables
  public Horizontal_Spectrum
  public HV_Spectrum
  public initialize_geometry
  public initialize_spectrum
  public dealias_truncate
  public get_mask
  public calculate_wavenumbers
  public incompressible_projection
  public reset_energy_data
  public finalize_geometry
  public redistribute_real_data

contains

! *****************************************************************************

  subroutine initialize_geometry()
    use Parallel, only : total_process
    use Misk, only : assert
    implicit none
    double precision, parameter :: PI = acos(-1.d0)
    integer :: number_of_grid_half(1:3)

    integer :: i1
    integer :: i1_global
    integer :: i2
    integer :: i2_global
    integer :: i3

    n1 = number_of_grid(1)
    n2 = number_of_grid(2)
    n3 = number_of_grid(3)
    number_of_grid_half(1:3) = number_of_grid(1:3) / 2
    n1_half = n1 / 2
    n2_half = n2 / 2
    n3_half = n3 / 2

    call assert(product(division)==total_process,  &
      &  'Product of division should be total process.')

    call assert(mod(n2, division(1))==0 .and. mod(n3, division(2))==0,  &
      &  'Number of grid (2, 3) should be multiple of division(1, 2).')

    call assert(mod(n1, division(2))==0,  &
      &  'Number of grid 1 should be multiple of division(2).')

    n1C_local = n1_half / division(1) + 1
    n2C_local = (n2 - 1) / division(2) + 1
    n1R_local = (n1 - 1) / division(2) + 1
    n2R_local = (n2 - 1) / division(1) + 1
    n3R_local = (n3 - 1) / division(2) + 1

    j1C_begin = min(n1C_local * rank_1, n1_half+1)
    j1C_end = min(n1C_local * (rank_1 + 1), n1_half+1)
    j2C_begin = min(n2C_local * rank_2, n2)
    j2C_end = min(n2C_local * (rank_2 + 1), n2)

    j2R_begin = min(n2R_local * rank_1, n2)
    j2R_end = min(n2R_local * (rank_1 + 1), n2)
    j3R_begin = min(n3R_local * rank_2, n3)
    j3R_end = min(n3R_local * (rank_2 + 1), n3)

    delta_K(1:3) = 2 * PI / domain_length(1:3)

    allocate(wavenumber_moving_frame_local(1)%val(0:n1C_local-1), source=0.d0)
    allocate(wavenumber_moving_frame_local(2)%val(0:n2C_local-1), source=0.d0)
    allocate(wavenumber_moving_frame_local(3)%val(0:n3-1), source=0.d0)

    allocate(  &
      &  state%velocity(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1),  &
      &  state%temperature(0:n1C_local-1, 0:n2C_local-1, 0:n3-1),  &
      &  source=(0.d0, 0.d0))

    allocate( &
      &  energy%kinetic_energy_production(  &
        &  0:n1C_local-1, 0:n2C_local-1, 0:n3-1), &
      &  energy%potential_energy_production(  &
        &  0:n1C_local-1, 0:n2C_local-1, 0:n3-1), &
      &  energy%energy_conversion(  &
        &  0:n1C_local-1, 0:n2C_local-1, 0:n3-1), &
      &  energy%kinetic_energy_dissipation(  &
        &  0:n1C_local-1, 0:n2C_local-1, 0:n3-1), &
      &  energy%potential_energy_dissipation(  &
        &  0:n1C_local-1, 0:n2C_local-1, 0:n3-1), &
      &  energy%kinetic_energy_transfer(  &
        &  0:n1C_local-1, 0:n2C_local-1, 0:n3-1), &
      &  energy%potential_energy_transfer(  &
        &  0:n1C_local-1, 0:n2C_local-1, 0:n3-1), &
      &  source=0.d0)

    allocate(energy_mask(0:n1C_local-1, 0:n2C_local-1, 0:n3-1),  &
      &  source=0.d0)

    do i1 = 0, n1C_local-1
      i1_global = rank_1 * n1C_local + i1
      wavenumber_moving_frame_local(1)%val(i1) =  &
        &  delta_K(1) * i1_global
    enddo

    do i2 = 0, n2C_local-1
      i2_global = rank_2 * n2C_local + i2
      if (i2_global <= n2_half) then
        wavenumber_moving_frame_local(2)%val(i2) =  &
         &  delta_K(2) * i2_global
      else
        wavenumber_moving_frame_local(2)%val(i2) =  &
         &  delta_K(2) * (i2_global - n2)
      endif
    enddo

    do i3 = 0, n3-1
      if (i3 <= n3_half) then
        wavenumber_moving_frame_local(3)%val(i3) =  &
         &  delta_K(3) * i3
      else
        wavenumber_moving_frame_local(3)%val(i3) =  &
         &  delta_K(3) * (i3 - n3)
      endif
    enddo

    n1_truncated = int(n1_half * 2.d0 / 3)
    n2_truncated = int(n2_half * 2.d0 / 3)
    n3_truncated = int(n3_half * 2.d0 / 3)

    K1_truncated = n1_truncated * delta_K(1)
    K2_truncated = n2_truncated * delta_K(2)
    K3_truncated = n3_truncated * delta_K(3)

    NH_max = n1_truncated
    NV_max = n3_truncated

    delta_KH = delta_K(1)
    delta_KV = delta_K(3)

    N_Thorpe = n3 * 2 - 1

    return
  end subroutine initialize_geometry

! *****************************************************************************

  subroutine initialize_spectrum()
    implicit none

    allocate( &
      &  spectrum_H%kinetic_energy_production(0:NH_max, 0:NH_min*2), &
      &  spectrum_H%potential_energy_production(0:NH_max, 0:NH_min*2), &
      &  spectrum_H%energy_conversion(0:NH_max, 0:NH_min*2), &
      &  spectrum_H%kinetic_energy_dissipation(0:NH_max, 0:NH_min*2), &
      &  spectrum_H%potential_energy_dissipation(0:NH_max, 0:NH_min*2), &
      &  spectrum_H%kinetic_energy_transfer(0:NH_max, 0:NH_min*2), &
      &  spectrum_H%potential_energy_transfer(0:NH_max, 0:NH_min*2), &
      &  source=0.d0)

    allocate( &
      &  spectrum_HV%kinetic_energy_production(0:NH_max, 0:NV_max), &
      &  spectrum_HV%potential_energy_production(0:NH_max, 0:NV_max), &
      &  spectrum_HV%energy_conversion(0:NH_max, 0:NV_max), &
      &  spectrum_HV%kinetic_energy_dissipation(0:NH_max, 0:NV_max), &
      &  spectrum_HV%potential_energy_dissipation(0:NH_max, 0:NV_max), &
      &  spectrum_HV%kinetic_energy_transfer(0:NH_max, 0:NV_max), &
      &  spectrum_HV%potential_energy_transfer(0:NH_max, 0:NV_max), &
      &  spectrum_HV%wave_production_1(0:NH_max, 0:NV_max), &
      &  spectrum_HV%wave_production_2(0:NH_max, 0:NV_max), &
      &  spectrum_HV%vortex_production(0:NH_max, 0:NV_max), &
      &  spectrum_HV%kinetic_flux_convergence(0:NH_max, 0:NV_max), &
      &  spectrum_HV%potential_flux_convergence(0:NH_max, 0:NV_max), &
      &  source=0.d0)

    return
  end subroutine initialize_spectrum

! *****************************************************************************

  subroutine dealias_truncate(U, T)
    implicit none
    complex(kind(0d0)), intent(inout) ::  &
      &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(inout) ::  &
      &  T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    integer :: mask(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    integer :: i

    call get_mask(mask)

    do i = 1, 3
      U(i,:,:,:) = U(i,:,:,:) * mask(:,:,:)
    enddo
    T(:,:,:) = T(:,:,:) * mask(:,:,:)

    return
  end subroutine dealias_truncate

! *****************************************************************************

  subroutine get_mask(mask)
    implicit none
    integer, intent(out) :: mask(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    integer :: i1
    integer :: i2
    integer :: i3
    double precision :: scaled_wavenumber_V
    double precision :: scaled_wavenumber_H

    mask(:,:,:) = 1

    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local - 1
        do i1 = 0, n1C_local - 1
          scaled_wavenumber_V  &
          &  = abs(wavenumber_moving_frame_local(3)%val(i3))  &
          &  / K3_truncated

          scaled_wavenumber_H  &
          &  = (wavenumber_moving_frame_local(1)%val(i1)  &
          &  / K1_truncated)**2  &
          &  + (wavenumber_moving_frame_local(2)%val(i2)  &
          &  / K2_truncated)**2
          scaled_wavenumber_H = sqrt(scaled_wavenumber_H)
          if (scaled_wavenumber_V > 1.d0 .or. scaled_wavenumber_H > 1.d0) then
            mask(i1, i2, i3) = 0
          endif
        enddo
      enddo
    enddo

    return
  end subroutine get_mask

! *****************************************************************************

  subroutine calculate_wavenumbers(A, K)
    implicit none
    double precision, intent(in) :: A(1:3, 1:3)
    double precision, intent(out) ::  &
      &  K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: tmp
    double precision :: KM(1:3) = 0.d0
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: i
    integer :: j

    do i3 = 0, n3 - 1
      do i2 = 0, n2C_local - 1
        do i1 = 0, n1C_local - 1
          KM(1) = wavenumber_moving_frame_local(1)%val(i1)
          KM(2) = wavenumber_moving_frame_local(2)%val(i2)
          KM(3) = wavenumber_moving_frame_local(3)%val(i3)
          do i = 1, 3
            tmp = 0.d0
            do j = 1, 3
              tmp = tmp + KM(j) * A(j,i)
            enddo
            K(i, i1, i2, i3) = tmp
          enddo
        enddo
      enddo
    enddo

    return
  end subroutine calculate_wavenumbers

! *****************************************************************************

  subroutine incompressible_projection(K_a, K_b, U, T, PK,  &
      &  viscous_factor, diffusive_factor)
    implicit none
    double precision, intent(in) ::  &
      &  K_a(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in) ::  &
      &  K_b(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(inout) ::  &
      &  U(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(inout) ::  &
      &  T(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)), intent(out) ::  &
      &  PK(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in), optional ::  &
      &  viscous_factor(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, intent(in), optional ::  &
      &  diffusive_factor(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision, allocatable :: K_squared(:,:,:)
    complex(kind(0d0)), allocatable :: P(:,:,:)
    integer :: i

    allocate(K_squared(0:n1C_local-1, 0:n2C_local-1, 0:n3-1), source=0.d0)
    allocate(P(0:n1C_local-1, 0:n2C_local-1, 0:n3-1), source=(0.d0, 0.d0))

    do i = 1, 3
      K_squared(:,:,:) = K_squared(:,:,:) + K_a(i,:,:,:) * K_b(i,:,:,:)
    enddo
    if (my_rank == 0) K_squared(0,0,0) = 1.0d0

    do i = 1, 3
      P(:,:,:) = P(:,:,:) + K_b(i,:,:,:) * U(i,:,:,:)
    enddo
    P(:,:,:) = P(:,:,:) / K_squared(:,:,:)

    do i = 1, 3
      PK(i,:,:,:) = K_a(i,:,:,:) * P(:,:,:)
    enddo

    if (present(viscous_factor)) then
      do i = 1, 3
        U(i,:,:,:) = (U(i,:,:,:) - PK(i,:,:,:)) / viscous_factor(:,:,:)
      enddo
    else
      U(:,:,:,:) = U(:,:,:,:) - PK(:,:,:,:)
    endif

    if (present(diffusive_factor)) then
      T(:,:,:) = T(:,:,:) / diffusive_factor(:,:,:)
    endif

    deallocate(K_squared)
    deallocate(P)

    return
  end subroutine incompressible_projection

! *****************************************************************************

  subroutine calculate_kinetic_energy(self, kinetic_energy)
    implicit none
    class(State_3D_Complex), intent(in) :: self
    double precision, intent(out) ::  &
      &  kinetic_energy(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    integer :: i

    kinetic_energy(:,:,:) = 0.d0
    do i = 1, 3
      kinetic_energy(:,:,:) = kinetic_energy(:,:,:)  &
        &  + abs(self%velocity(i,:,:,:))**2 * 0.5d0
    end do

    return
  end subroutine calculate_kinetic_energy

! *****************************************************************************

  subroutine calculate_potential_energy(self, potential_energy)
    implicit none
    class(State_3D_Complex), intent(in) :: self
    double precision, intent(out) ::  &
      &  potential_energy(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    potential_energy = abs(self%temperature(:,:,:))**2 * 0.5d0

    return
  end subroutine calculate_potential_energy

! *****************************************************************************

  subroutine reset_energy_data()
    implicit none

    spectrum_HV%kinetic_energy_production(:,:) = 0.d0
    spectrum_HV%potential_energy_production(:,:) = 0.d0
    spectrum_HV%energy_conversion(:,:) = 0.d0
    spectrum_HV%kinetic_energy_dissipation(:,:) = 0.d0
    spectrum_HV%potential_energy_dissipation(:,:) = 0.d0
    spectrum_HV%kinetic_energy_transfer(:,:) = 0.d0
    spectrum_HV%potential_energy_transfer(:,:) = 0.d0
    spectrum_HV%wave_production_1(:,:) = 0.d0
    spectrum_HV%wave_production_2(:,:) = 0.d0
    spectrum_HV%vortex_production(:,:) = 0.d0
    spectrum_HV%kinetic_flux_convergence(:,:) = 0.d0
    spectrum_HV%potential_flux_convergence(:,:) = 0.d0

    spectrum_H%kinetic_energy_production(:,:) = 0.d0
    spectrum_H%potential_energy_production(:,:) = 0.d0
    spectrum_H%energy_conversion(:,:) = 0.d0
    spectrum_H%kinetic_energy_dissipation(:,:) = 0.d0
    spectrum_H%potential_energy_dissipation(:,:) = 0.d0
    spectrum_H%kinetic_energy_transfer(:,:) = 0.d0
    spectrum_H%potential_energy_transfer(:,:) = 0.d0

    energy%kinetic_energy_production(:,:,:) = 0.d0
    energy%potential_energy_production(:,:,:) = 0.d0
    energy%energy_conversion(:,:,:) = 0.d0
    energy%kinetic_energy_dissipation(:,:,:) = 0.d0
    energy%potential_energy_dissipation(:,:,:) = 0.d0
    energy%kinetic_energy_transfer(:,:,:) = 0.d0
    energy%potential_energy_transfer(:,:,:) = 0.d0

    return
  end subroutine reset_energy_data

! *****************************************************************************

  subroutine finalize_geometry()
    implicit none
    integer :: i

    do i = 1, 3
      deallocate(wavenumber_moving_frame_local(i)%val)
    enddo

    deallocate(state%velocity)
    deallocate(state%temperature)

    deallocate(energy%kinetic_energy_production)
    deallocate(energy%potential_energy_production)
    deallocate(energy%energy_conversion)
    deallocate(energy%kinetic_energy_dissipation)
    deallocate(energy%potential_energy_dissipation)
    deallocate(energy%kinetic_energy_transfer)
    deallocate(energy%potential_energy_transfer)

    return
  end subroutine finalize_geometry

! *****************************************************************************

  subroutine redistribute_real_data(R, R_redist)
    use mpi, only : MPI_DOUBLE_PRECISION
    implicit none
    double precision, intent(in) ::  &
      &  R(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    double precision, intent(out) ::  &
      &  R_redist(0:n1R_local-1, 0:n2R_local-1, 0:n3-1)
    double precision :: R_trans(0:n3R_local-1, 0:n1-1)
    double precision :: R_tmp(0:n3R_local-1, 0:n1R_local-1, 0:division(2)-1)
    integer :: i1
    integer :: i2
    integer :: j
    integer :: ierr

    do i2 = 0, n2R_local - 1
      R_trans(0:n3R_local-1, 0:n1-1)  &
        &  = transpose(R(0:n1-1, i2, 0:n3R_local-1))

      call MPI_Alltoall(R_trans(0, 0), n3R_local * n1R_local,  &
        &  MPI_DOUBLE_PRECISION, R_tmp(0, 0, 0), n3R_local * n1R_local,  &
        &  MPI_DOUBLE_PRECISION, comm_2, ierr)

      do j = 0, division(2)-1
        do i1 = 0, n1R_local - 1
          R_redist(i1, i2, j*n3R_local:(j+1)*n3R_local-1)  &
            &  = R_tmp(:, i1, j)
        enddo
      enddo
    enddo

    return
  end subroutine redistribute_real_data

end module Geometry

module Parameters
  implicit none
  private
  ! Basic parameters
  double precision, public :: buoyancy_frequency = 1.d0
  double precision, public :: angular_velocity(1:3) = 0.d0

  ! Dimensionless parameters
  double precision, public :: initial_noise_level = 1.d4

  ! Dissipation parameters
  logical, public :: LES_on = .True.
  double precision, public :: LES_coefficient_H = 0.5d0
  double precision, public :: LES_coefficient_V = 0.5d0
  double precision, public :: Prandtl_number = 1.d0
  double precision, public :: viscosity = 0.5d0

  logical, public :: wave_on

contains

end module Parameters

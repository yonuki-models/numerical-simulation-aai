module Inertia_Gravity_wave
  use Parameters, only : buoyancy_frequency
  use Parameters, only : angular_velocity
  implicit none
  private

  double precision, public :: frequency_divided_by_N
  double precision, public :: max_shear_divided_by_N
  double precision, public :: output_interval_wave_period
  double precision, public :: time_end_wave_period

  double precision, public :: frequency
  double precision, public :: Omega
  double precision, public :: max_shear
  double precision, public :: max_temperature_gradient
  double precision, public :: initial_phase = 0.d0

  double precision, public :: rotation_tensor(1:3, 1:3)
  double precision, public :: transform_tensor(1:3, 1:3)

  public :: wave_initial_setting
  public :: derive_temperature_gradient
  public :: derive_transform_tensor
  public :: derive_velocity_gradient

contains

! *****************************************************************************

  subroutine wave_initial_setting()
    use Time_Module, only : output_interval
    use Time_Module, only : time_end
    implicit none
    double precision, parameter :: PI = acos(-1.d0)
    double precision :: cos_a = 1.d0
    double precision :: sin_a = 0.d0
    double precision :: cos_b = 1.d0
    double precision :: sin_b = 0.d0
    double precision :: wave_period
    integer :: i

    max_shear = max_shear_divided_by_N * buoyancy_frequency
    frequency = frequency_divided_by_N * buoyancy_frequency

    angular_velocity(1) = 0.d0
    angular_velocity(2) = 0.d0

    cos_b = sqrt(frequency**2 - 4 * angular_velocity(3)**2  &
      &  / (buoyancy_frequency**2 - 4 * angular_velocity(3)**2))
    sin_b = sqrt(1 - cos_b**2)

    call derive_rotation_tensor(cos_a, sin_a, cos_b, sin_b)

    Omega = 0.d0
    do i = 1, 3
      Omega = Omega + angular_velocity(i) * rotation_tensor(1,i)
    enddo

    max_temperature_gradient = buoyancy_frequency * cos_b / frequency  &
      &  * max_shear

    ! Time settings
    wave_period = 2 * PI / frequency
    output_interval = wave_period * output_interval_wave_period
    time_end = wave_period * time_end_wave_period

    return
  end subroutine wave_initial_setting

! *****************************************************************************

  subroutine derive_rotation_tensor(cos_a, sin_a, cos_b, sin_b)
    implicit none
    double precision, intent(in) :: cos_a
    double precision, intent(in) :: sin_a
    double precision, intent(in) :: cos_b
    double precision, intent(in) :: sin_b

    rotation_tensor(1, :) = (/  cos_a * cos_b,  sin_a * cos_b, sin_b /)
    rotation_tensor(2, :) = (/         -sin_a,          cos_a,  0.d0 /)
    rotation_tensor(3, :) = (/ -cos_a * sin_b, -sin_a * sin_b, cos_b /)

    return
  end subroutine derive_rotation_tensor

! *****************************************************************************

  subroutine derive_transform_tensor(time, A)
    implicit none
    double precision, intent(in) :: time
    double precision, intent(out) :: A(1:3, 1:3)
    double precision :: D_integrate(1:3, 1:3)
    double precision :: S_integrate(1:3)
    double precision :: S_a_integrate
    double precision :: S_b_integrate
    integer :: i

    S_a_integrate = max_shear / frequency  &
      &  * sin(frequency * time + initial_phase)
    S_b_integrate = 2 * Omega / frequency**2 * max_shear  &
      &  * cos(frequency * time + initial_phase)

    do i = 1, 3
      S_integrate(i) = rotation_tensor(2,i) * S_a_integrate  &
        &  + rotation_tensor(3,i) * S_b_integrate
    enddo

    do i = 1, 3
      D_integrate(i,:) = S_integrate(i) * rotation_tensor(1,:)
    enddo

    A(:,:) = rotation_tensor(:,:) - D_integrate(:,:)

    return
  end subroutine derive_transform_tensor

! *****************************************************************************

  subroutine derive_velocity_gradient(time, D)
    implicit none
    double precision, intent(in) :: time
    double precision, intent(out) :: D(1:3, 1:3)
    double precision :: S(1:3)
    double precision :: S_a
    double precision :: S_b
    integer :: i

    S_a = max_shear * cos(frequency * time + initial_phase)
    S_b = - 2 * Omega / frequency * max_shear  &
      &  * sin(frequency * time + initial_phase)

    do i = 1, 3
      S(i) = rotation_tensor(2,i) * S_a + rotation_tensor(3,i) * S_b
    enddo

    do i = 1, 3
      D(i, 1:3) = S(i) * rotation_tensor(1, 1:3)
    enddo

    return
  end subroutine derive_velocity_gradient

! *****************************************************************************

  subroutine derive_temperature_gradient(time, M)
    implicit none
    double precision, intent(in) :: time
    double precision, intent(out) :: M(1:3)
    double precision :: M_absolute

    M_absolute = max_temperature_gradient  &
      &  * sin(frequency * time + initial_phase)

    M(1:3) = M_absolute * rotation_tensor(1, 1:3)

    return
  end subroutine derive_temperature_gradient

end module Inertia_Gravity_wave

module External_Field
  implicit none
  private

  double precision, public :: alpha_plus
  double precision, public :: alpha_minus
  integer, public :: output_during_a_period

  integer, public :: sign_vortex !  1 for anticyclonic vortex
                                 ! -1 for cyclonic vortex
  double precision, public :: frequency
  double precision, public :: e
  double precision, public :: transform_tensor(1:3, 1:3)

  public initial_setting
  public setting_output_time
  public derive_transform_tensor
  public derive_velocity_gradient
  public derive_temperature_gradient

contains

! *****************************************************************************

  subroutine initial_setting()
    use misk, only : assert
    use Geometry, only : n1_truncated
    use Geometry, only : n2_truncated
    use Geometry, only : n3_truncated
    use Geometry, only : NH_max
    use Geometry, only : NH_min
    use Geometry, only : NV_max
    implicit none

    call assert(alpha_plus * alpha_minus > 0.d0,  &
      &  'alpha_plus and alpha_minus mush have a same sign')

    sign_vortex = int(dsign(1.d0, alpha_plus))
    frequency = sqrt(alpha_plus * alpha_minus)
    e = sqrt(alpha_plus / alpha_minus)

    NH_max = int(sqrt(e) * n1_truncated + 0.5d0)
    NH_min = int(n2_truncated / sqrt(e) + 0.5d0)
    NV_max = n3_truncated

    return
  end subroutine initial_setting

! *****************************************************************************

  subroutine setting_output_time()
    use Time_Module, only : output_interval
    implicit none
    double precision :: PI = acos(-1.d0)

    output_interval = 2 * PI / frequency / output_during_a_period

    return
  end subroutine setting_output_time

! *****************************************************************************

  subroutine derive_transform_tensor(time, A)
    implicit none
    double precision, intent(in) :: time
    double precision, intent(out) :: A(1:3, 1:3)
    double precision :: root_e
    double precision :: re_ro_e
    double precision :: cos_wt
    double precision :: sin_wt

    root_e = sqrt(e)
    re_ro_e = 1.d0 / sqrt(e)
    cos_wt = cos(frequency * time)
    sin_wt = sin(frequency * time)

    A(1, :) = (/ re_ro_e * cos_wt, - sign_vortex * root_e * sin_wt, 0.d0 /)
    A(2, :) = (/ sign_vortex * re_ro_e * sin_wt, root_e * cos_wt, 0.d0 /)
    A(3, :) = (/ 0.d0, 0.d0, 1.d0 /)

    return
  end subroutine derive_transform_tensor

! *****************************************************************************

  subroutine derive_velocity_gradient(time, D)
    use misk, only : assert
    implicit none
    double precision, intent(in) :: time
    double precision, intent(out) :: D(1:3, 1:3)

    call assert(time >= 0.d0, 'time should be a non-negative number.')

    D(:,:) = 0.d0
    D(1,2) = alpha_plus
    D(2,1) = - alpha_minus

    return
  end subroutine derive_velocity_gradient

! *****************************************************************************

  subroutine derive_temperature_gradient(time, M)
    use misk, only : assert
    implicit none
    double precision, intent(in) :: time
    double precision, intent(out) :: M(1:3)

    call assert(time >= 0.d0, 'time shoule be a non-negative number.')

    M(1:3) = 0.d0

    return
  end subroutine derive_temperature_gradient

end module External_Field

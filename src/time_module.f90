module Time_Module
  implicit none
  private

  double precision, public :: delta_time_max
  double precision, public :: output_interval
  double precision, public :: time_left_until_output
  double precision, public :: current_time = 0.d0
  double precision, public :: time_end
  integer, public :: current_step = 0
  integer, public :: file_number = 0
  logical, public :: output

  integer, public :: output_state_interval

  public :: initialize_time
  public :: evaluate_delta_time
  public :: continuation

contains

! *****************************************************************************

  subroutine initialize_time()
    implicit none

    output = .False.
    time_left_until_output = output_interval
    file_number = file_number + 1

    return
  end subroutine initialize_time

! *****************************************************************************

  double precision function evaluate_delta_time(delta_time_CFL)
    implicit none
    double precision, intent(in) :: delta_time_CFL
    double precision :: delta_time_tmp

    delta_time_tmp = min(delta_time_max, delta_time_CFL)
    if (delta_time_tmp < time_left_until_output) then
      evaluate_delta_time = delta_time_tmp
      time_left_until_output = time_left_until_output - delta_time_tmp
    else
      evaluate_delta_time = time_left_until_output
      output = .True.
      time_left_until_output = output_interval
    endif

    return
  end function evaluate_delta_time

! *****************************************************************************

  logical function continuation()
    implicit none

    continuation = (current_time < time_end)

    return
  end function continuation

end module Time_Module

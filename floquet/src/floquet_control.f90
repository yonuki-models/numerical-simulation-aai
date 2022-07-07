module Floquet_Control
  implicit none
  private

  integer, parameter :: UNIT_FILE_DEFAULT = 10
  double precision, public :: buoyancy_frequency
  double precision, public :: angular_velocity(1:3) = (/ 0.d0, 0.d0, 0.5d0 /)
  logical, public :: Output_max
  logical, public :: N_log_scale
  logical, public :: M_log_scale
  integer, public :: N_time
  integer, public :: NN
  integer, public :: NM
  integer, public :: NR
  integer, public :: NE
  double precision, public :: N_min
  double precision, public :: N_max
  double precision, public :: M_min
  double precision, public :: M_max
  double precision, public :: Ro_min
  double precision, public :: Ro_max
  double precision, public :: e_min
  double precision, public :: e_max

  character(len=200), public :: data_directory = "../data/floquet/"
  character(len=200) :: control_directory = "../control/"
  character(len=200) :: experiment_name = ""

  public :: read_control_file

contains

! *****************************************************************************

  subroutine read_control_file()
    implicit none
    intrinsic get_environment_variable
    double precision, parameter :: PI = acos(-1.d0)
    integer :: len
    integer :: status
    character(len=200) :: control_file_path = ''

    namelist /nml_parameter/ Output_max, N_log_scale, M_log_scale, N_time,  &
    &  NN, NM, NR, NE,  &
    &  N_min, N_max, M_min, M_max, Ro_min, Ro_max, e_min, e_max

    call get_environment_variable("experiment", status=status, length=len,  &
    &  value=experiment_name)

    control_file_path = trim(control_directory)//trim(experiment_name)//'.in'
    data_directory = "../data/"//trim(experiment_name)

    open(unit=UNIT_FILE_DEFAULT, file=control_file_path,  &
    &  status='old', action='read')

    read(unit=UNIT_FILE_DEFAULT, nml=nml_parameter)

    close(UNIT_FILE_DEFAULT)

    control_file_path =  &
    &  trim(data_directory)//'/'//trim(experiment_name)//'.out'

    open(unit=UNIT_FILE_DEFAULT, file=control_file_path,  &
    &  status='replace', action='write')

    write(unit=UNIT_FILE_DEFAULT, nml=nml_parameter)

    close(UNIT_FILE_DEFAULT)

    return
  end subroutine read_control_file

end module Floquet_Control

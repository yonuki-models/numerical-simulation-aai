module Control
  implicit none
  private

  character(len=100), public :: experiment_name = 'exDefault'

  double precision :: Coriolis_parameter
  double precision :: Rossby_number
  double precision :: ellipticity

  public :: read_experiment_name
  public :: read_control_file

contains

! *****************************************************************************

  subroutine read_experiment_name()
    use Input_Output, only : data_directory
    implicit none
    integer :: len
    integer :: status

    call get_environment_variable("experiment", status=status, length=len,  &
      &  value=experiment_name)
    call get_environment_variable("output_dir", status=status, length=len,  &
      &  value=data_directory)

    return
  end subroutine read_experiment_name

! *****************************************************************************

  subroutine read_control_file()
    use Parameters, only : buoyancy_frequency
    use Parameters, only : angular_velocity
    use Parameters, only : initial_noise_level
    use Parameters, only : LES_on
    use Parameters, only : LES_coefficient_H
    use Parameters, only : LES_coefficient_V
    use Parameters, only : Prandtl_number
    use Parameters, only : viscosity
    use Geometry, only : number_of_grid
    use Geometry, only : division
    use Geometry, only : domain_length
    use Geometry, only : restart
    use Geometry, only : save_initial_data
    use Geometry, only : file_number_initial
    use Time_Module, only : time_end
    use Time_Module, only : delta_time_max
    use Time_Module, only : output_interval
    use Time_Module, only : output_state_interval
    use Input_Output, only : control_directory
    use Input_Output, only : UNIT_FILE_DEFAULT
    use External_Field, only : alpha_plus
    use External_Field, only : alpha_minus
    use External_Field, only : output_during_a_period
    implicit none
    intrinsic get_environment_variable
    double precision, parameter :: PI = acos(-1.d0)
    integer :: len
    integer :: status
    character(len=200) :: control_file_path = ''

    namelist /nml_geometry/ number_of_grid, division, domain_length

    namelist /nml_time/ time_end, delta_time_max

    namelist /nml_parameter/ buoyancy_frequency, Coriolis_parameter,  &
      &  initial_noise_level

    namelist /nml_dissipation/ LES_on, Prandtl_number,  &
      &  viscosity, LES_coefficient_H, LES_coefficient_V

    namelist /nml_output/ output_interval, output_state_interval

    namelist /nml_restart/ restart, save_initial_data, file_number_initial

    namelist /nml_vortex/ Rossby_number, ellipticity, output_during_a_period

    call get_environment_variable("experiment", status=status, length=len,  &
      &  value=experiment_name)

    control_file_path =  &
      &  trim(control_directory)//'/'//trim(experiment_name)//'.in'

    open(unit=UNIT_FILE_DEFAULT, file=control_file_path,  &
      &  status='old', action='read')

    read(unit=UNIT_FILE_DEFAULT, nml=nml_geometry)
    read(unit=UNIT_FILE_DEFAULT, nml=nml_time)
    read(unit=UNIT_FILE_DEFAULT, nml=nml_parameter)
    read(unit=UNIT_FILE_DEFAULT, nml=nml_dissipation)
    read(unit=UNIT_FILE_DEFAULT, nml=nml_output)
    read(unit=UNIT_FILE_DEFAULT, nml=nml_restart)
    read(unit=UNIT_FILE_DEFAULT, nml=nml_vortex)

    close(UNIT_FILE_DEFAULT)

    angular_velocity = (/ 0.d0, 0.d0, Coriolis_parameter/2 /)
    alpha_plus = Rossby_number * Coriolis_parameter  &
      &  / (1 + (1 - ellipticity)**2)
    alpha_minus = Rossby_number * Coriolis_parameter * (1 - ellipticity)**2  &
    &  / (1 + (1 - ellipticity)**2)

    return
  end subroutine read_control_file

end module Control

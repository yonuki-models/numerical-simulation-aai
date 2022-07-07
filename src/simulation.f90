module Simulation
  use Parallel, only : my_rank
  use Input_Output, only : LEN_MESSAGE
  use Input_Output, only : print_main
  use Geometry, only : state
  use Geometry, only : energy
  use Geometry, only : spectrum_H
  use Geometry, only : spectrum_HV
  use Geometry, only : reset_energy_data
  use Time_Module, only : current_time
  use Time_Module, only : current_step
  use Time_Module, only : file_number
  use Log, only : log_standard_output
  implicit none
  private
  character(LEN=LEN_MESSAGE) :: message

  public :: initial_process
  public :: step_forward
  public :: final_process

contains

  subroutine initial_process()
    use Parallel, only : start_mpi
    use FFTW_Interface, only : fftw_init_parallel
    use Geometry, only : file_number_initial
    use Geometry, only : initialize_geometry
    use Geometry, only : initialize_spectrum
    use Geometry, only : restart
    use Geometry, only : save_initial_data
    use Geometry, only : number_of_grid
    use Geometry, only : division
    use Time_Module, only : initialize_time
    use Time_Module, only : output_interval
    use Control, only : read_experiment_name
    use Control, only : read_control_file
    use Log, only : start_log
    use Log, only : save_parameters
    use Input_Output, only : output_geometry
    use Input_Output, only : output_energy_sum
    use Input_Output, only : output_state
    use Input_Output, only : output_energy
    use Input_Output, only : output_snapshot_all
    use Input_Output, only : output_Richardson_distribution
    use Input_Output, only : output_Froude_spectrum
    use Input_Output, only : output_Thorpe_distribution
    use Input_Output, only : input_state
    use Input_Output, only : output_spectrum
    use Input_Output, only : initialize_reference_energy
    use Initial_Condition, only : white_noise
    use Initial_Condition, only : white_noise_low_modes
    use Initial_Condition, only : preconditioning
    use Initial_Condition, only : eliminate_vortical_mode
    use Initial_Condition, only : adjust_to_power_law
    use External_Field, only : initial_setting
    use External_Field, only : setting_output_time
    implicit none

    call start_log()
    call read_experiment_name()
    call read_control_file()

    call start_mpi(division)
    call initialize_geometry()
    call fftw_init_parallel(number_of_grid, division)
    call initial_setting()
    call initialize_spectrum()
    call setting_output_time()

    if (restart) then
      file_number = file_number_initial
      current_time = file_number * output_interval
      call input_state()
      call preconditioning()
      if (save_initial_data) then
        call output_geometry()
        call output_energy_sum()
        call output_spectrum()
        call output_Richardson_distribution()
        call output_Froude_spectrum()
      endif
      call initialize_reference_energy()
    else
      call white_noise_low_modes()
      call output_state()
      call output_energy()
      call output_energy_sum()
      call output_spectrum()
      call output_snapshot_all()
      call output_Richardson_distribution()
      call output_Froude_spectrum()
      call output_Thorpe_distribution()
      call output_geometry()
    endif

    if (my_rank == 0) call save_parameters()

    write(message, '("Integration is starting at file number: ", i4.4)')  &
      &  file_number
    call print_main(trim(message))
    call log_standard_output()
    call initialize_time()

    return
  end subroutine initial_process

! *****************************************************************************

  subroutine step_forward()
    use mpi, only : MPI_COMM_WORLD
    use Time_Module, only : output_state_interval
    use Time_Module, only : output
    use Initial_Condition, only : preconditioning
    use Input_Output, only : output_geometry
    use Input_Output, only : output_energy_sum
    use Input_Output, only : output_state
    use Input_Output, only : output_energy
    use Input_Output, only : output_snapshot_all
    use Input_Output, only : output_spectrum
    use Input_Output, only : output_Richardson_distribution
    use Input_Output, only : output_Froude_spectrum
    use Input_Output, only : output_Thorpe_distribution
    use Integration, only : runge_kutta_3
    implicit none
    double precision :: delta_time

    call runge_kutta_3(delta_time)
    current_time = current_time + delta_time
    current_step = current_step + 1

    if (output) then
      call preconditioning()

      if (mod(file_number, output_state_interval) == 0) then
        write(message, '("*** Output state variables ***")')
        call print_main(trim(message))
        call output_state()
        call output_energy()
      endif

      write(message, '("Output at file number: ", i6.6)') file_number
      call print_main(trim(message))
      call output_energy_sum()
      call output_spectrum()
      call output_snapshot_all()
      call output_Richardson_distribution()
      call output_Froude_spectrum()
      call output_Thorpe_distribution()
      call reset_energy_data()
      call output_geometry()
      call log_standard_output()

      file_number = file_number + 1
      output = .False.
    endif

    return
  end subroutine step_forward

! *****************************************************************************

  subroutine final_process()
    use Parallel, only : end_mpi
    use FFTW_Interface, only : fftw_finalize_parallel
    use Geometry, only : finalize_geometry
    implicit none

    call finalize_geometry()
    call fftw_finalize_parallel()

    write(message, '("Integration has finished normally")')
    call print_main(trim(message))
    call end_mpi()

    return
  end subroutine final_process

end module Simulation

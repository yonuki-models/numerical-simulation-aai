module Log
  use Parallel, only : my_rank
  use Input_Output, only : UNIT_STDOUT
  use Input_Output, only : data_directory
  use Geometry, only : state
  use Control, only : experiment_name
  implicit none
  private

  integer, parameter :: UNIT_LOG = 20

  real :: start_cpu_time
  integer :: start_system_count
  integer :: system_count
  integer :: system_count_per_second
  integer :: system_count_max

  integer :: count_cycle
  logical :: larger_than_start_count

  public :: start_log
  public :: log_standard_output
  public :: log_time
  public :: log_energy
  public :: save_parameters

contains

! *****************************************************************************

  subroutine start_log()
    implicit none

    call cpu_time(start_cpu_time)
    call system_clock(count=start_system_count,  &
      &  count_rate=system_count_per_second,  &
      &  count_max=system_count_max)
    count_cycle = 0
    larger_than_start_count = .True.

    return
  end subroutine start_log

! *****************************************************************************

  subroutine log_standard_output()
    implicit none

    if (my_rank == 0) call log_time()
    call log_energy()

    return
  end subroutine log_standard_output

! *****************************************************************************

  subroutine log_time()
    use Time_Module, only : current_time
    use Time_Module, only : time_end
    use Time_Module, only : current_step
    use Time_Module, only : file_number
    implicit none
    integer :: hour = 0
    integer :: minute = 0
    integer :: second = 0
    real :: cpu_time_ = 0.d0

    write(UNIT_STDOUT, '(A, I0)') 'Current step: ', current_step

    call system_clock(count=system_count)
    if (larger_than_start_count) then
      if (system_count < start_system_count) then
        count_cycle = count_cycle + 1
        larger_than_start_count = .False.
      endif
    else
      if (system_count > start_system_count) then
        larger_than_start_count = .True.
      endif
    endif
    system_count = system_count + system_count_max * count_cycle

    second = (system_count - start_system_count) /  &
      &  system_count_per_second
    minute = second / 60
    hour = minute / 60
    second = second - 60 * minute
    minute = minute - 60 * hour
    write(UNIT_STDOUT, '(A, I0, ":" I2.2, ":", I2.2)') 'Real time: ',  &
      &  hour, minute, second

    call cpu_time(cpu_time_)
    second = int(cpu_time_ - start_cpu_time)
    minute = second / 60
    hour = minute / 60
    second = second - 60 * minute
    minute = minute - 60 * hour
    write(UNIT_STDOUT, '(A, I0, ":" I2.2, ":", I2.2)') ' CPU time: ',  &
      &  hour, minute, second

    write(UNIT_STDOUT, '(A, f12.6)') 'Simulation time (days): ',  &
      &  current_time / 86400.d0
    write(UNIT_STDOUT, '(A, f12.6)') '                 Until: ',  &
      &  time_end / 86400.d0
    write(UNIT_STDOUT, '(A, I0)') 'File number: ', file_number

    return
  end subroutine log_time

! *****************************************************************************

  subroutine log_energy()
    use mpi, only : MPI_COMM_WORLD
    use Geometry, only : n1C_local, n2C_local, n3
    use Analysis, only : sum_real_data
    implicit none
    double precision :: sum_KE
    double precision :: sum_PE
    double precision :: kinetic_energy(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: potential_energy(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    integer :: ierr

    call state%calculate_kinetic_energy(kinetic_energy)
    call state%calculate_potential_energy(potential_energy)
    sum_KE = sum_real_data(kinetic_energy)
    sum_PE = sum_real_data(potential_energy)

    if (my_rank == 0) then
      write(UNIT_STDOUT, '(A, f24.12)') 'Total energy: ', sum_KE + sum_PE
      write(UNIT_STDOUT, *)
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return
  end subroutine log_energy

! *****************************************************************************

  subroutine save_parameters()
    use Parallel, only : total_process
    use Geometry, only : division
    use Parameters, only : buoyancy_frequency
    use Parameters, only : angular_velocity
    use Geometry, only : number_of_grid
    use Geometry, only : division
    use Geometry, only : domain_length
    use Geometry, only : n1_low
    use Geometry, only : n2_low
    use Geometry, only : n3_low
    use Geometry, only : N_Ri
    use Geometry, only : Ri_max
    use Geometry, only : Ri_min
    use Geometry, only : N_Thorpe
    use Time_Module, only : time_end
    use Time_Module, only : output_interval
    use Time_Module, only : output_state_interval
    use External_Field, only : alpha_plus
    use External_Field, only : alpha_minus
    use Input_Output, only : print_main
    implicit none
    character(len=120) :: control_file_output

    namelist /params/ number_of_grid, division, domain_length,  &
      &  time_end, output_interval, output_state_interval,  &
      &  buoyancy_frequency, angular_velocity, alpha_plus, alpha_minus,  &
      &  N_Ri, Ri_min, Ri_max, N_Thorpe

    control_file_output =  &
      &  trim(data_directory)//'/'//trim(experiment_name)//'.out'

    open(unit=unit_log, file=control_file_output,  &
        &  status='replace', form='formatted')
    write(unit_log, nml=params)
    close(unit_log)

    call print_main("Saved the parameter file: "//control_file_output)

    return
  end subroutine save_parameters

end module Log

program main
  use mpi, only : MPI_THREAD_FUNNELED
  use mpi, only : MPI_COMM_WORLD
  use mpi, only : MPI_DOUBLE_PRECISION
  use Parallel, only : total_process
  use Parallel, only : my_rank
  use Parallel, only : start_mpi
  use Parallel, only : end_mpi
  use Floquet_Control, only : read_experiment_name
  use Floquet_Control, only : read_control_file
  use Floquet_Control, only : experiment_name
  use Floquet_Control, only : M_max
  use Floquet_Control, only : M_min
  use External_Field, only : alpha_plus
  use External_Field, only : alpha_minus
  use Floquet_Solver, only : elliptic_instability_solver
  implicit none

  double precision, parameter :: PI = acos(-1.d0)

  integer, parameter :: NZ = 4
  integer, parameter :: NM = 1000
  integer, parameter :: NA = 1000
  integer, parameter :: UNIT_DATA = 10

  integer :: NM_local
  integer :: im
  integer :: im_local

  double precision :: M
  double precision :: log_M_min = -1.d0
  double precision :: log_M_max = 3.d0
  double precision :: log_M

  double precision :: growth_rate
  double precision :: growth_rate_global(1:NM)
  double precision, allocatable :: growth_rate_local(:)
  integer :: ierr

  character(200) :: growth_rate_file

  call start_mpi()

  call read_experiment_name()
  call read_control_file()

  NM_local = NM / total_process
  log_M_min = log(M_min)
  log_M_max = log(M_max)

  allocate(growth_rate_local(1:NM_local))

  do im = 1, NM_local
    im_local = NM_local * my_rank + im
    log_M = log_M_min + (log_M_max - log_M_min) / NM * im_local
    M = exp(log_M)
    call elliptic_instability_solver(alpha_plus, alpha_minus, M, growth_rate)
    growth_rate_local(im) = growth_rate
    write(*,*) im, M, growth_rate
  enddo

  call MPI_Gather(growth_rate_local(1), NM_local,  &
    &  MPI_DOUBLE_PRECISION, growth_rate_global(1), NM_local,  &
    &  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  if (my_rank .EQ. 0) then
    growth_rate_file = '../data/'//trim(experiment_name)//'/growth_rate.out'
    open(UNIT_DATA, file=growth_rate_file, access='direct',  &
      &  status='unknown', form='unformatted', recl=8*NM)
    write(UNIT_DATA, REC=1) growth_rate_global
    close(UNIT_DATA)
  end if

  call end_mpi()

  stop
end program main

program main
  use mpi, only : MPI_THREAD_FUNNELED
  use mpi, only : MPI_COMM_WORLD
  use mpi, only : MPI_DOUBLE_PRECISION
  use mpi, only : MPI_DOUBLE_COMPLEX
  use mpi, only : MPI_MAX
  use Parallel, only : total_process
  use Parallel, only : my_rank
  use Parallel, only : start_mpi
  use Parallel, only : end_mpi
  use Floquet_Control, only : buoyancy_frequency
  use Floquet_Control, only : angular_velocity
  use Floquet_Control, only : read_control_file
  use Floquet_Control, only : Output_max
  use Floquet_Control, only : N_log_scale
  use Floquet_Control, only : M_log_scale
  use Floquet_Control, only : NN
  use Floquet_Control, only : NM
  use Floquet_Control, only : NR
  use Floquet_Control, only : NE
  use Floquet_Control, only : N_max
  use Floquet_Control, only : N_min
  use Floquet_Control, only : M_max
  use Floquet_Control, only : M_min
  use Floquet_Control, only : Ro_max
  use Floquet_Control, only : Ro_min
  use Floquet_Control, only : e_max
  use Floquet_Control, only : e_min
  use Floquet_Control, only : data_directory
  use Floquet_Solver, only : elliptic_instability_solver
  implicit none

  integer, parameter :: NZ = 4
  integer, parameter :: UNIT_DATA = 10

  integer :: in
  double precision :: log_N_min
  double precision :: log_N_max
  double precision :: log_N

  integer :: NM_local
  integer :: im
  integer :: im_local
  integer :: ir
  integer :: ie

  double precision :: M
  double precision :: log_M_min
  double precision :: log_M_max
  double precision :: log_M

  double precision :: Ro

  double precision :: e

  double precision :: alpha_plus
  double precision :: alpha_minus

  double precision :: growth_rate
  double precision, allocatable :: growth_rate_global(:,:,:)
  double precision, allocatable :: growth_rate_local(:,:,:)
  double precision, allocatable :: growth_rate_max(:,:)
  double precision, allocatable :: growth_rate_max_local(:,:)

  double precision :: unstable_eigen_vector(1:NZ)
  double precision :: neutral_eigen_vector(1:NZ)
  double precision, allocatable :: unstable_global(:,:,:,:)
  double precision, allocatable :: unstable_local(:,:,:,:)
  double precision, allocatable :: neutral_global(:,:,:,:)
  double precision, allocatable :: neutral_local(:,:,:,:)

  integer :: ierr

  character(200) :: growth_rate_file
  character(200) :: growth_rate_max_file
  character(200) :: unstable_file
  character(200) :: neutral_file

  call start_mpi()
  ! write(*,*) my_rank
  call read_control_file()

  NM_local = NM / total_process
  log_M_min = log(M_min)
  log_M_max = log(M_max)

  log_N_min = log(N_min)
  log_N_max = log(N_max)

  allocate(growth_rate_local(1:NR, 1:NE, 1:NM_local))
  allocate(growth_rate_max_local(1:NR, 1:NE))
  allocate(unstable_local(1:NZ, 1:NR, 1:NE, 1:NM_local))
  allocate(neutral_local(1:NZ, 1:NR, 1:NE, 1:NM_local))

  if (my_rank == 0) then
    if (Output_max) then
      allocate(growth_rate_max(1:NR, 1:NE))
    else
      allocate(growth_rate_global(1:NR, 1:NE, 1:NM))
      allocate(unstable_global(1:NZ, 1:NR, 1:NE, 1:NM))
      allocate(neutral_global(1:NZ, 1:NR, 1:NE, 1:NM))
    endif
  endif

  do in = 1, NN
    if (N_log_scale) then
      log_N = log_N_min + (log_N_max - log_N_min) / (NN - 1) * (in - 1)
      buoyancy_frequency = exp(log_N)
    else
      buoyancy_frequency = N_min + (N_max - N_min) / (NN - 1) * (in - 1)
    endif

    do ir = 1, NR
      do ie = 1, NE
        do im_local = 1, NM_local
          im = NM_local * my_rank + im_local
          if (M_log_scale) then
            log_M = log_M_min + (log_M_max - log_M_min) / (NM - 1) * (im - 1)
            M = exp(log_M)
          else
            M = M_min + (M_max - M_min) / (NM - 1) * (im - 1)
          endif
          Ro = Ro_min + (Ro_max - Ro_min) / (NR - 1) * (ir - 1)
          e = e_min + (e_max - e_min) / (NE - 1) * (ie - 1)
          alpha_minus = Ro * (1 - e)**2 / (1 + (1 - e)**2)
          alpha_plus = Ro / (1 + (1 - e)**2)
          call elliptic_instability_solver(  &
            &  alpha_plus, alpha_minus, M,  &
            &  growth_rate, unstable_eigen_vector, neutral_eigen_vector)
          growth_rate_local(ir, ie, im_local) = growth_rate
          unstable_local(1:NZ, ir, ie, im_local) = unstable_eigen_vector(1:NZ)
          neutral_local(1:NZ, ir, ie, im_local) = neutral_eigen_vector(1:NZ)
        enddo
      end do
      if (my_rank == 0) write(*,*) in, ir
    enddo

    if (Output_max) then

      growth_rate_max_local(:, :) = maxval(growth_rate_local, 3)
      call MPI_Reduce(growth_rate_max_local(1,1),  &
      &  growth_rate_max(1,1), NR * NE, MPI_DOUBLE_PRECISION,  &
      &  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      if (my_rank == 0) then
        write(growth_rate_max_file, '(A, i6.6, A)')  &
        &  trim(data_directory)//'/growth_rate_max', in-1, '.out'
        open(UNIT_DATA, file=growth_rate_max_file, access='direct',  &
        &  status='unknown', form='unformatted', recl=8*NR*NE)
        write(UNIT_DATA, REC=1) growth_rate_max(:, :)
        close(UNIT_DATA)
      endif

    else

      call MPI_Gather(growth_rate_local(1,1,1), NR * NE * NM_local,  &
      &  MPI_DOUBLE_PRECISION, growth_rate_global(1,1,1),  &
      &  NR * NE * NM_local, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,  &
      &  ierr)

      call MPI_Gather(unstable_local(1,1,1,1), NZ * NR * NE * NM_local,  &
      &  MPI_DOUBLE_PRECISION, unstable_global(1,1,1,1),  &
      &  NZ * NR * NE * NM_local, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,  &
      &  ierr)

      call MPI_Gather(neutral_local(1,1,1,1), NZ * NR * NE * NM_local,  &
      &  MPI_DOUBLE_PRECISION, neutral_global(1,1,1,1),  &
      &  NZ * NR * NE * NM_local, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,  &
      &  ierr)

      if (my_rank == 0) then
        write(growth_rate_file, '(A, i6.6, A)')  &
        &  trim(data_directory)//'/growth_rate_R', in-1, '.out'
        do ir = 1, NR
          open(UNIT_DATA, file=growth_rate_file, access='direct',  &
          &  status='unknown', form='unformatted', recl=8*NM*NE)
          write(UNIT_DATA, REC=ir) growth_rate_global(ir, :, :)
          close(UNIT_DATA)
        enddo

        write(growth_rate_file, '(A, i6.6, A)')  &
        &  trim(data_directory)//'/growth_rate_E', in-1, '.out'
        do ie = 1, NE
          open(UNIT_DATA, file=growth_rate_file, access='direct',  &
          &  status='unknown', form='unformatted', recl=8*NM*NR)
          write(UNIT_DATA, REC=ie) growth_rate_global(:, ie, :)
          close(UNIT_DATA)
        enddo

        write(unstable_file, '(A, i6.6, A)')  &
        &  trim(data_directory)//'/unstable_R', in-1, '.out'
        do ir = 1, NR
          open(UNIT_DATA, file=unstable_file, access='direct',  &
          &  status='unknown', form='unformatted', recl=8*NZ*NM*NE)
          write(UNIT_DATA, REC=ir) unstable_global(:, ir, :, :)
          close(UNIT_DATA)
        enddo

        write(unstable_file, '(A, i6.6, A)')  &
        &  trim(data_directory)//'/unstable_E', in-1, '.out'
        do ie = 1, NE
          open(UNIT_DATA, file=unstable_file, access='direct',  &
          &  status='unknown', form='unformatted', recl=8*NZ*NM*NR)
          write(UNIT_DATA, REC=ie) unstable_global(:, :, ie, :)
          close(UNIT_DATA)
        enddo

        write(neutral_file, '(A, i6.6, A)')  &
        &  trim(data_directory)//'/neutral_R', in-1, '.out'
        do ir = 1, NR
          open(UNIT_DATA, file=neutral_file, access='direct',  &
          &  status='unknown', form='unformatted', recl=8*NZ*NM*NE)
          write(UNIT_DATA, REC=ir) neutral_global(:, ir, :, :)
          close(UNIT_DATA)
        enddo

        write(neutral_file, '(A, i6.6, A)')  &
        &  trim(data_directory)//'/neutral_E', in-1, '.out'
        do ie = 1, NE
          open(UNIT_DATA, file=neutral_file, access='direct',  &
          &  status='unknown', form='unformatted', recl=8*NZ*NM*NR)
          write(UNIT_DATA, REC=ie) neutral_global(:, :, ie, :)
          close(UNIT_DATA)
        enddo
      end if

    endif
  end do

  call end_mpi()

  stop
end program main

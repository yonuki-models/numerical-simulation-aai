module Input_Output
  use mpi, only : MPI_COMM_WORLD
  use mpi, only : MPI_REAL
  use mpi, only : MPI_COMPLEX
  use mpi, only : MPI_SUM
  use mpi, only : MPI_STATUS_SIZE
  use Parallel, only : my_rank
  use Parallel, only : total_process
  use Parallel, only : comm_1
  use Parallel, only : comm_2
  use Parallel, only : rank_1
  use Parallel, only : rank_2
  use Parallel, only : size_1
  use Parallel, only : size_2
  use Time_Module, only : current_time
  use Time_Module, only : file_number
  use Time_Module, only : output_state_interval
  use Geometry, only : n1C_local, n2C_local, n3
  use Geometry, only : n1, n2R_local, n3R_local
  use Geometry, only : n2
  use Geometry, only : n1_truncated, n2_truncated, n3_truncated
  use Geometry, only : state
  use Geometry, only : energy
  use Geometry, only : spectrum_H
  use Geometry, only : spectrum_HV
  implicit none
  private

  integer, parameter :: LEN_FILE_NAME = 64
  integer, parameter :: LEN_DIRECTORY = 256
  integer, parameter :: LEN_FILE_PATH = 256
  integer, public, parameter :: LEN_MESSAGE = 1024

  integer, public, parameter :: UNIT_STDERR = 0
  integer, public, parameter :: UNIT_STDIN = 5
  integer, public, parameter :: UNIT_STDOUT = 6
  integer, public, parameter :: UNIT_FILE_DEFAULT = 10

  integer :: record_factor = 4
             ! depends on compiler
             ! = 4 for gfortran or ifort with option "-assume byterecl"
             ! = 1 for ifort without option
  integer, public, parameter :: FACTOR_SINGLE_PRECISION = 1
  integer, public, parameter :: FACTOR_DOUBLE_PRECISION = 2
  integer, public, parameter :: FACTOR_SINGLE_COMPLEX = 2
  integer, public, parameter :: FACTOR_DOUBLE_COMPLEX = 4

  double precision, save :: reference_energy = 0.d0
  double precision, allocatable :: HV_low_TKE_prev(:,:)
  double precision, allocatable :: HV_low_TPE_prev(:,:)

  character(len=LEN_DIRECTORY), public ::  &
    &  control_directory = '../control'
  character(len=LEN_DIRECTORY), public ::  &
    &  data_directory = '../data/exDefault'

  public print_main
  public output_geometry
  public initialize_reference_energy
  public output_energy_sum
  public output_state
  public output_energy
  public output_snapshot_all
  public input_state
  public output_spectrum
  public output_Richardson_distribution
  public output_Froude_spectrum
  public output_Thorpe_distribution

contains

! *****************************************************************************

  subroutine print_main(message, value)
    implicit none
    character(*), intent(in) :: message
    double precision, intent(in), optional :: value

    if (my_rank == 0) then
      write(unit_stdout, *) message
      if (present(value)) write(unit_stdout, *) value
    end if

    return
  end subroutine print_main

! *****************************************************************************

  subroutine output_geometry()
    use Geometry, only : domain_length
    use External_Field, only : derive_transform_tensor
    implicit none
    double precision :: transform_tensor(1:3, 1:3)
    integer :: number_of_column = 16
    character(len=LEN_FILE_NAME) :: file_name
    character(len=LEN_FILE_PATH) :: file_path

    integer :: access

    call derive_transform_tensor(current_time, transform_tensor)

    if (my_rank == 0) then

      ! Output CSV file
      file_name = "geometry.csv"
      file_path = trim(data_directory)//'/'//trim(file_name)

      if (file_number==0 .or. access(file_path, " ")/=0) then
        open(unit=UNIT_FILE_DEFAULT, file=file_path,  &
          &  status='replace', form='formatted')
        write(UNIT_FILE_DEFAULT, '(A)')  &
          &  'number,time,n1,n2,n3,l1,l2,l3,&
          &A11,A12,A13,A21,A22,A23,A31,A32,A33'
      else
        open(unit=UNIT_FILE_DEFAULT, file=file_path,  &
          &  status='old', form='formatted', position='append')
      endif

      write(UNIT_FILE_DEFAULT, '(I0, ",", e16.8, ",", I0, ",",  &
        &I0, ",", I0 ",", e16.8, ",", e16.8, ",", e16.8, ",",  &
        &e16.8, ",", e16.8, ",", e16.8, ",", e16.8, ",", e16.8, ",",  &
        &e16.8, ",", e16.8, ",", e16.8, ",", e16.8)')  &
        &  file_number,            &
        &  current_time,           &
        &  n1,n2,n3,               &
        &  domain_length(1),       &
        &  domain_length(2),       &
        &  domain_length(3),       &
        &  transform_tensor(1,1),  &
        &  transform_tensor(1,2),  &
        &  transform_tensor(1,3),  &
        &  transform_tensor(2,1),  &
        &  transform_tensor(2,2),  &
        &  transform_tensor(2,3),  &
        &  transform_tensor(3,1),  &
        &  transform_tensor(3,2),  &
        &  transform_tensor(3,3)

      close(UNIT_FILE_DEFAULT)


      ! Output binary file
      file_name = "geometry.out"
      file_path = trim(data_directory)//'/'//trim(file_name)

      open(unit=UNIT_FILE_DEFAULT, file=file_path,  &
        &  status='unknown', form='unformatted', access='direct',  &
        &  recl=record_factor * number_of_column)

      write(unit=UNIT_FILE_DEFAULT, rec=file_number+1)  &
        &  real(current_time),           &
        &  real(n1),                     &
        &  real(n2),                     &
        &  real(n3),                     &
        &  real(domain_length(1)),       &
        &  real(domain_length(2)),       &
        &  real(domain_length(3)),       &
        &  real(transform_tensor(1,1)),  &
        &  real(transform_tensor(1,2)),  &
        &  real(transform_tensor(1,3)),  &
        &  real(transform_tensor(2,1)),  &
        &  real(transform_tensor(2,2)),  &
        &  real(transform_tensor(2,3)),  &
        &  real(transform_tensor(3,1)),  &
        &  real(transform_tensor(3,2)),  &
        &  real(transform_tensor(3,3))

      close(UNIT_FILE_DEFAULT)

    endif

    return
  end subroutine output_geometry

! *****************************************************************************

  subroutine initialize_reference_energy()
    use Analysis, only : sum_real_data
    implicit none
    double precision :: sum_KE
    double precision :: sum_PE
    double precision :: kinetic_energy(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: potential_energy(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    call state%calculate_kinetic_energy(kinetic_energy)
    call state%calculate_potential_energy(potential_energy)

    sum_KE = sum_real_data(kinetic_energy)
    sum_PE = sum_real_data(potential_energy)
    reference_energy = sum_KE + sum_PE

    return
  end subroutine initialize_reference_energy

! *****************************************************************************

  subroutine output_energy_sum()
    use Analysis, only : derive_wave_vortex_energy
    use Analysis, only : calculate_PV
    use Analysis, only : sum_real_data
    use Analysis, only : get_primary_energy
    use Dynamics, only : decay_coefficient_H
    use Dynamics, only : decay_coefficient_V
    implicit none
    double precision :: sum_KE
    double precision :: sum_PE
    double precision :: total_energy
    double precision :: wave_energy
    double precision :: primary_energy
    double precision :: vortex_energy
    double precision :: sum_KE_transfer
    double precision :: sum_PE_transfer
    double precision :: sum_KE_production
    double precision :: sum_PE_production
    double precision :: sum_WE_production_1
    double precision :: sum_WE_production_2
    double precision :: sum_VE_production
    double precision :: sum_conversion
    double precision :: sum_KE_dissipation
    double precision :: sum_PE_dissipation
    double precision :: enstrophy2
    double precision :: enstrophy3
    double precision :: enstrophy4
    double precision :: wave_spectrum(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: vortex_spectrum(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: kinetic_energy(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: potential_energy(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: Q1(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex(kind(0d0)) :: Q2(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    integer :: access
    integer :: number_of_column = 22
    character(len=LEN_FILE_NAME) :: file_name
    character(len=LEN_FILE_PATH) :: file_path

    double precision :: energy_budget

    call state%calculate_kinetic_energy(kinetic_energy)
    call state%calculate_potential_energy(potential_energy)
    call derive_wave_vortex_energy(state%velocity, state%temperature,  &
      &  current_time, wave_spectrum, vortex_spectrum)
    call calculate_PV(state%velocity, state%temperature, current_time, Q1, Q2)

    primary_energy = get_primary_energy(wave_spectrum)
    wave_energy = sum_real_data(wave_spectrum)
    vortex_energy = sum_real_data(vortex_spectrum)
    enstrophy2 = sum_real_data(abs(Q1)**2)
    enstrophy3 = 2 * sum_real_data(dble(conjg(Q1) * Q2))
    enstrophy4 = sum_real_data(abs(Q2)**2)

    sum_KE = sum_real_data(kinetic_energy)
    sum_PE = sum_real_data(potential_energy)
    sum_KE_production = sum(spectrum_HV%kinetic_energy_production)
    sum_PE_production = sum(spectrum_HV%potential_energy_production)
    sum_WE_production_1 = sum(spectrum_HV%wave_production_1)
    sum_WE_production_2 = sum(spectrum_HV%wave_production_2)
    sum_VE_production = sum(spectrum_HV%vortex_production)
    sum_KE_dissipation = sum(spectrum_HV%kinetic_energy_dissipation)
    sum_PE_dissipation = sum(spectrum_HV%potential_energy_dissipation)
    sum_conversion = sum(spectrum_HV%energy_conversion)
    sum_KE_transfer = sum(spectrum_HV%kinetic_energy_transfer)
    sum_PE_transfer = sum(spectrum_HV%potential_energy_transfer)

    total_energy = sum_KE + sum_PE

    if (my_rank == 0) then
      ! Output CSV file
      file_name = "budget.csv"
      file_path = trim(data_directory)//'/'//trim(file_name)

      if (file_number==0 .or. access(file_path, " ")/=0) then
        open(unit=UNIT_FILE_DEFAULT, file=file_path,  &
          &  status='replace', form='formatted')
        write(UNIT_FILE_DEFAULT, '(A)') 'T,KE,PE,E,WE,VE,PWE,&
          &TKE,TPE,PKE,PPE,PW1,PW2,PVE,CKP,DKE,DPE,EN2,EN3,EN4,DCH,DCV'
      else
        open(unit=UNIT_FILE_DEFAULT, file=file_path,  &
          &  status='old', form='formatted', position='append')
      endif

      write(UNIT_FILE_DEFAULT, '(21(e16.8, ","), e16.8)')  &
        &  current_time,          &
        &  sum_KE,                &
        &  sum_PE,                &
        &  total_energy,          &
        &  wave_energy,           &
        &  vortex_energy,         &
        &  primary_energy,        &
        &  sum_KE_transfer,       &
        &  sum_PE_transfer,       &
        &  sum_KE_production,     &
        &  sum_PE_production,     &
        &  sum_WE_production_1,   &
        &  sum_WE_production_2,   &
        &  sum_VE_production,     &
        &  sum_conversion,        &
        &  sum_KE_dissipation,    &
        &  sum_PE_dissipation,    &
        &  enstrophy2,            &
        &  enstrophy3,            &
        &  enstrophy4,            &
        &  decay_coefficient_H,   &
        &  decay_coefficient_V

      close(UNIT_FILE_DEFAULT)

      energy_budget = (sum_KE + sum_PE - reference_energy  &
        &  - sum_KE_transfer - sum_PE_transfer  &
        &  - sum_KE_production - sum_PE_production  &
        &  + sum_KE_dissipation + sum_PE_dissipation)  &
        &  / (abs(sum_KE + sum_PE - reference_energy)  &
        &  + abs(sum_KE_transfer + sum_PE_transfer)  &
        &  + abs(sum_KE_production + sum_PE_production)  &
        &  + abs(sum_KE_dissipation + sum_PE_dissipation))

      call print_main("Precision of energy budget: ", energy_budget)
      reference_energy = sum_KE + sum_PE

      ! Output binary file
      file_name = "budget.out"
      file_path = trim(data_directory)//'/'//trim(file_name)

      open(unit=UNIT_FILE_DEFAULT, file=file_path,  &
        &  status='unknown', form='unformatted', access='direct',  &
        &  recl=record_factor * number_of_column)

      write(unit=UNIT_FILE_DEFAULT, rec=file_number+1)  &
        &  real(current_time),          &
        &  real(sum_KE),                &
        &  real(sum_PE),                &
        &  real(total_energy),          &
        &  real(wave_energy),           &
        &  real(vortex_energy),         &
        &  real(primary_energy),        &
        &  real(sum_KE_transfer),       &
        &  real(sum_PE_transfer),       &
        &  real(sum_KE_production),     &
        &  real(sum_PE_production),     &
        &  real(sum_WE_production_1),   &
        &  real(sum_WE_production_2),   &
        &  real(sum_VE_production),     &
        &  real(sum_conversion),        &
        &  real(sum_KE_dissipation),    &
        &  real(sum_PE_dissipation),    &
        &  real(enstrophy2),            &
        &  real(enstrophy3),            &
        &  real(enstrophy4),            &
        &  real(decay_coefficient_H),   &
        &  real(decay_coefficient_V)

      close(UNIT_FILE_DEFAULT)

    endif

    return
  end subroutine output_energy_sum

! *****************************************************************************

  subroutine output_state()
    implicit none

    call output_complex_data('U')
    call output_complex_data('V')
    call output_complex_data('W')
    call output_complex_data('T')

    return
  end subroutine output_state

! *****************************************************************************

  subroutine output_energy()
    implicit none

    call output_real_data('NKE')
    call output_real_data('NPE')
    call output_real_data('PKE')
    call output_real_data('PPE')
    call output_real_data('CKP')
    call output_real_data('DKE')
    call output_real_data('DPE')

    return
  end subroutine output_energy

! *****************************************************************************

  subroutine output_complex_data(variable_name)
    implicit none
    character(LEN=*), intent(in) :: variable_name

    complex, allocatable :: C_Gat_trans(:,:,:)  ! Single precision
    complex, allocatable :: C_Gat(:,:,:)
    complex, allocatable :: C_Out(:,:,:)
    complex :: C_tmp(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex :: C_trans(0:n2C_local-1, 0:n3-1, 0:n1C_local-1)
    integer :: record_length
    integer :: nt1  ! short name of n1_truncated
    integer :: nt2  ! short name of n2_truncated
    integer :: nt3  ! short name of n3_truncated
    integer :: n1_ext

    character(LEN=LEN_FILE_PATH) :: file_path
    integer :: i1
    integer :: ierr

    integer :: send_size
    integer :: record_number
    logical :: IO_process

    nt1 = n1_truncated
    nt2 = n2_truncated
    nt3 = n3_truncated
    n1_ext = n1C_local * size_1

    IO_process = (rank_1 == 0)
    if (IO_process) then
      allocate(C_Gat_trans(0:n2C_local-1, 0:n3-1, 0:n1_ext-1),  &
      &  source=(0.0, 0.0))
      allocate(C_Gat(0:n1_ext-1, 0:n2C_local-1, 0:n3-1), source=(0.0, 0.0))
      allocate(C_Out(0:nt1, 0:n2C_local-1, 0:2*nt3), source=(0.0, 0.0))
    end if

    if (variable_name == 'U') C_tmp(:,:,:) = cmplx(state%velocity(1,:,:,:))
    if (variable_name == 'V') C_tmp(:,:,:) = cmplx(state%velocity(2,:,:,:))
    if (variable_name == 'W') C_tmp(:,:,:) = cmplx(state%velocity(3,:,:,:))
    if (variable_name == 'T') C_tmp(:,:,:) = cmplx(state%temperature(:,:,:))

    do i1 = 0, n1C_local - 1
      C_trans(:, :, i1) = C_tmp(i1, :, :)
    enddo

    send_size = n1C_local * n2C_local * n3
    call MPI_Gather(C_trans(0,0,0), send_size, MPI_COMPLEX,  &
    &  C_Gat_trans(0,0,0), send_size, MPI_COMPLEX, 0, comm_1, ierr)

    if (IO_process) then

      if (rank_2 * n2C_local > nt2  &
      &  .and. (rank_2 + 1) * n2C_local - 1 < (n2 - nt2)) then
        deallocate(C_Gat_trans, C_Gat, C_Out)
        return
      endif

      do i1 = 0, n1_ext - 1
        C_Gat(i1, :, :) = C_Gat_trans(:, :, i1)
      enddo

      C_Out(0:nt1, 0:n2C_local-1, 0:nt3)  &
        &  = C_Gat(0:nt1, 0:n2C_local-1, 0:nt3)
      C_Out(0:nt1, 0:n2C_local-1, nt3+1:2*nt3)  &
        &  = C_Gat(0:nt1, 0:n2C_local-1, n3-nt3:n3-1)

      write(file_path, '(A, "/state/", A, i4.4, ".out")')  &
        &  trim(data_directory), trim(variable_name), rank_2

      record_length = (nt1 + 1) * n2C_local * (2 * nt3 + 1)
      record_number = file_number / output_state_interval + 1

      open(unit=UNIT_FILE_DEFAULT, file=file_path, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * record_length * FACTOR_SINGLE_COMPLEX)
      write(UNIT_FILE_DEFAULT, rec=record_number) C_Out
      close(UNIT_FILE_DEFAULT)

      deallocate(C_Gat_trans, C_Gat, C_Out)
    end if

    return
  end subroutine output_complex_data

! *****************************************************************************

  subroutine output_real_data(variable_name)
    implicit none
    character(LEN=*), intent(in) :: variable_name

    real, allocatable :: R_Gat_trans(:,:,:)  ! Single precision
    real, allocatable :: R_Gat(:,:,:)
    real, allocatable :: R_Out(:,:,:)
    real :: R_tmp(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    real :: R_trans(0:n2C_local-1, 0:n3-1, 0:n1C_local-1)
    integer :: record_length
    integer :: nt1  ! short name of n1_truncated
    integer :: nt2  ! short name of n2_truncated
    integer :: nt3  ! short name of n3_truncated
    integer :: n1_ext

    character(LEN=LEN_FILE_PATH) :: file_path
    integer :: i1
    integer :: ierr

    integer :: send_size
    integer :: record_number
    logical :: IO_process

    nt1 = n1_truncated
    nt2 = n2_truncated
    nt3 = n3_truncated
    n1_ext = n1C_local * size_1

    IO_process = (rank_1 == 0)
    if (IO_process) then
      allocate(R_Gat_trans(0:n2C_local-1, 0:n3-1, 0:n1_ext-1),  &
        &  source=0.0)
      allocate(R_Gat(0:n1_ext-1, 0:n2C_local-1, 0:n3-1), source=0.0)
      allocate(R_Out(0:nt1, 0:n2C_local-1, 0:2*nt3), source=0.0)
    end if

    if (variable_name == 'NKE')  &
      &  R_tmp(:,:,:) = real(energy%kinetic_energy_transfer(:,:,:))
    if (variable_name == 'NPE')  &
      &  R_tmp(:,:,:) = real(energy%potential_energy_transfer(:,:,:))
    if (variable_name == 'PKE')  &
      &  R_tmp(:,:,:) = real(energy%kinetic_energy_production(:,:,:))
    if (variable_name == 'PPE')  &
      &  R_tmp(:,:,:) = real(energy%potential_energy_production(:,:,:))
    if (variable_name == 'CKP')  &
      &  R_tmp(:,:,:) = real(energy%energy_conversion(:,:,:))
    if (variable_name == 'DKE')  &
      &  R_tmp(:,:,:) = real(energy%kinetic_energy_dissipation(:,:,:))
    if (variable_name == 'DPE')  &
      &  R_tmp(:,:,:) = real(energy%potential_energy_dissipation(:,:,:))

    do i1 = 0, n1C_local - 1
      R_trans(:, :, i1) = R_tmp(i1, :, :)
    enddo

    send_size = n1C_local * n2C_local * n3
    call MPI_Gather(R_trans(0,0,0), send_size, MPI_REAL,  &
      &  R_Gat_trans(0,0,0), send_size, MPI_REAL, 0, comm_1, ierr)

    if (IO_process) then

      if (rank_2 * n2C_local > nt2  &
        &  .and. (rank_2 + 1) * n2C_local - 1 < (n2 - nt2)) then
        deallocate(R_Gat_trans, R_Gat, R_Out)
        return
      endif

      do i1 = 0, n1_ext - 1
        R_Gat(i1, :, :) = R_Gat_trans(:, :, i1)
      enddo

      R_Out(0:nt1, 0:n2C_local-1, 0:nt3)  &
        &  = R_Gat(0:nt1, 0:n2C_local-1, 0:nt3)
      R_Out(0:nt1, 0:n2C_local-1, nt3+1:2*nt3)  &
        &  = R_Gat(0:nt1, 0:n2C_local-1, n3-nt3:n3-1)

      write(file_path, '(A, "/energy/", A, i4.4, ".out")')  &
        &  trim(data_directory), trim(variable_name), rank_2

      record_length = (nt1 + 1) * n2C_local * (2 * nt3 + 1)
      record_number = file_number / output_state_interval + 1

      open(unit=UNIT_FILE_DEFAULT, file=file_path, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * record_length * FACTOR_SINGLE_PRECISION)
      write(UNIT_FILE_DEFAULT, rec=record_number) R_Out
      close(UNIT_FILE_DEFAULT)

      deallocate(R_Gat_trans, R_Gat, R_Out)
    end if

    return
  end subroutine output_real_data

! *****************************************************************************

  subroutine output_snapshot_all()
    implicit none

    call output_snapshot('U')
    call output_snapshot('V')
    call output_snapshot('W')
    call output_snapshot('T')

    return
  end subroutine output_snapshot_all

! *****************************************************************************

  subroutine output_snapshot(variable_name)
    use mpi, only: MPI_REAL
    use mpi, only: MPI_STATUS_SIZE
    use FFTW_Interface, only : fftw_backward_parallel
    use Geometry, only : n1R_local
    implicit none
    character(LEN=*), intent(in) :: variable_name

    complex(kind(0d0)) :: C_data(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: R_data(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
    real :: data_XY_local(0:n1-1, 0:n2R_local-1)
    real :: data_XZ_local(0:n1-1, 0:n3R_local-1)
    real, allocatable :: data_XY(:,:)
    real, allocatable :: data_XZ_send(:,:)
    real, allocatable :: data_XZ(:,:)
    integer :: record_length

    integer :: recv_rank
    integer :: send_rank
    integer :: ista(MPI_STATUS_SIZE)

    character(LEN=LEN_FILE_PATH) :: file_path
    integer :: ierr

    if (rank_1 == 0) allocate(data_XY(0:n1-1, 0:n2-1))

    if (variable_name == 'U') C_data(:,:,:) = cmplx(state%velocity(1,:,:,:))
    if (variable_name == 'V') C_data(:,:,:) = cmplx(state%velocity(2,:,:,:))
    if (variable_name == 'W') C_data(:,:,:) = cmplx(state%velocity(3,:,:,:))
    if (variable_name == 'T') C_data(:,:,:) = cmplx(state%temperature(:,:,:))

    call fftw_backward_parallel(C_data, R_data)

    data_XY_local(:,:) = real(R_data(:,:,0))
    call MPI_Gather(data_XY_local(0,0), n1 * n2R_local,  &
      &  MPI_REAL, data_XY(0,0), n1 * n2R_local,  &
      &  MPI_REAL, 0, comm_1, ierr)

    if (rank_1 == 0) then

      record_length = n1 * n2
      write(file_path, '(A, "/snapshot/XY_", A, i4.4, ".out")')  &
        &  trim(data_directory), trim(variable_name), rank_2
      open(unit=UNIT_FILE_DEFAULT, file=file_path, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * record_length * FACTOR_SINGLE_PRECISION)
      write(UNIT_FILE_DEFAULT, rec=file_number+1) data_XY
      close(UNIT_FILE_DEFAULT)

      deallocate(data_XY)

    endif

    if (size_1 /= size_2) return

    if (rank_1 == 0) allocate(data_XZ(0:n1-1, 0:n3-1))
    if (rank_2 == 0) allocate(data_XZ_send(0:n1-1, 0:n3-1))

    data_XZ_local(:,:) = real(R_data(:,0,:))
    call MPI_Gather(data_XZ_local(0,0), n1 * n3R_local, MPI_REAL,  &
    & data_XZ_send(0,0), n1 * n3R_local, MPI_REAL, 0, comm_2, ierr)

    if (rank_1 == 0 .and. rank_2 == 0) then
      data_XZ(:,:) = data_XZ_send(:,:)
    elseif (rank_2 == 0) then
      recv_rank = rank_1 * size_2
      call MPI_Send(data_XZ_send(0,0), n1 * n3, MPI_REAL,  &
      &  recv_rank, rank_1, MPI_COMM_WORLD, ierr)
    elseif (rank_1 == 0) then
      send_rank = rank_2
      call MPI_Recv(data_XZ(0,0), n1 * n3, MPI_REAL,  &
      &  send_rank, rank_2, MPI_COMM_WORLD, ista, ierr)
    endif

    if (rank_1 == 0) then

      record_length = n1 * n3
      write(file_path, '(A, "/snapshot/XZ_", A, i4.4, ".out")')  &
        &  trim(data_directory), trim(variable_name), rank_2
      open(unit=UNIT_FILE_DEFAULT, file=file_path, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * record_length * FACTOR_SINGLE_PRECISION)
      write(UNIT_FILE_DEFAULT, rec=file_number+1) data_XZ
      close(UNIT_FILE_DEFAULT)

    endif

    if (rank_1 == 0) deallocate(data_XZ)
    if (rank_2 == 0) deallocate(data_XZ_send)

    return
  end subroutine output_snapshot

! *****************************************************************************

  subroutine input_state()
    implicit none

    call input_complex_data('U')
    call input_complex_data('V')
    call input_complex_data('W')
    call input_complex_data('T')

    return
  end subroutine input_state

! *****************************************************************************

  subroutine input_complex_Data(variable_name)
    implicit none
    character(LEN=*), intent(in) :: variable_name

    complex, allocatable :: C_Gat_trans(:,:,:)  ! Single precision
    complex, allocatable :: C_Gat(:,:,:)  ! Single precision
    complex, allocatable :: C_In(:,:,:)
    complex :: C_tmp(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    complex :: C_trans(0:n2C_local-1, 0:n3-1, 0:n1C_local-1)
    integer :: record_length
    integer :: nt1  ! short name of n1_truncated
    integer :: nt2  ! short name of n2_truncated
    integer :: nt3  ! short name of n3_truncated
    integer :: n1_ext

    character(LEN=LEN_FILE_PATH) :: file_path
    integer :: i1
    integer :: ierr

    integer :: send_size
    integer :: record_number
    logical :: IO_process

    nt1 = n1_truncated
    nt2 = n2_truncated
    nt3 = n3_truncated
    n1_ext = n1C_local * size_1

    IO_process = (rank_1 == 0)
    if (IO_process) then
      allocate(C_Gat_trans(0:n2C_local-1, 0:n3-1, 0:n1_ext-1),  &
        &  source=(0.0, 0.0))
      allocate(C_Gat(0:n1_ext-1, 0:n2C_local-1, 0:n3-1), source=(0.0, 0.0))
      allocate(C_In(0:nt1, 0:n2C_local-1, 0:2*nt3), source=(0.0, 0.0))
    end if

    if (IO_process) then

      if (.not. (rank_2 * n2C_local > nt2  &
      &  .and. (rank_2 + 1) * n2C_local - 1 < (n2 - nt2))) then

        write(file_path, '(A, "/state/", A, i4.4, ".out")')  &
        &  trim(data_directory), trim(variable_name), rank_2

        record_length = (nt1 + 1) * n2C_local * (2 * nt3 + 1)
        record_number = file_number / output_state_interval + 1

        open(unit=UNIT_FILE_DEFAULT, file=file_path, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * record_length * FACTOR_SINGLE_COMPLEX)
        read(UNIT_FILE_DEFAULT, rec=record_number) C_In
        close(UNIT_FILE_DEFAULT)

        C_Gat(0:nt1, 0:n2C_local-1, 0:nt3)  &
        &  = C_In(0:nt1, 0:n2C_local-1, 0:nt3)
        C_Gat(0:nt1, 0:n2C_local-1, n3-nt3:n3-1)  &
        &  = C_In(0:nt1, 0:n2C_local-1, nt3+1:2*nt3)

      endif

      do i1 = 0, n1_ext - 1
        C_Gat_trans(:, :, i1) = C_Gat(i1, :, :)
      enddo

    end if

    send_size = n1C_local * n2C_local * n3
    call MPI_Scatter(C_Gat_trans(0,0,0), send_size, MPI_COMPLEX,  &
      &  C_trans(0,0,0), send_size, MPI_COMPLEX, 0, comm_1, ierr)

    do i1 = 0, n1C_local - 1
      C_tmp(i1, :, :) = C_trans(:, :, i1)
    enddo

    if (variable_name == 'U') state%velocity(1,:,:,:) = C_tmp(:,:,:)
    if (variable_name == 'V') state%velocity(2,:,:,:) = C_tmp(:,:,:)
    if (variable_name == 'W') state%velocity(3,:,:,:) = C_tmp(:,:,:)
    if (variable_name == 'T') state%temperature(:,:,:) = C_tmp(:,:,:)

    if (IO_process) deallocate(C_Gat_trans, C_Gat, C_In)

    return
  end subroutine input_complex_Data

! *****************************************************************************

  subroutine output_spectrum()
    implicit none

    call output_spectrum_HV('KE')
    call output_spectrum_HV('PE')
    call output_spectrum_HV('NKE')
    call output_spectrum_HV('NPE')
    call output_spectrum_HV('PKE')
    call output_spectrum_HV('PPE')
    call output_spectrum_HV('CKP')
    call output_spectrum_HV('DKE')
    call output_spectrum_HV('DPE')
    call output_spectrum_HV('VE')
    call output_spectrum_HV('WE')
    call output_spectrum_HV('EN2')
    call output_spectrum_HV('EN3')
    call output_spectrum_HV('EN4')
    call output_spectrum_HV('KFC')
    call output_spectrum_HV('PFC')

    call output_spectrum_horizontal('KE')
    call output_spectrum_horizontal('PE')
    call output_spectrum_horizontal('NKE')
    call output_spectrum_horizontal('NPE')
    call output_spectrum_horizontal('PKE')
    call output_spectrum_horizontal('PPE')
    call output_spectrum_horizontal('CKP')
    call output_spectrum_horizontal('DKE')
    call output_spectrum_horizontal('DPE')

    return
  end subroutine output_spectrum

! *****************************************************************************

  subroutine output_spectrum_HV(variable_name)
    use Geometry, only : calculate_wavenumbers
    use Geometry, only : NH_max
    use Geometry, only : NV_max
    use Analysis, only : spectrum_integrate_azimuth
    use Analysis, only : derive_wave_vortex_energy
    use Analysis, only : calculate_PV
    use External_Field, only : derive_transform_tensor
    implicit none
    character(LEN=*), intent(in) :: variable_name

    double precision :: A(1:3, 1:3)
    double precision :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

    double precision, allocatable :: tmp(:,:,:)
    complex(kind(0d0)), allocatable :: tmpC(:,:,:)
    complex(kind(0d0)), allocatable :: tmpC2(:,:,:)
    double precision :: data_(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: spectrum_(0:NH_max, 0:NV_max)
    integer :: i

    integer :: record_length
    character(LEN=LEN_FILE_PATH) :: file_path

    record_length = (NH_max + 1) * (NV_max + 1)
    call derive_transform_tensor(current_time, A)
    call calculate_wavenumbers(A, K)

    if (variable_name == 'KE') then
      data_(:,:,:) = 0.d0
      do i = 1, 3
        data_(:,:,:) = data_(:,:,:) + abs(state%velocity(i,:,:,:))**2
      enddo
      data_(:,:,:) = data_(:,:,:) * 0.5d0
      call spectrum_integrate_azimuth(data_, K, spectrum_)
    endif

    if (variable_name == 'PE') then
      data_(:,:,:) = abs(state%temperature(:,:,:))**2 * 0.5d0
      call spectrum_integrate_azimuth(data_, K, spectrum_)
    endif

    if (variable_name == 'WE') then
      allocate(tmp(0:n1C_local-1, 0:n2C_local-1, 0:n3-1))
      call derive_wave_vortex_energy(state%velocity, state%temperature,  &
        &  current_time, data_, tmp)
      call spectrum_integrate_azimuth(data_, K, spectrum_)
    endif

    if (variable_name == 'VE') then
      allocate(tmp(0:n1C_local-1, 0:n2C_local-1, 0:n3-1))
      call derive_wave_vortex_energy(state%velocity, state%temperature,  &
        &  current_time, tmp, data_)
      call spectrum_integrate_azimuth(data_, K, spectrum_)
    endif

    if (variable_name == 'EN2') then
      allocate(tmpC(0:n1C_local-1, 0:n2C_local-1, 0:n3-1))
      allocate(tmpC2(0:n1C_local-1, 0:n2C_local-1, 0:n3-1))
      call calculate_PV(state%velocity, state%temperature,  &
        &  current_time, tmpC, tmpC2)
      data_(:,:,:) = abs(tmpC)**2
      call spectrum_integrate_azimuth(data_, K, spectrum_)
    endif

    if (variable_name == 'EN3') then
      allocate(tmpC(0:n1C_local-1, 0:n2C_local-1, 0:n3-1))
      allocate(tmpC2(0:n1C_local-1, 0:n2C_local-1, 0:n3-1))
      call calculate_PV(state%velocity, state%temperature,  &
        &  current_time, tmpC, tmpC2)
      data_(:,:,:) = 2 * dble(conjg(tmpC) * tmpC2)
      call spectrum_integrate_azimuth(data_, K, spectrum_)
    endif

    if (variable_name == 'EN4') then
      allocate(tmpC(0:n1C_local-1, 0:n2C_local-1, 0:n3-1))
      allocate(tmpC2(0:n1C_local-1, 0:n2C_local-1, 0:n3-1))
      call calculate_PV(state%velocity, state%temperature,  &
        &  current_time, tmpC, tmpC2)
      data_(:,:,:) = abs(tmpC2)**2
      call spectrum_integrate_azimuth(data_, K, spectrum_)
    endif

    if (variable_name == 'NKE')  &
      &  spectrum_(:,:) = spectrum_HV%kinetic_energy_transfer(:,:)
    if (variable_name == 'NPE')  &
      &  spectrum_(:,:) = spectrum_HV%potential_energy_transfer(:,:)
    if (variable_name == 'PKE')  &
      &  spectrum_(:,:) = spectrum_HV%kinetic_energy_production(:,:)
    if (variable_name == 'PPE')  &
      &  spectrum_(:,:) = spectrum_HV%potential_energy_production(:,:)
    if (variable_name == 'PW1')  &
      &  spectrum_(:,:) = spectrum_HV%wave_production_1(:,:)
    if (variable_name == 'PW2')  &
      &  spectrum_(:,:) = spectrum_HV%wave_production_2(:,:)
    if (variable_name == 'PVE')  &
      &  spectrum_(:,:) = spectrum_HV%vortex_production(:,:)
    if (variable_name == 'CKP')  &
      &  spectrum_(:,:) = spectrum_HV%energy_conversion(:,:)
    if (variable_name == 'DKE')  &
      &  spectrum_(:,:) = spectrum_HV%kinetic_energy_dissipation(:,:)
    if (variable_name == 'DPE')  &
      &  spectrum_(:,:) = spectrum_HV%potential_energy_dissipation(:,:)
    if (variable_name == 'KFC')  &
      &  spectrum_(:,:) = spectrum_HV%kinetic_flux_convergence(:,:)
    if (variable_name == 'PFC')  &
      &  spectrum_(:,:) = spectrum_HV%potential_flux_convergence(:,:)

    if (my_rank == 0) then

      write(file_path, '(A, "/Spec_HV_", A, ".out")')  &
        &  trim(data_directory), trim(variable_name)

      open(unit=UNIT_FILE_DEFAULT, file=file_path, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * record_length * FACTOR_SINGLE_PRECISION)

      write(UNIT_FILE_DEFAULT, rec=file_number+1) real(spectrum_)

      close(UNIT_FILE_DEFAULT)

    endif

    return
  end subroutine output_spectrum_HV

! *****************************************************************************

  subroutine output_spectrum_horizontal(variable_name)
    use Geometry, only : NH_max
    use Geometry, only : NH_min
    use Geometry, only : calculate_wavenumbers
    use External_Field, only : derive_transform_tensor
    use Analysis, only : spectrum_project_horizontal
    implicit none
    character(LEN=*), intent(in) :: variable_name

    double precision :: A(1:3, 1:3)
    double precision :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: data_tmp(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    double precision :: data_(0:NH_max, 0:NH_min*2)

    double precision :: data_sum
    integer :: i

    integer :: record_length
    character(LEN=LEN_FILE_PATH) :: file_path

    record_length = (NH_max + 1) * (NH_min*2 + 1)

    if (variable_name == 'KE') then
      call derive_transform_tensor(current_time, A)
      call calculate_wavenumbers(A, K)
        data_tmp(:,:,:) = 0.d0
      do i = 1, 3
        data_tmp(:,:,:) = data_tmp(:,:,:) + abs(state%velocity(i,:,:,:))**2
      enddo
      data_tmp(:,:,:) = data_tmp(:,:,:) * 0.5d0
      call spectrum_project_horizontal(data_tmp, K, data_)
    endif

    if (variable_name == 'PE') then
      call derive_transform_tensor(current_time, A)
      call calculate_wavenumbers(A, K)
      data_tmp(:,:,:) = abs(state%temperature(:,:,:))**2 * 0.5d0
      call spectrum_project_horizontal(data_tmp, K, data_)
    endif

    if (variable_name == 'NKE')  &
      &  data_(:,:) = spectrum_H%kinetic_energy_transfer(:,:)
    if (variable_name == 'NPE')  &
      &  data_(:,:) = spectrum_H%potential_energy_transfer(:,:)
    if (variable_name == 'PKE')  &
      &  data_(:,:) = spectrum_H%kinetic_energy_production(:,:)
    if (variable_name == 'PPE')  &
      &  data_(:,:) = spectrum_H%potential_energy_production(:,:)
    if (variable_name == 'CKP')  &
      &  data_(:,:) = spectrum_H%energy_conversion(:,:)
    if (variable_name == 'DKE')  &
      &  data_(:,:) = spectrum_H%kinetic_energy_dissipation(:,:)
    if (variable_name == 'DPE')  &
      &  data_(:,:) = spectrum_H%potential_energy_dissipation(:,:)

    if (my_rank == 0) then

      write(file_path, '(A, "/Spec_H2_", A, ".out")')  &
        &  trim(data_directory), trim(variable_name)

      open(unit=UNIT_FILE_DEFAULT, file=file_path, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * record_length * FACTOR_SINGLE_PRECISION)

      write(UNIT_FILE_DEFAULT, rec=file_number+1) real(data_)

      close(UNIT_FILE_DEFAULT)

      data_sum = sum(data_(0, :)) + sum(data_(1:, :))*2

    endif

    return
  end subroutine output_spectrum_horizontal

! *****************************************************************************

  subroutine output_Richardson_distribution()
    use Analysis, only : analysis_Richardson_distribution
    use Geometry, only : N_Ri
    implicit none
    double precision :: Richardson_distribution(0:N_Ri-1)
    character(LEN=LEN_FILE_PATH) :: file_path

    call analysis_Richardson_distribution(state%velocity(:,:,:,:),  &
      &  state%temperature(:,:,:), current_time, Richardson_distribution(:))

    if (my_rank == 0) then

      write(file_path, '(A, "/Richardson_distribution.out")')  &
        &  trim(data_directory)

      open(unit=UNIT_FILE_DEFAULT, file=file_path, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * N_Ri * FACTOR_SINGLE_PRECISION)

      write(UNIT_FILE_DEFAULT, rec=file_number+1) real(Richardson_distribution)

      close(UNIT_FILE_DEFAULT)

    endif

    return
  end subroutine output_Richardson_distribution

! *****************************************************************************

  subroutine output_Froude_spectrum()
    use mpi, only: MPI_DOUBLE_PRECISION
    use Parameters, only : N => buoyancy_frequency
    use Geometry, only : calculate_wavenumbers
    use Geometry, only : n3_truncated
    use External_Field, only : derive_transform_tensor
    implicit none
    double precision :: Froude_spectrum_local(0:n3_truncated)
    double precision :: Froude_spectrum(0:n3_truncated)

    double precision :: A(1:3, 1:3)
    double precision :: K(1:3, 0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
    character(LEN=LEN_FILE_PATH) :: file_path

    integer :: i
    integer :: ierr

    call derive_transform_tensor(current_time, A)
    call calculate_wavenumbers(A, K)

    Froude_spectrum_local(0) = 0.d0
    do i = 1, n3_truncated
      Froude_spectrum_local(i) = sum(K(3,:,:,i)**2  &
        &  * (abs(state%velocity(1,:,:,i))**2  &
        &  + abs(state%velocity(2,:,:,i))**2 &
        &  + abs(state%velocity(1,:,:,n3-i))**2  &
        &  + abs(state%velocity(2,:,:,n3-i))**2)) / N**2
    enddo

    call MPI_Reduce(Froude_spectrum_local(0), Froude_spectrum(0),  &
      &  n3_truncated+1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,  &
      &  ierr)

    if (my_rank == 0) then

      write(file_path, '(A, "/Froude_spectrum.out")')  &
        &  trim(data_directory)

      open(unit=UNIT_FILE_DEFAULT, file=file_path, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * (n3_truncated+1) * FACTOR_SINGLE_PRECISION)

      write(UNIT_FILE_DEFAULT, rec=file_number+1) real(Froude_spectrum)

      close(UNIT_FILE_DEFAULT)

    endif

    return
  end subroutine output_Froude_spectrum

! *****************************************************************************

  subroutine output_Thorpe_distribution()
    use Geometry, only : N_Thorpe
    use Analysis, only : Thorpe_displacements
    implicit none
    double precision :: PDF(0:N_Thorpe-1)
    character(LEN=LEN_FILE_PATH) :: file_path

    call Thorpe_displacements(state%temperature, PDF)

    if (my_rank == 0) then

      write(file_path, '(A, "/Thorpe_displacement.out")')  &
        &  trim(data_directory)

      open(unit=UNIT_FILE_DEFAULT, file=file_path, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * N_Thorpe * FACTOR_SINGLE_PRECISION)

      write(UNIT_FILE_DEFAULT, rec=file_number+1) real(PDF)

      close(UNIT_FILE_DEFAULT)

    endif

    return
  end subroutine output_Thorpe_distribution

end module Input_Output

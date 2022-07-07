program main
  use Floquet_Control, only : N_time
  use Floquet_Control, only : N => buoyancy_frequency
  use Floquet_Control, only : Omega => angular_velocity
  use Floquet_Governing_Equation, only : equation
  use Floquet_Control, only : read_control_file
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
  implicit none

  double precision, parameter :: PI = acos(-1.d0)
  double precision, parameter :: K0 = 1.d0
  integer, parameter :: NZ = 4
  integer, parameter :: UNIT_DATA = 10

  double precision, allocatable :: input_data(:,:,:)

  double precision :: zeta(1:NZ)
  double precision :: zeta_tmp(1:NZ) = 0.d0
  double precision :: term(1:NZ, 4) = 0.d0

  double precision :: zeta_start(1:NZ)
  double precision :: zeta_end(1:NZ)

  double precision :: log_M_min
  double precision :: log_M_max
  double precision :: log_M

  double precision :: Ro
  double precision :: e
  double precision :: alpha_plus
  double precision :: alpha_minus
  double precision :: ax_ratio

  double precision :: current_time = 0.d0
  double precision :: wavenumber_frequency
  double precision :: end_time
  double precision :: tmp_time = 0.d0
  double precision :: delta_time

  double precision :: D(1:3, 1:3)
  double precision :: M(1:3)
  double precision :: K(1:3)

  integer :: in
  integer :: ir
  integer :: ie
  integer :: im
  integer :: it

  ! integer :: im_array(1:3) = (/ 43937, 36154, 31680 /) ! Unstable modes
  integer :: im_array(1:3) = (/ 29189, 29944, 29819 /) ! Neutral modes
  ! integer :: ir_array(1:3) = (/ 15, 15, 15 /) ! Unstable modes
  integer :: ir_array(1:3) = (/ 15, 17, 18 /) ! Neutral modes
  integer :: j

  logical :: choose_unstable = .False.

  character(200) :: input_file
  character(200) :: output_file

  call read_control_file()
  allocate(input_data(1:NZ, 1:NE, 1:NM), source=0.d0)

  in = 9 + 1

  ! Unstable modes
  ! ie = 12 + 1

  ! Neutral modes
  ie = 16 + 1

  do j = 1, 3
  im = im_array(j) + 1
  ir = ir_array(j) + 1

  N = 10.d0

  log_M_min = log(M_min)
  log_M_max = log(M_max)
  log_M = log_M_min + (log_M_max - log_M_min) / (NM - 1) * (im - 1)
  K(3) = exp(log_M)

  Ro = Ro_min + (Ro_max - Ro_min) / (NR - 1) * (ir - 1)
  e = e_min + (e_max - e_min) / (NE - 1) * (ie - 1)
  alpha_minus = Ro * (1 - e)**2 / (1 + (1 - e)**2)
  alpha_plus = Ro / (1 + (1 - e)**2)
  ax_ratio = sign(sqrt(alpha_plus / alpha_minus), alpha_plus)

  D(1:3, 1:3) = 0.d0
  D(2,1) = alpha_plus
  D(1,2) = - alpha_minus
  M(1:3) = 0.d0

  wavenumber_frequency = sqrt(alpha_minus * alpha_plus)
  end_time = 2 * PI / wavenumber_frequency
  delta_time = end_time / N_time
  current_time = 0.d0

  write(*,*)
  write(*,*) in, ir, ie, im, choose_unstable

  if (choose_unstable) then

    write(input_file, '(A, i6.6, A)')  &
      &  '../data/floquet_sparse_2/unstable_R', in-1, '.out'
    open(UNIT_DATA, file=input_file, access='direct',  &
      &  status='old', form='unformatted',  &
      &  recl=8*NZ*NM*NE)
    read(UNIT_DATA, REC=ir) input_data
    close(UNIT_DATA)
    zeta(1:NZ) = input_data(1:NZ, ie, im)

    write(output_file,  &
      &  '(A, i6.6, "_", i6.6, "_", i6.6, "_", i6.6, A)')  &
      &  '../data/floquet_sparse_2/unstable_series_',  &
      &  (in-1), (ir-1), (ie-1), (im-1), '.out'
    open(UNIT_DATA, file=output_file, access='direct',  &
      &  status='unknown', form='unformatted', recl=8*NZ)

  else

    write(input_file, '(A, i6.6, A)')  &
      &  '../data/floquet_sparse_2/neutral_R', in-1, '.out'
    open(UNIT_DATA, file=input_file, access='direct',  &
      &  status='old', form='unformatted',  &
      &  recl=8*NZ*NM*NE)
    read(UNIT_DATA, REC=ir) input_data
    close(UNIT_DATA)
    zeta(1:NZ) = input_data(1:NZ, ie, im)

    write(output_file,  &
      &  '(A, i6.6, "_", i6.6, "_", i6.6, "_", i6.6, A)')  &
      &  '../data/floquet_sparse_2/neutral_series_',  &
      &  (in-1), (ir-1), (ie-1), (im-1), '.out'
    open(UNIT_DATA, file=output_file, access='direct',  &
      &  status='unknown', form='unformatted', recl=8*NZ)

  endif

  write(*,*) Omega, N, e, Ro, end_time
  zeta_start(:) = zeta(:)

  do it = 0, N_time - 1
    write(UNIT_DATA, REC=it+1) zeta

    zeta_tmp(:) = zeta(:)
    tmp_time = current_time
    K(1) = K0 * cos(wavenumber_frequency * tmp_time)
    K(2) = - ax_ratio * K0 * sin(wavenumber_frequency * tmp_time)
    call equation(N, Omega, D, M, K, zeta_tmp(1:3), zeta_tmp(4), term(:,1))

    zeta_tmp(:) = zeta(:) + term(:,1) * delta_time * 0.5
    tmp_time = current_time + delta_time * 0.5
    K(1) = K0 * cos(wavenumber_frequency * tmp_time)
    K(2) = - ax_ratio * K0 * sin(wavenumber_frequency * tmp_time)
    call equation(N, Omega, D, M, K, zeta_tmp(1:3), zeta_tmp(4), term(:,2))

    zeta_tmp(:) = zeta(:) + term(:,2) * delta_time * 0.5
    tmp_time = current_time + delta_time * 0.5
    K(1) = K0 * cos(wavenumber_frequency * tmp_time)
    K(2) = - ax_ratio * K0 * sin(wavenumber_frequency * tmp_time)
    call equation(N, Omega, D, M, K, zeta_tmp(1:3), zeta_tmp(4), term(:,3))

    zeta_tmp(:) = zeta(:) + term(:,3) * delta_time
    tmp_time = current_time + delta_time
    K(1) = K0 * cos(wavenumber_frequency * tmp_time)
    K(2) = - ax_ratio * K0 * sin(wavenumber_frequency * tmp_time)
    call equation(N, Omega, D, M, K, zeta_tmp(1:3), zeta_tmp(4), term(:,4))

    zeta(:) = zeta(:)  &
      &  + (term(:,1) + 2 * term(:,2) + 2 * term(:,3) + term(:,4))  &
      &  * delta_time / 6

    current_time = current_time + delta_time

  end do
  write(UNIT_DATA, REC=N_time+1) zeta

  close(UNIT_DATA)
  zeta_end(:) = zeta(:)

  write(*,*) end_time
  write(*,*) zeta_end(:) / zeta_start(:)

  enddo

  stop
end program main

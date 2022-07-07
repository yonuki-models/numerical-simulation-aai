module FFTW_Interface
  use mpi, only : MPI_COMM_WORLD
  use mpi, only : MPI_DOUBLE_PRECISION
  use omp_lib
  use Parallel, only : total_process
  use Parallel, only : comm_1
  use Parallel, only : comm_2
  use Parallel, only : rank_1
  use Parallel, only : rank_2
  use Parallel, only : size_1
  use Parallel, only : size_2
  use misk, only : assert
  implicit none
  private

  integer :: kx1
  integer :: kx2
  integer :: kx2p
  integer :: kx3p
  integer :: kz1
  integer :: kz2
  integer :: kz3
  integer :: nz1b
  integer :: nz1e
  integer :: nz2b
  integer :: nz2e
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: isn
  integer :: idir
  integer :: icon
  integer(8) :: nw

  integer :: count = 0

  integer, public :: n1C_local
  integer, public :: n2C_local
  integer, public :: n2R_local
  integer, public :: n3R_local

  complex(kind(0d0)), parameter :: UI = (0.d0, 1.d0)

  double precision, allocatable :: x(:,:,:)
  double precision, allocatable :: w(:)
  complex(kind(0d0)), allocatable :: z(:,:,:)

  public fftw_init_parallel
  public fftw_backward_parallel
  public fftw_forward_parallel
  public fftw_finalize_parallel

contains

! *****************************************************************************

subroutine fftw_init_parallel(number_of_grid, division)
  implicit none
  integer, intent(in) :: number_of_grid(3)
  integer, intent(in) :: division(2)

  n1 = number_of_grid(1)
  n2 = number_of_grid(2)
  n3 = number_of_grid(3)
  n1C_local = (n1 / 2) / division(1) + 1
  n2C_local = (n2 - 1) / division(2) + 1
  n2R_local = (n2 - 1) / division(1) + 1
  n3R_local = (n3 - 1) / division(2) + 1

  kx1 = (n1 / 2 + 1) * 2
  kx2p = n2R_local
  kx2 = n2R_local
  kx3p = n3R_local

  allocate(x(1:kx1, 1:kx2, 1:kx3p))
  nw = 0
  call ds_v3drcf2x(x, kx1, kx2, kx2p, kx3p, z, kz1, kz2, kz3,  &
    &  nz1b, nz1e, nz2b, nz2e, n1, n2, n3, w, nw, isn, idir,  &
    &  comm_1, comm_2, icon)

  allocate(z(1:kz1, 1:kz2, 1:kz3))
  allocate(w(1:nw))

  return
end subroutine fftw_init_parallel

! *****************************************************************************

subroutine fftw_forward_parallel(Q1, Q2)
  implicit none
  double precision, intent(in) :: Q1(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)
  complex(kind(0d0)), intent(out) :: Q2(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

  x(:,:,:) = 0.d0
  z(:,:,:) = (0.d0, 0.d0)

  x(1:n1, 1:n2R_local, 1:n3R_local)  &
    &  = Q1(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)

  idir = 1
  isn = 1
  call ds_v3drcf2x(x, kx1, kx2, kx2p, kx3p, z, kz1, kz2, kz3,  &
    &  nz1b, nz1e, nz2b, nz2e, n1, n2, n3, w, nw, isn, idir,  &
    &  comm_1, comm_2, icon)

  count = count + 1

  Q2(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)  &
    &  = z(1:n1C_local, 1:n2C_local, 1:n3) / (dble(n1) * n2 * n3)
    ! Avoid integer overflow

  return
end subroutine fftw_forward_parallel

! *****************************************************************************

subroutine fftw_backward_parallel(Q1, Q2)
  implicit none
  complex(kind(0d0)), intent(in) :: Q1(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)
  double precision, intent(out) :: Q2(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)

  x(:,:,:) = 0.d0
  z(:,:,:) = (0.d0, 0.d0)

  z(1:n1C_local, 1:n2C_local, 1:n3)  &
    &  = Q1(0:n1C_local-1, 0:n2C_local-1, 0:n3-1)

  idir = -1
  isn = -1
  call ds_v3drcf2x(x, kx1, kx2, kx2p, kx3p, z, kz1, kz2, kz3,  &
    &  nz1b, nz1e, nz2b, nz2e, n1, n2, n3, w, nw, isn, idir,  &
    &  comm_1, comm_2, icon)

  count = count + 1

  Q2(0:n1-1, 0:n2R_local-1, 0:n3R_local-1)  &
    &  = x(1:n1, 1:n2R_local, 1:n3R_local)

  return
end subroutine fftw_backward_parallel

! *****************************************************************************

subroutine fftw_finalize_parallel()
  implicit none

  deallocate(x, z, w)

  return
end subroutine fftw_finalize_parallel

end module FFTW_Interface

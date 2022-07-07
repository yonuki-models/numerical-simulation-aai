module Parallel
  use mpi, only : MPI_THREAD_FUNNELED
  use mpi, only : MPI_COMM_WORLD
  use Misk, only : assert
  implicit none
  private
  integer, public :: total_process
  integer, public :: my_rank
  integer, public :: comm_1
  integer, public :: comm_2
  integer, public :: rank_1
  integer, public :: rank_2
  integer, public :: size_1
  integer, public :: size_2

  public :: start_mpi
  public :: end_mpi

  contains

! *****************************************************************************

  subroutine start_mpi(division)
    implicit none
    integer, intent(in) :: division(1:2)
    integer :: provided_mpi_thread
    integer :: ierr

    call MPI_Init_thread(MPI_THREAD_FUNNELED, provided_mpi_thread, ierr)
    call assert(ierr==0, "MPI_Init_thread falied.")

    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call assert(ierr==0, "MPI initialization falied.")
    call MPI_Comm_size(MPI_COMM_WORLD, total_process, ierr)
    call assert(ierr==0, "MPI initialization falied.")

    rank_2 = my_rank / division(1)
    rank_1 = mod(my_rank, division(1))

    call MPI_Comm_split(MPI_COMM_WORLD, rank_2, my_rank, comm_1, ierr)
    call assert(ierr==0, "MPI initialization falied.")
    call MPI_Comm_size(comm_1, size_1, ierr)
    call assert(ierr==0, "MPI initialization falied.")

    call MPI_Comm_split(MPI_COMM_WORLD, rank_1, my_rank, comm_2, ierr)
    call assert(ierr==0, "MPI initialization falied.")
    call MPI_Comm_size(comm_2, size_2, ierr)
    call assert(ierr==0, "MPI initialization falied.")

    return
  end subroutine start_mpi

! *****************************************************************************

  subroutine end_mpi()
    implicit none
    integer :: ierr

    call MPI_FINALIZE(ierr)
    call assert(ierr==0, "MPI finalization falied.")

    return
  end subroutine end_mpi

end module Parallel

module Parallel
  use mpi, only : MPI_THREAD_FUNNELED
  use mpi, only : MPI_COMM_WORLD
  use Misk, only : assert
  implicit none
  private
  integer, public :: total_process
  integer, public :: my_rank

  public :: start_mpi
  public :: end_mpi

  contains

! *****************************************************************************

  subroutine start_mpi()
    implicit none
    integer :: provided_mpi_thread
    integer :: ierr

    call MPI_Init_thread(MPI_THREAD_FUNNELED, provided_mpi_thread, ierr)
    call assert(ierr==0, "MPI_Init_thread falied.")

    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call assert(ierr==0, "MPI initialization falied.")
    call MPI_Comm_size(MPI_COMM_WORLD, total_process, ierr)
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

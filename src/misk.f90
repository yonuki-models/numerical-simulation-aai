module Misk
  implicit none
  private
  integer, parameter :: UNIT_STDERR = 0
  integer, parameter :: UNIT_STDOUT = 6

  public :: assert

  contains

! *****************************************************************************

  subroutine assert(should_be_true, message, continuation)
    implicit none
    logical, intent(in) :: should_be_true
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: continuation

    if (.not. should_be_true) then
      if (present(continuation)) then
        if (continuation) then
          write(UNIT_STDERR, *) 'Warning: ', message
        endif
      else
        write(UNIT_STDERR, *) 'Fatal error: ', message
        stop
      endif
    endif

    return
  end subroutine assert

end module Misk

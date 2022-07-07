program main
  use Time_Module, only : continuation
  use Simulation, only : initial_process
  use Simulation, only : step_forward
  use Simulation, only : final_process
  implicit none

  call initial_process()

  do while(continuation())
    call step_forward()
  enddo

  call final_process()

  stop
end program main

module Floquet_Eigen_Problem
  implicit none

  integer :: NM = 2

contains

  subroutine eigen_solver(NZ, zeta, lambda, eigen_zeta)
    implicit none
    integer, intent(in) :: NZ
    double precision, intent(in) :: zeta(1:NZ, 1:NZ)
    complex(kind(0d0)), intent(out) :: lambda(1:NZ)
    complex(kind(0d0)), intent(out) :: eigen_zeta(1:NZ, 1:NZ)

    complex(kind(0d0)) :: zeta_complex(1:NZ, 1:NZ)
    complex(kind(0d0)) :: VL(1:NZ,1:NZ)
    complex(kind(0d0)) :: WORK(1:2*NZ)
    double precision :: RWORK(1:2*NZ)

    integer :: info

    zeta_complex(:,:) = zeta(:,:)
    call zgeev('N', 'V', NZ, zeta_complex, NZ, lambda, VL, NZ, eigen_zeta,  &
      &  NZ, WORK, 2*NZ, RWORK, info)

    return
  end subroutine eigen_solver

end module Floquet_Eigen_Problem

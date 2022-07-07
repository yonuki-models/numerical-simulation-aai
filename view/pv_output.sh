#!/bin/sh

#PJM -L rscgrp=regular-o
#PJM -L node=16
#PJM --mpi proc=768
#PJM --omp thread=1
#PJM -L elapse=6:00:00
#PJM -g gs53
#PJM -j
#PJM -m b
#PJM -m e
#PJM -m r
#PJM -s

module load fj/1.2.34
module load fjmpi/1.2.34
module load paraview/5.9.1
module load llvm
module load mesa

# export PV_USE_TRANSMIT=1
mpiexec -np ${PJM_MPI_PROC} pvbatch draw_Rho_case2.py
# pvbatch draw_Q.py

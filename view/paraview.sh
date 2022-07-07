#!/bin/sh

#PJM -L rscgrp=regular-o
#PJM -L node=16
#PJM --mpi proc=768
#PJM --omp thread=1
#PJM -L elapse=48:00:00
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

# mpiexec -np ${PJM_MPI_PROC} -of log.pvserver pvserver
mpiexec -np ${PJM_MPI_PROC} pvserver

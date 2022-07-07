#!/bin/sh

#PJM -L rscgrp=regular-a
#PJM -L node=3
#PJM --mpi proc=96
#PJM --omp thread=1
#PJM -L elapse=6:00:00
#PJM -g gs53
#PJM -j
#PJM -m b
#PJM -m e
#PJM -m r
#PJM -s

module load gcc
module load ompi

mpirun -machinefile $PJM_O_NODEINF -n ${PJM_MPI_PROC} -npernode 32 python frequency_analysis.py

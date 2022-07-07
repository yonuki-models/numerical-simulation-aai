#!/bin/sh

#PJM -L rscgrp=regular-a
#PJM -L node=2
#PJM --mpi proc=48
#PJM --omp thread=1
#PJM -L elapse=24:00:00
#PJM -g gs53
#PJM -j
#PJM -m b
#PJM -m e
#PJM -m r
#PJM -s

module load gcc
module load ompi

mpirun -machinefile $PJM_O_NODEINF -n ${PJM_MPI_PROC} -npernode 32 python export_vtk_case1.py

#!/bin/sh

#PJM -L rscgrp=regular-o
#PJM -L node=48:torus
#PJM --mpi proc=2304
#PJM --omp thread=1
#PJM -L elapse=24:00:00
#PJM -g gs53
#PJM -j
#PJM -m b
#PJM -m e
#PJM -m r
#PJM -s

module load mpi-fftw

data_dir='../data/'
export experiment="exRo10N3e06_LES"
export output_dir=${data_dir}${experiment}
output_dir=${data_dir}${experiment}
mkdir -p ${output_dir}
mkdir -p ${output_dir}/state
mkdir -p ${output_dir}/energy
mkdir -p ${output_dir}/snapshot
echo Output directory: ${output_dir}
export dir_name_result=${output_dir}
mpiexec ./main.exe

#!/bin/bash

#PJM -L rscgrp=regular-o
#PJM -L node=42
#PJM --mpi proc=2000
#PJM --omp thread=1
#PJM -L elapse=24:00:00
#PJM -g gs53
#PJM -j
#PJM -m b
#PJM -m e
#PJM -m r
#PJM -s

export experiment="floquet"
output_dir="../data/"${experiment}
mkdir ${output_dir}
echo Output directory: ${output_dir}
mpiexec ./floquet_sweep.exe

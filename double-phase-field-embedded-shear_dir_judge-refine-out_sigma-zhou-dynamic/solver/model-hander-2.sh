#!/bin/bash
set +e
echo '======================================================================='
echo 'd-modle-1 start'
mpirun -n 20 ./solver 3 parameters/embedded-3.prm
mpirun -n 20 ./solver 3 parameters/embedded-4.prm
mpirun -n 20 ./solver 3 parameters/embedded-5.prm
echo '======================================================================='

#!/bin/bash
set +e
echo '======================================================================='
echo 'd-modle-1 start'
mpirun -n 10 ./solver 3 parameters/embedded-5.prm
mpirun -n 10 ./solver 3 parameters/embedded-6.prm
echo '======================================================================='

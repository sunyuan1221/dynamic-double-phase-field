#!/bin/bash
set +e
echo '======================================================================='
echo 'd-modle-1 start'
mpirun -n 20 ./solver 3 parameters/embedded-1.prm
mpirun -n 20 ./solver 3 parameters/embedded-2.prm
echo '======================================================================='

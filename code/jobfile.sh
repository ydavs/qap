#!/bin/bash

#Job Name 
#PBS -N qap_tai30a

#Walltime 
#PBS -l walltime=24:00:00
#Queue 
#PBS -q p-queue

### Merging output and error files 
#PBS -j oe
#PBS -o output.log

##Selecting nodes and cpu's and processors 
#PBS -l select=1:ncpus=32

#Sending email 
#PBS -m abe
#PBS -M sourav17@iiserb.ac.in

#Copying all the required stuff to the scratch 
SCRATCH_DIR=/home2/paurora/scratch 
cp /home2/paurora/qap/* ${SCRATCH_DIR}
cd ${SCRATCH_DIR}

make
./main model.mps ins.dat > tmp

cp ${SCRATCH_DIR}/tmp /home2/paurora/qap/
cp ${SCRATCH_DIR}/output.log /home2/paurora/qap/

rm ${SCRATCH_DIR}/*

exit 0 

#!/bin/bash

WS=/tmp/workspace/


# source amber.sh to be able to find some program from the current running shell 
source /tmp/amber20/amber.sh 

# add gromacs path to be able to find gromacs 
export GMX=/usr/local/gromacs
export PATH="$GMX/bin:$PATH"

cd $WS/scripts/pythonScript

python3 MDBatch.py -f $WS/MD/ALDOSE -g gmx_mpi -a /tmp/amber20/miniconda/bin/gmx_MMPBSA -m 2 -d


exit;

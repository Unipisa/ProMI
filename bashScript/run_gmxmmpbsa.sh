#!/bin/bash

# Set some environment variables 
WS=/tmp/workspace/

cd $WS/MD/1hknew

echo "The working directory is:"
pwd

# source amber.sh to be able to find some program from the current running shell 
source /tmp/amber20/amber.sh 

# add gromacs path to be able to find gromacs 
export GMX=/usr/local/gromacs
export PATH="$GMX/bin:$PATH"

echo $PATH

# rm -rf 1HK1
# rm -rf 1HK2 
# rm -rf 1HK3

# /tmp/amber20/miniconda/bin/gmx_MMPBSA MPI -O -i mmpbsa.in -cs md.tpr -ci index.ndx -cg 1 13 -ct md_fit.xtc  -cp topol.top 

# /tmp/amber20/miniconda/bin/gmx_MMPBSA -h
# gmx_mpi -h 

exit;

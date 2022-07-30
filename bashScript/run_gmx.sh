#!/bin/bash

# Set some environment variables 
WS=/tmp/workspace

cd $WS/HIV_PROTEASE_1/MD/V82N/

echo "Current directory:"
pwd

export GMX=/usr/local/gromacs
export PATH="$GMX/bin:$PATH"

# find / -iname *gmx*
# gmx_mpi grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
# gmx_mpi mdrun -v -deffnm em

gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt -nb gpu 

gmx_mpi grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 2
gmx_mpi mdrun -deffnm npt -nb gpu

gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr -maxwarn 2
gmx_mpi mdrun -deffnm md -nb gpu

exit;

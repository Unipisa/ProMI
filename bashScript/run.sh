#!/bin/bash

# Set some environment variables 
WS=/tmp/workspace/

cd $WS/PredictedStructures/DIHYDROFOLATE/WT

# find / -iname chimerax

/usr/bin/chimerax -h

/usr/bin/chimerax  --nogui ranked_0.pdb

# cd $WS/ligands

# rm -rf NDP


exit;

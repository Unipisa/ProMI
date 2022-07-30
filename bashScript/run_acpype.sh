#!/bin/bash

# Set some environment variables 
WS=/tmp/workspace/

cd $WS/scripts/pythonScript/test/4enz/MD/Share
# acpype -f -i CA1.pdb -n 0 -o gmx -d

acpype -f -i OXY.pdb -o gmx -d 


# cd $WS/ligands/NAP/
# echo "Current directory:"
# pwd
# acpype -i NAP.pdb -n 0 -o gmx -d

# echo "============================="

# cd $WS/ligands/47D/
# echo "Current directory:"
# pwd
# acpype -i 47D.pdb -n 0 -o gmx -d
# echo "============================="

exit;

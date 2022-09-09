#!/bin/sh
WS=/tmp/workspace/
cd $WS
acpype -f -i NDP.pdb -o gmx -d
acpype -f -i XCF.pdb -o gmx -d

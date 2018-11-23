#!/bin/bash
CURRENT=`pwd`
#MPMFILE="$CURRENT/dem2mpm.dat"
MPMFILE="$CURRENT/$1"
cd /home/fabricio/eclipse-workspace-oxigen-cdt/MPM-Generator/Release
time ./MPM-Generator $MPMFILE 
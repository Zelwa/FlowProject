#!bin/bash

# MODIFIER LE CHEMIN VERS ORTOOL

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/blop/Documents/ORtoolsDir/lib"
export LIBRARY_PATH="$LIBRARY_PATH:/home/blop/Documents/ORtoolsDir/lib"
#make build ORTOOLS=/home/blop/Documents/ORtoolsDir SOURCE=src/main.cc
make run ORTOOLS=/home/blop/Documents/ORtoolsDir SOURCE=src/main.cc ARGS=$1

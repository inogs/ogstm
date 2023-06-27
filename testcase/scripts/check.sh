#! /usr/bin/env bash

ROOT="$CINECA_SCRATCH/BENCHMARKS"
BUILD=PELCHEM_BUILD

cd $ROOT || exit

./build.sh --fast --verbose --build-path="$BUILD" || exit

cd "$ROOT/ogstm/testcase/$1" || exit

source "$ROOT/ogstm/compilers/machine_modules/m100.hpc-sdk"
module unload numpy python
export RANKS_PER_NODE=1

mpirun -gpu -np 1 ./ogstm.xx || exit

cd "$ROOT/ogstm/testcase" || exit

conda run -n esiwace python scripts/comparedatasets.py --rtol=0.0 $1 REFERENCE AVE_FREQ_2



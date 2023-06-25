#! /usr/bin/env bash

ROOT="$CINECA_SCRATCH/BENCHMARKS"
BUILD=PELCHEM_BUILD
TEST=TEST

cd $ROOT || exit

./build.sh --fast --verbose --build-path="$BUILD"

cd "$ROOT/ogstm/testcase/$TEST" || exit

source "$ROOT/ogstm/compilers/machine_modules/m100.hpc-sdk"
module unload numpy python
export RANKS_PER_NODE=1

mpirun -gpu -np 1 ./ogstm.xx

cd "$ROOT/ogstm/testcase" || exit

conda run -n esiwace python comparedatasets.py --rtol=0.0 TEST REFERENCE AVE_FREQ_2



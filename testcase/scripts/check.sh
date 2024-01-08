#! /usr/bin/env bash

ROOT="$CINECA_SCRATCH/ModelBuild"
BUILD=GPU_BUILD

cd $ROOT || exit

./build.sh --fast --verbose --build-path="$BUILD" || exit

cd "$ROOT/ogstm/testcase/$1" || exit

source "$ROOT/ogstm/compilers/machine_modules/leonardo.nvhpc"
export RANKS_PER_NODE=1

ln -sf "$ROOT/$BUILD/ogstm.xx" || exit

mpirun -np 1 ./ogstm.xx || exit

cd "$ROOT/ogstm/testcase" || exit

conda run -n ogstm-bfm python scripts/comparedatasets.py --atol=0.0 --rtol=0.0 $1 REFERENCE AVE_FREQ_2



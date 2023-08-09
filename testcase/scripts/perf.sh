#! /usr/bin/env bash

ROOT="$CINECA_SCRATCH/BENCHMARKS"
BUILD=GPU_BUILD

cd "$ROOT/ogstm/testcase/$1" || exit

source "$ROOT/ogstm/compilers/machine_modules/leonardo.nvhpc"
export RANKS_PER_NODE=1

ln -sf "$ROOT/$BUILD/ogstm.xx" || exit

mpirun -gpu -np 1 perf record -g ./ogstm.xx || exit

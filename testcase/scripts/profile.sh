ROOT=$CINECA_SCRATCH/BENCHMARKS

cd $ROOT/ogstm/testcase/$1 || exit
source $ROOT/ogstm/compilers/machine_modules/m100.hpc-sdk

module unload numpy python
module load cuda/11.4

export RANKS_PER_NODE=1

rm -rf /tmp/nvidia
ln -s $TMPDIR /tmp/nvidia
mpirun -gpu -np 1 nsys profile -o $PWD/$(date -Iminute) -f true --stats=true --trace=openacc ./ogstm.xx
rm -rf /tmp/nvidia

#! /usr/bin/env bash

cd $CINECA_SCRATCH/ModelBuild/ogstm/testcase/GPU_BENCHMARK || exit

#find . -name "namelist.init" -exec sed -i "s/nwritetrc = 10000/nwritetrc = 96/g" {} +

for DIR in $(find . -mindepth 1 -maxdepth 1 -type d); do
    echo "Processing $DIR"
    if [[ -d $DIR/logs ]]; then
        rm -r $DIR/logs || exit
    fi
    mkdir $DIR/logs && cd $DIR/logs || exit
    TIME=$(( $(basename $DIR) ** 2 / 3000 ))
    TIME=$(($TIME < 1 ? 1 : $TIME))
    TIME=$(($TIME < 1440 ? $TIME : 1440))
    sbatch --array=1-10 --time $TIME --output %A_%a.log ../../../scripts/job.slurm ..
    cd ../..
done


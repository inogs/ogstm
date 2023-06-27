#! /usr/bin/env bash

find . -name "namelist.init" -exec sed -i "s/nwritetrc = 10000/nwritetrc = 128/g" {} +

for DIR in $(find . -mindepth 1 -maxdepth 1 -type d); do
    rm -r $DIR/logs
    mkdir $DIR/logs && cd $DIR/logs
    #cd $DIR/logs
    sbatch --array=1-10 --time $((5 + $(basename $DIR) ** 2 / 120)) --output %A_%a.log ../../../job.slurm ..
    cd ../..
done


#! /usr/bin/env bash

find . -name "namelist.init" -exec sed -i "s/nwritetrc = 10000/nwritetrc = 512/g" {} +

for DIR in $(find . -mindepth 1 -maxdepth 1 -type d); do
    #rm -r $DIR/logs
    mkdir $DIR/logs && cd $DIR/logs
    #cd $DIR/logs
    sbatch --array=1-10 --time $((2 * $(basename $DIR))) --output %A_%a.log ../../../job.slurm ..
    cd ../..
done


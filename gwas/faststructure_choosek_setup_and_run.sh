#!/usr/bin/env bash

while getopts i:f: flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
        f) CHOOSEK=${OPTARG};;
    esac
done

source /data/$USER/conda/etc/profile.d/conda.sh && conda activate base && conda activate fastStructure

python $CHOOSEK\
 --input=$INPUT\
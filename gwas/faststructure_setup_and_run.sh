#!/usr/bin/env bash

while getopts i:o:f:k: flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        f) FASTSTRUCTURE=${OPTARG};;
        k) K=${OPTARG};;
    esac
done

source /data/$USER/conda/etc/profile.d/conda.sh && conda activate base && conda activate fastStructure

python $FASTSTRUCTURE\
 -K $K\
 --input=$INPUT\
 --output=$OUTPUT
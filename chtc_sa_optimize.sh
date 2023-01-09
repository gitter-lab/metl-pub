#!/bin/bash
Cluster=$1
Process=$2
config=$3

echo "$Cluster"
echo "$Process"
echo "$config"

tar -zxf sa_optimize_files.tar.gz
mkdir output
cp "$config" output/

python run_sa.py $config

tar -zcf sa_optimize_output_"$Cluster"_"$Process".tar.gz ./output

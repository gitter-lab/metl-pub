#!/bin/bash


# recommended to export the home directory as 'cd' command no longer works
export HOME=$PWD

echo "$Cluster"
echo "$Process"
echo "$config"


echo "Launching $RUN_SCRIPT"
python3 "$RUN_SCRIPT" @"args/${PROCESS}.txt" --cluster="$CLUSTER" --process="$PROCESS"




# make sure output files
# tar -czf "$Cluster"_"$Process".tar.gz ./output



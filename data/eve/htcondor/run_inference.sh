#!/bin/bash

set -e

cp /staging/bcjohnson7/eve.tar.gz ./

ENVNAME=eve
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME

# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
echo "un taring environment file"
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate


echo " environment activated "


echo "removing environment of interest"
rm $ENVNAME.tar.gz 


tar -xzf EVE_GFP_AEQVI.tar.gz  


# change these
export uniprot='GFP_AEQVI'
export model_name_suffix='protein_gym'
export protein_index=4


echo "moving file over that is necessary"

NAME_OF_CSV=$(ls *.csv)
mv $NAME_OF_CSV ./EVE/data/mutations/${uniprot}.csv

cd EVE

export MSA_data_folder='./data/MSA'
export MSA_list='./data/mappings/example_mapping.csv'
export MSA_weights_location='./data/weights'
export VAE_checkpoint_location='./results/VAE_parameters'

export model_parameters_location='./EVE/default_model_params.json'
export training_logs_location='./logs/'
export mutations_location='./data/mutations/'
export computation_mode='mutations_location'
export output_evol_indices_location='./results/evol_indices'
export num_samples_compute_evol_indices=20000
export batch_size=2048



python compute_evol_indices.py \
    --MSA_data_folder ${MSA_data_folder} \
    --MSA_list ${MSA_list} \
    --protein_index ${protein_index} \
    --MSA_weights_location ${MSA_weights_location} \
    --VAE_checkpoint_location ${VAE_checkpoint_location} \
    --model_name_suffix ${model_name_suffix} \
    --model_parameters_location ${model_parameters_location} \
    --computation_mode ${computation_mode} \
    --mutations_location ${mutations_location}  \
    --output_evol_indices_location ${output_evol_indices_location} \
    --num_samples_compute_evol_indices ${num_samples_compute_evol_indices} \
    --batch_size ${batch_size}
# 1>dev/null

cd .. 

mv ./EVE/results/evol_indices/${uniprot}_${num_samples_compute_evol_indices}_samples.csv  results_$NAME_OF_CSV
# tar -czf EVE_${uniprot}_inference.tar.gz EVE/

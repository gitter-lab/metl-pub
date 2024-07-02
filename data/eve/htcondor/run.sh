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


tar -xzf EVE.tar.gz  


cd EVE

export MSA_data_folder='./data/MSA'
export MSA_list='./data/mappings/example_mapping.csv'
export MSA_weights_location='./data/weights'
export VAE_checkpoint_location='./results/VAE_parameters'
export model_parameters_location='./EVE/default_model_params.json'
export training_logs_location='./logs/'


#what has been changed
export protein_index=4
export model_name_suffix='protein_gym'
export uniprot='GFP_AEQVI'


python train_VAE.py \
    --MSA_data_folder ${MSA_data_folder} \
    --MSA_list ${MSA_list} \
    --protein_index ${protein_index} \
    --MSA_weights_location ${MSA_weights_location} \
    --VAE_checkpoint_location ${VAE_checkpoint_location} \
    --model_name_suffix ${model_name_suffix} \
    --model_parameters_location ${model_parameters_location} \
    --training_logs_location ${training_logs_location} 

cd ..

tar -czf EVE_${uniprot}.tar.gz EVE/
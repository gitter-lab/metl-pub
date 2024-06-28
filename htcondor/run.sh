#!/bin/bash

cp /staging/{chtc_id}/metl.tar.gz ./

ENVNAME=metl
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


tar -xzf code.tar.gz  


# process ID needs to go right their
python src/run.py @${UUID}.txt







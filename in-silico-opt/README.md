# Simulated annealing

METL implementation of simulated annealing adapted from Sarah Fahlberg.

### Algorithm Overview: 
1. Run simulated annealing 10,000 times for 10,000 steps for each design criteria (Unobserved 5, Unobserved 10, Observed 5, and Observed 10)
   (look to `src` directory to replicate and run simulated annealing)
2. For each design criteria, cluster the top performing variants from simulated annealing to reduce the chosen designs to only 5 variants. (`analysis/cluster.py`)

***Note***: If you are looking to replicate the clustering portion of the METL manuscript to select design criteria look to the `analysis` directory.   

## Installation 

Run the follow command from the terminal. 

```
conda create --name in_silico_opt python numpy pandas
conda activate in_silico_opt
conda install -c anaconda seaborn=>0.11.*
pip install tqdm
```

### Overview (Simulated Annealing)

- Step 1 : Specify a sequence to fitness function
which takes as an argument a string of the mutation,
along with wild type. The file must be named the same
as the function name (seq2fitness_example below).
```python
# seq2fitness_example.py
import numpy as np
def seq2fitness_example(mutant:str,WT): 
    '''
    my example sequence to fitness function
    :param mutant: string of mutants, "," as specified by arguments "E3K,G102S"
    :param WT: Wildtype sequence 
    :return: a scalar value of fitness
    '''
    
    # note your split character must align with that specified in args document
    return np.random.poisson(len(mutant.split(',')))
```

- Step 2: From inside the in_silico_opt directory run```src/run.py``` 
with proper command line arguments. 
The three required arguments are the wildtype sequence (```--wild_type```), along 
with the sequence to fitness filename (```--seq2fitness_filename```), and the unique
identifier for the run (```--uuid```). These are the only required inputs, 
many other input parameters can be specified depending on your specifications (i.e amino acid
constraints ```--AA_constraints``` or specified sampler from random, nested sampling and simulated annealing [default], ```--sampler```). For 
help on all the inputs run ```python src/run.py -h``` in the ```in_silico_opt```
directory, or look to ```src/run.py```. 

Input arguments can come from the command line or a document.

#### Specify args from command line:  
```shell
python src/run.py --seq2fitness_file seq2fitness_example
                  --uuid example_run 
                  --wild_type aypsvbetbUI 
```
Note in this circumstance residues UI would be held constant as they 
are capital. 

#### Specify args from document: 

```shell
python src/run.py @args/example_args.txt
```

## HTCondor

1. To generate argument files for a new HTCondor run, please go to `htcondor/submit.py` and change the input arguments in function `generate_args`. For each variable that is not defined, the default argument shown by the list `defaults` at the top of the file is used. The argument corresponds to the variable in the `inputs` list. All argument files will be saved according the `save_path`. 
2. Change output directory in `htcondor/submit.sub` to your run name and create those directories on submit node. Make an argument file which contains the unique ids for each run which is consistent with the names or the argument files. If looking to reproduce paper results this will be a text file with each line containing an integer 1 to 10,000.
3. Place the name of the argument file in `<argument.txt file here>`. 
```shell
# submit.sub
transfer_input_files = code.tar.gz, args/<your_run_name>/$(uuid).txt, run.sh
transfer_output_remaps= "$(uuid).csv = results/<your_run_name>/$(uuid).csv"

queue uuid from <argument.txt file here>
```
4. Replace the location of your conda environment. To make same conda environment as below, use `conda_pack`. (Step 3 in this guide: [https://chtc.cs.wisc.edu/uw-research-computing/conda-installation](URL)) Make sure to name your environment `metl`, or replace it on `run.sh`. If you are working on CHTC at UW-Madison, place tar'd conda environment on shared file system `staging` and replace below line with your `chtc_id`. 

```shell
# run.sh
cp /staging/{chtc_id}/metl.tar.gz ./

ENVNAME=metl #or your environment name
```
3. Transfer `htcondor/run.sh`, `htcondor/submit.sub`, and argument file to CHTC, make results directory consistent with run name. Submit job for simulated annealing using `condor_submit`. 

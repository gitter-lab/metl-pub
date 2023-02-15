# In_silico_optimization 

METL implementation of in_silico_opt heavily adapted from
Sarah and Chase intial implementation. 

## Installation 

Run the follow command from the terminal. 

```conda create --name in_silico_opt python numpy pandas```

```conda activate in_silico_opt```

```pip install tqdm```

## Overview 

- Step 1 : Specify a sequence to fitness function
which takes as an arugment a string of the mutation,
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

## CHTC
Still to do. Look to htcondor/run.sh, htcondor/submit.sub

## Ensembling
Still to do! 

 


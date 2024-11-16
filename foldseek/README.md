# FoldSeek clustering of PDB files
This directory contains results from running [FoldSeek](https://github.com/steineggerlab/foldseek/tree/9-427df8a) to cluster the [PDB files](../data/pdb_files) in this repository to check for overlap between the training data for METL-Global and the structures used when modeling the deep mutational scanning data.

FoldSeek was run in the miniforge Docker image [`condaforge/miniforge-pypy3:24.9.2-0`](https://hub.docker.com/layers/condaforge/miniforge-pypy3/24.9.2-0/images/sha256-7550be2f2b35a166be8db35f93244795ade3f0dd117aa6f90c650184b9277cb9) and installed with
```
# conda install -c conda-forge -c bioconda foldseek
```
This installed version `9.427df8a`:
```
# conda list | grep foldseek
foldseek                  9.427df8a       pl5321hb365157_1    bioconda

```

FoldSeek `easy-cluster` was run with default settings (see below) except for the coverage threshold, which was set to 0.5 and 0.9:
```
# foldseek easy-cluster data/pdb_files/ metl_0.5 tmp_0.5 -c 0.5 > foldseek_0.5.out
# foldseek easy-cluster data/pdb_files/ metl_0.9 tmp_0.9 -c 0.9 > foldseek_0.9.out
```

Coverage threshold 0.5 generated 130 clusters from the 159 input structures.
Coverage threshold 0.9 generated 151 clusters.
The FoldSeek documentation describes the [output](https://github.com/steineggerlab/foldseek/tree/9-427df8a?tab=readme-ov-file#output-cluster) file formats.

## Default settings
```
# foldseek easy-cluster --help
usage: foldseek easy-cluster <i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]> <o:clusterPrefix> <tmpDir> [options]
 By Martin Steinegger <martin.steinegger@snu.ac.kr>
options: prefilter:
 --seed-sub-mat TWIN              Substitution matrix file for k-mer generation [aa:3di.out,nucl:3di.out]
 -s FLOAT                         Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [4.000]
 -k INT                           k-mer length (0: automatically set to optimum) [0]
 --target-search-mode INT         target search mode (0: regular k-mer, 1: similar k-mer) [0]
 --k-score TWIN                   k-mer threshold for generating similar k-mer lists [seq:2147483647,prof:2147483647]
 --max-seqs INT                   Maximum results per query sequence allowed to pass the prefilter (affects sensitivity) [300]
 --split INT                      Split input into N equally distributed chunks. 0: set the best split automatically [0]
 --split-mode INT                 0: split target db; 1: split query db; 2: auto, depending on main memory [2]
 --split-memory-limit BYTE        Set max memory per split. E.g. 800B, 5K, 10M, 1G. Default (0) to all available system memory [0]
 --comp-bias-corr INT             Correct for locally biased amino acid composition (range 0-1) [1]
 --comp-bias-corr-scale FLOAT     Correct for locally biased amino acid composition (range 0-1) [1.000]
 --diag-score BOOL                Use ungapped diagonal scoring during prefilter [1]
 --exact-kmer-matching INT        Extract only exact k-mers for matching (range 0-1) [0]
 --mask INT                       Mask sequences in k-mer stage: 0: w/o low complexity masking, 1: with low complexity masking [1]
 --mask-prob FLOAT                Mask sequences is probablity is above threshold [0.900]
 --mask-lower-case INT            Lowercase letters will be excluded from k-mer search 0: include region, 1: exclude region [1]
 --min-ungapped-score INT         Accept only matches with ungapped alignment score above threshold [30]
 --spaced-kmer-mode INT           0: use consecutive positions in k-mers; 1: use spaced k-mers [1]
 --spaced-kmer-pattern STR        User-specified spaced k-mer pattern []
 --local-tmp STR                  Path where some of the temporary files will be created []
align:
 -c FLOAT                         List matches above this fraction of aligned (covered) residues (see --cov-mode) [0.000]
 --cov-mode INT                   0: coverage of query and target
                                  1: coverage of target
                                  2: coverage of query
                                  3: target seq. length has to be at least x% of query length
                                  4: query seq. length has to be at least x% of target length
                                  5: short seq. needs to be at least x% of the other seq. length [0]
 --sort-by-structure-bits INT     sort by bits*sqrt(alnlddt*alntmscore) [1]
 -a BOOL                          Add backtrace string (convert to alignments with mmseqs convertalis module) [0]
 --alignment-mode INT             How to compute the alignment:
                                  0: automatic
                                  1: only score and end_pos
                                  2: also start_pos and cov
                                  3: also seq.id [0]
 --alignment-output-mode INT      How to compute the alignment:
                                  0: automatic
                                  1: only score and end_pos
                                  2: also start_pos and cov
                                  3: also seq.id
                                  4: only ungapped alignment
                                  5: score only (output) cluster format [0]
 -e DOUBLE                        List matches below this E-value (range 0.0-inf) [1.000E+01]
 --min-seq-id FLOAT               List matches above this sequence identity (for clustering) (range 0.0-1.0) [0.000]
 --min-aln-len INT                Minimum alignment length (range 0-INT_MAX) [0]
 --seq-id-mode INT                0: alignment length 1: shorter, 2: longer sequence [0]
 --alt-ali INT                    Show up to this many alternative alignments [0]
 --max-rejected INT               Maximum rejected alignments before alignment calculation for a query is stopped [2147483647]
 --max-accept INT                 Maximum accepted alignments before alignment calculation for a query is stopped [2147483647]
 --gap-open TWIN                  Gap open cost [aa:10,nucl:10]
 --gap-extend TWIN                Gap extension cost [aa:1,nucl:1]
clust:
 --cluster-mode INT               0: Set-Cover (greedy)
                                  1: Connected component (BLASTclust)
                                  2,3: Greedy clustering by sequence length (CDHIT) [0]
 --max-iterations INT             Maximum depth of breadth first search in connected component clustering [1000]
 --similarity-type INT            Type of score used for clustering. 1: alignment score 2: sequence identity [2]
 --single-step-clustering BOOL    Switch from cascaded to simple clustering workflow [0]
 --cluster-steps INT              Cascaded clustering steps from 1 to -s [3]
 --cluster-reassign BOOL          Cascaded clustering can cluster sequence that do not fulfill the clustering criteria.
                                  Cluster reassignment corrects these errors [0]
kmermatcher:
 --weights STR                    Weights used for cluster priorization []
 --cluster-weight-threshold FLOAT Weight threshold used for cluster priorization [0.900]
 --kmer-per-seq INT               k-mers per sequence [21]
 --kmer-per-seq-scale TWIN        Scale k-mer per sequence based on sequence length as kmer-per-seq val + scale x seqlen [aa:0.000,nucl:0.200]
 --adjust-kmer-len BOOL           Adjust k-mer length based on specificity (only for nucleotides) [0]
 --hash-shift INT                 Shift k-mer hash initialization [67]
 --include-only-extendable BOOL   Include only extendable [0]
 --ignore-multi-kmer BOOL         Skip k-mers occurring multiple times (>=2) [0]
misc:
 --tmscore-threshold FLOAT        accept alignments with a tmsore > thr [0.0,1.0] [0.000]
 --lddt-threshold FLOAT           accept alignments with a lddt > thr [0.0,1.0] [0.000]
 --alignment-type INT             How to compute the alignment:
                                  0: 3di alignment
                                  1: TM alignment
                                  2: 3Di+AA [2]
 --exact-tmscore INT              turn on fast exact TMscore (slow), default is approximate [0]
 --tmalign-hit-order INT          order hits by 0: (qTM+tTM)/2, 1: qTM, 2: tTM, 3: min(qTM,tTM) 4: max(qTM,tTM) [0]
 --tmalign-fast INT               turn on fast search in TM-align [1]
 --rescore-mode INT               Rescore diagonals with:
                                  0: Hamming distance
                                  1: local alignment (score only)
                                  2: local alignment
                                  3: global alignment
                                  4: longest alignment fulfilling window quality criterion [0]
 --mask-bfactor-threshold FLOAT   mask residues for seeding if b-factor < thr [0,100] [0.000]
 --input-format INT               Format of input structures:
                                  0: Auto-detect by extension
                                  1: PDB
                                  2: mmCIF
                                  3: mmJSON
                                  4: ChemComp
                                  5: Foldcomp [0]
 --file-include STR               Include file names based on this regex [.*]
 --file-exclude STR               Exclude file names based on this regex [^$]
common:
 --sub-mat TWIN                   Substitution matrix file [aa:3di.out,nucl:3di.out]
 --max-seq-len INT                Maximum sequence length [65535]
 --db-load-mode INT               Database preload mode 0: auto, 1: fread, 2: mmap, 3: mmap+touch [0]
 --threads INT                    Number of CPU-cores used (all by default) [8]
 --compressed INT                 Write compressed output [0]
 -v INT                           Verbosity level: 0: quiet, 1: +errors, 2: +warnings, 3: +info [3]
 --remove-tmp-files BOOL          Delete temporary files [1]
 --force-reuse BOOL               Reuse tmp filse in tmp/latest folder ignoring parameters and version changes [0]
 --mpi-runner STR                 Use MPI on compute cluster with this MPI command (e.g. "mpirun -np 42") []
 --prostt5-model STR              Path to ProstT5 model []
expert:
 --taxon-list STR                 Taxonomy ID, possibly multiple values separated by ',' []
 --zdrop INT                      Maximal allowed difference between score values before alignment is truncated  (nucleotide alignment only) [40]
 --filter-hits BOOL               Filter hits by seq.id. and coverage [0]
 --sort-results INT               Sort results: 0: no sorting, 1: sort by E-value (Alignment) or seq.id. (Hamming) [0]
 --chain-name-mode INT            Add chain to name:
                                  0: auto
                                  1: always add
                                   [0]
 --write-mapping INT              write _mapping file containing mapping from internal id to taxonomic identifier [0]
 --coord-store-mode INT           Coordinate storage mode:
                                  1: C-alpha as float
                                  2: C-alpha as difference (uint16_t) [2]
 --write-lookup INT               write .lookup file containing mapping from internal id, fasta id and file number [1]
```

#! /bin/bash

#SBATCH -A MST108226
#SBATCH -J PSIBLAST         # Job name
#SBATCH -p ct56             # Partiotion name
#SBATCH -n 4               # Number of MPI tasks (i.e. processes)
#SBATCH -c 1                # Number of cores per MPI task
#SBATCH -N 1                # Maximum number of nodes to be allocated
#SBATCH -o %j.out           # Path to the standard output file
#SBATCH -e %j.err           # Path to the standard error ouput file
#SBATCH --mail-type=END
#SBATCH --mail-user=whchen0630@gmail.com


hint cnv -m ContactMatrix/GSE63525_K562_combined.hic -f juicer --refdir reference/hg19 -r 50 -g hg19 -n K562_DpnII --bicseq ../BICseq2-seg_v0.7.3 -e DpnII -p 0

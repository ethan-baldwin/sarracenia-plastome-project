#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ethan.baldwin@uga.edu
#SBATCH --ntasks=4
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/eab77806/logs/IQTREE.%j.out
#SBATCH --error=/scratch/eab77806/logs/IQTREE.%j.err

# create output directory if it does not exist, and change directory into it
OUTDIR="/scratch/eab77806/iqtree"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

ml IQ-TREE/2.2.2.6-gompi-2022a
iqtree2 -s alignment.fasta -nt AUTO -bb 1000 -m GTR+F+R4

#!/bin/bash
#SBATCH --job-name=pga
#SBATCH --partition=batch
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eab77806@uga.edu
#SBATCH --ntasks=8
#SBATCH --mem=12gb
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/eab77806/logs/pgalog.%j.out
#SBATCH --error=/scratch/eab77806/logs/pgalog.%j.err

# create output directory if it does not exist, and change directory into it
OUTDIR="/scratch/eab77806/genbank/pga"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

# load softwar; adding blastn to $PATH was necessary for PGA to run
ml BLAST+/2.10.1-gompi-2019b

# run PGA
perl /home/eab77806/Chloroplast/PGA/PGA.pl -r /home/eab77806/sarracenia-plastome-project/data/pgaRef -t /home/eab77806/cp_genomes_fa -o gb

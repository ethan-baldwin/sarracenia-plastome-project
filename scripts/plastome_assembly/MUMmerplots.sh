#!/bin/bash
#SBATCH --job-name=mummer
#SBATCH --partition=batch
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eab77806@uga.edu
#SBATCH --ntasks=6
#SBATCH --mem=24gb
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/eab77806/logs/mummerlog.%j.out
#SBATCH --error=/scratch/eab77806/logs/mummerlog.%j.err

# create output directory if it does not exist, and change directory into it
OUTDIR="/scratch/eab77806/assemblies"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

# load software
ml MUMmer/4.0.0beta2-foss-2019b

# path to reference
ref="/home/eab77806/referenceGenomes/RRCAN01.fasta"
# fasta file of the assembly you want to align
fasta="Uncircularized_assemblies_1_m003.fasta"

# generate mummer plot
nucmer $ref $fasta -p merged_out
delta-filter -1 merged_out.delta  >  merged_out.1delta
mummerplot --size large --layout --color -f --png merged_out.1delta -p merged_out

# remove intermediate files
# rm *delta
# rm *plot
# rm *filter
# rm *gp

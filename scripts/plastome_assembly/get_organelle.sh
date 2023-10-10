#!/bin/bash
#SBATCH --job-name=get_organelle
#SBATCH --partition=batch
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eab77806@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24gb
#SBATCH --time=024:00:00
#SBATCH --output=/scratch/eab77806/logs/get_organelle.%j.out
#SBATCH --error=/scratch/eab77806/logs/get_organelle.%j.err

# load software
ml GetOrganelle/1.7.5.2-foss-2020b

# create output directory if it does not exist, and change directory into it
OUTDIR="/scratch/eab77806/sarracenia_plastome"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

# input directory with reads
INDIR="/scratch/eab77806/trimmed_reads"
# name of sample (prefix of fastq files)
SAMPLE="m001"
# path to seed for get_organelle
SEED="/home/eab77806/references/NC_041129.fasta"

# run get_organelle 
get_organelle_from_reads.py -1 $INDIR/${SAMPLE}_P_R1.fastq.gz -2 $INDIR/${SAMPLE}_P_R2.fastq.gz  -o ${SAMPLE}.plastome_output -s $SEED -k 21,45,65,85,97,105,121 -t 12 -F embplant_pt

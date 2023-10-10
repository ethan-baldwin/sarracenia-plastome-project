#!/bin/bash
#SBATCH --job-name=gb2gff
#SBATCH --partition=batch
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eab77806@uga.edu
#SBATCH --ntasks=2
#SBATCH --mem=20gb
#SBATCH --time=8:00:00
#SBATCH --output=/scratch/eab77806/logs/gb2gff.%j.out
#SBATCH --error=/scratch/eab77806/logs/gb2gff.%j.err

# create output directory if it does not exist, and change directory into it
OUTDIR="/scratch/eab77806/assemblies/consensus"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

# load software
ml BioPerl/1.7.2-GCCcore-8.3.0
ml seqkit/0.12.1


# extract cds from gb, remove duplicate sequences, and rename headers
for i in *gb; do {
  bp_extract_feature_seq.pl -i $i -o $i.cds.fasta
  cat $i.cds.fasta | seqkit rmdup -s | sed 's/\s.*$//' > $i.fixed.cds.fasta
} done

# remove fastas with wrong headers and duplicates
rm *gb.cds.fasta

# rename fasta files
for i in *fixed.cds.fasta; do mv "$i" "${i/fixed.cds./}"; done

# run poolseq script, which takes a directory of sample fasta files (each header is a gene) and creates gene fasta files (each file is a gene and each header is a sample)
perl /home/eab77806/sarracenia-plastome-project/scripts/phylogeny_estimation/poolseq.pl

# set up directory for gene alignments
mkdir pooled
cp *fa pooled/
cd pooled/

# create genelist.txt file from the names of each gene fasta file
for i in *; do echo $i | sed 's/.fa//'; done > genelist.txt

# run script that creates scripts for each gene and aligns them with mafft
perl /home/eab77806/sarracenia-plastome-project/scripts/phylogeny_estimation/mafft.pl $OUTDIR/pooled genelist.txt

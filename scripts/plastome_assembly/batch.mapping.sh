#!/bin/bash

# create output directory if it does not exist, and change directory into it
OUTDIR="/scratch/eab77806/trimmed_reads"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

# create array of sample names; in this case all of my sample names are the first four characters of the file names
mapfile -t SAMPLE < <(for f in m*gz; do echo $f | cut -b 1-4; done | uniq)

# write and submit a script for each of the samples
for i in "${SAMPLE[@]}"
do
    {
      echo '#!/bin/bash
#SBATCH --job-name='$i'_mapping
#SBATCH --partition=batch
#SBATCH --ntasks=6
#SBATCH --mem=24gb
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/eab77806/logs/%j.out
#SBATCH --error=/scratch/eab77806/logs/%j.err

cd $SLURM_SUBMIT_DIR

# load software
ml BWA/0.7.17-GCC-8.3.0
ml SAMtools/1.10-GCC-8.3.0
ml BCFtools/1.10.2-GCC-8.3.0
ml BEDTools/2.29.2-GCC-8.3.0
ml SPAdes/3.14.1-GCC-8.3.0-Python-3.7.4
ml Fast-Plast/1.2.8-foss-2019b-Perl-5.30.0
ml SSPACE_Basic/2.1.1-foss-2019b-Perl-5.30.0

# map reads to reference
bwa mem -t 6 /home/eab77806/referenceGenomes/RRCAN01.fasta '$i'_R1_P.fastq.gz '$i'_R2_P.fastq.gz | samtools view -O BAM -o /scratch/eab77806/assemblies/'$i'.bam
cd /scratch/eab77806/assemblies
# sort bam file
samtools sort --threads 6 '$i'.bam -o '$i'.sorted.bam
# index bam file
samtools index '$i'.sorted.bam
# create bam file with only mapped reads
samtools view -b -F 4 '$i'.bam > '$i'.mapped.bam
# convert mapped bam file to fastq reads
bedtools bamtofastq -i '$i'.mapped.bam -fq '$i'.mapped_1.fq -fq2 '$i'.mapped_2.fq

# assemble mapped reads with SPAdes
spades.py -k 55,87,121 --only-assembler --careful -t 6 --memory 24 --pe1-1 '$i'.mapped_1.fq --pe1-2 '$i'.mapped_2.fq -o '$i'.mapped_spades

# extend contigs with afin
afin -c ./'$i'.mapped_spades/scaffolds.fasta -r /scratch/eab77806/trimmed_reads/'$i'_R1_P.fastq.gz /scratch/eab77806/trimmed_reads/'$i'_R2_P.fastq.gz -o '$i'_s_afin'
    } >> $i.mapping.sh
    sbatch $i.mapping.sh
    echo $i'.mapping submitted!'
done

#!/bin/bash

#create output directory
OUTDIR="/scratch/eab77806/trimmed_reads"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

#create array SAMPLE with all sample numbers
mapfile -t SAMPLE < <(for f in m*gz; do echo $f | cut -b 1-4; done | uniq)

for i in "${SAMPLE[@]}"
do
    {
      echo '#!/bin/bash'
      echo '#SBATCH --job-name='$i'_mapping'
      echo '#SBATCH --partition=batch'
      echo '#SBATCH --ntasks=6'
      echo '#SBATCH --mem=24gb'
      echo '#SBATCH --time=08:00:00'
      echo '#SBATCH --output=/scratch/eab77806/logs/%j.out'
      echo '#SBATCH --error=/scratch/eab77806/logs/%j.err'
      echo ' '
      echo 'cd $SLURM_SUBMIT_DIR'
      echo ' '
      echo '#load modules'
      echo 'ml BWA/0.7.17-GCC-8.3.0'
      echo 'ml SAMtools/1.10-GCC-8.3.0'
      echo 'ml BCFtools/1.10.2-GCC-8.3.0'
      echo 'ml BEDTools/2.29.2-GCC-8.3.0'
      echo 'ml SPAdes/3.14.1-GCC-8.3.0-Python-3.7.4'
      echo 'ml Fast-Plast/1.2.8-foss-2019b-Perl-5.30.0'
      echo 'ml SSPACE_Basic/2.1.1-foss-2019b-Perl-5.30.0'
      echo ' '
      echo '#map'
      echo 'bwa mem -t 6 /home/eab77806/referenceGenomes/RRCAN01.fasta '$i'_R1_P.fastq.gz '$i'_R2_P.fastq.gz | samtools view -O BAM -o /scratch/eab77806/assemblies.5/'$i'.bam'
      echo 'cd /scratch/eab77806/assemblies.5'
      echo 'samtools sort --threads 6 '$i'.bam -o '$i'.sorted.bam'
      echo 'samtools index '$i'.sorted.bam'
      echo 'samtools view -b -F 4 '$i'.bam > '$i'.mapped.bam'
      echo 'bedtools bamtofastq -i '$i'.mapped.bam -fq '$i'.mapped_1.fq -fq2 '$i'.mapped_2.fq'
      echo ' '
      echo '#assemble mapped reads with SPAdes'
      echo 'spades.py -k 55,87,121 --only-assembler --careful -t 6 --memory 24 --pe1-1 '$i'.mapped_1.fq --pe1-2 '$i'.mapped_2.fq -o '$i'.mapped_spades'
      echo ' '
      echo '#extend contigs with afin'
    #  echo 'afin -c ./'$i'.mapped_spades/contigs.fasta -r /scratch/eab77806/trimmed_reads/'$i'_R1_P.fastq.gz /scratch/eab77806/trimmed_reads/'$i'_R2_P.fastq.gz -o '$i'_c_afin'
      echo 'afin -c ./'$i'.mapped_spades/scaffolds.fasta -r /scratch/eab77806/trimmed_reads/'$i'_R1_P.fastq.gz /scratch/eab77806/trimmed_reads/'$i'_R2_P.fastq.gz -o '$i'_s_afin'
      # echo 'ml MUMmer/4.0.0beta2-foss-2019b'
      # echo 'ml MUMmer/3.23_conda'
      # echo ' '
      # echo '#mummer plots'
      # echo ''
      # echo 'module purge'
      # echo ''
      # echo 'ml MUMmer/4.0.0beta2-foss-2019b'
      # echo 'ml MUMmer/3.23_conda'
      # echo ''
      # # echo 'nucmer /home/eab77806/referenceGenomes/NC_041129.fasta ./'$i'.mapped_spades/scaffolds.fasta -p '$i'mapped_spades'
      # # echo 'delta-filter -1 '$i'mapped_spades.delta  >  '$i'mapped_spades.1delta'
      # # echo 'mummerplot --size large -layout --color -f --png '$i'mapped_spades.1delta -p '$i'mapped_spades'
      # echo ' '
      # echo 'nucmer /home/eab77806/referenceGenomes/NC_041129.fasta ./'$i'_afin_iter0.fa -p '$i'mapped_afin'
      # echo 'delta-filter -1 '$i'mapped_afin.delta  >  '$i'mapped_afin.1delta'
      # echo 'mummerplot --size large -layout --color -f --png '$i'mapped_afin.1delta -p '$i'mapped_afin'
      # echo ''
    } >> $i.mapping.sh
    sbatch $i.mapping.sh
    echo $i'.mapping submitted!'
done

#!/bin/bash
#SBATCH --job-name=get_organelle
#SBATCH --partition=batch
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eab77806@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=024:00:00
#SBATCH --output=/scratch/eab77806/logs/get_organelle.%j.out
#SBATCH --error=/scratch/eab77806/logs/get_organelle.%j.err


INDIR="/scratch/eab77806/trim_reads/trimmed_reads"
if [ ! -d $INDIR ]
then
    mkdir -p $INDIR
fi
cd $INDIR

ml GetOrganelle/1.7.5.2-foss-2020b

mapfile -t SAMPLE < <(for f in *gz; do echo $f | cut -b 1-4; done | uniq)

OUTDIR="/scratch/eab77806/sarracenia_plastome"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

for i in "${SAMPLE[@]}"
do get_organelle_from_reads.py -1 $INDIR/${i}_P_R1.fastq.gz -2 $INDIR/${i}_P_R2.fastq.gz  -o ${i}.plastome_output -k 21,45,65,85,97,105,121 -F embplant_pt
done


#get_organelle_from_reads.py -1 ${OUTDIR}/Spu-E-R1.fastq.gz -2 ${OUTDIR}/Spu-E-R2.fastq.gz -o ${OUTDIR}/plastome_output -s /home/eab77806/referenceGenomes/NC_041129.fasta -k 71,85,91,105 -F embplant_pt --reduce-reads-for-coverage inf --max-reads 3E7

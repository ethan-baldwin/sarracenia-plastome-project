#!/usr/bin/perl
use strict;
#Karolina Heyduk - heyduk@uga.edu - 2014

#run in directory where reads are located

my $list = $ARGV[0]; #read to ID index file, see manual for example
my %read1;
my %read2;
my @libIDs;

#read library index file, store R1 and R2 in hashes
open IN, "<$list";
while (<IN>) {
	chomp;
	my ($libID, $readID) = split /\t/;
	chop $readID;
	 if ($readID =~ '_1_') {
                $read1{$libID} = $readID;
                print "$libID\t$readID\n";
                if ($libID ~~ @libIDs) {
                        next;
                        }
                else {
                        push(@libIDs, $libID);
                        }
                }
        elsif ($readID =~ '_2_') {
                $read2{$libID} = $readID;
                print "$libID\t$readID\n";
                }
		else {
			print "Failed\n";
		}
	}
close IN;



for my $libID (@libIDs) {
    open OUT, ">$libID.trim.sh";

	#change path to trimmomatic for your usage. Alter phred score (64 or 33) as needed.
	print OUT "#!/bin/bash
#SBATCH --job-name=$libID\_trim
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=40gb
#SBATCH --time=020:00:00
#SBATCH --output=/scratch/eab77806/logs/$libID\_trim.%j.out
#SBATCH --error=/scratch/eab77806/logs/$libID\_trim.%j.err

cd \$SLURM_SUBMIT_DIR
module load Trimmomatic/0.39-Java-1.8.0_144
time java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -threads 2 -trimlog $libID\_trim.log $read1{$libID} $read2{$libID} $libID\_R1_P.fastq.gz $libID\_R1_U.fastq.gz $libID\_R2_P.fastq.gz $libID\_R2_U.fastq.gz ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
	system "sbatch $libID.trim.sh"
        }
	#change the submission line below to meet your system's needs. If threading is not available, remove the -th flag above.

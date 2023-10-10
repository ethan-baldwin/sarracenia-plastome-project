#!/usr/bin/perl
use strict;
use Cwd;
use Data::Dumper;
#Shawn Thomas - shawnt@uga.edu - 2019
#Philip Bentz adapted this for MAFFT - 2020
#USAGE: perl mafft.pl <pasta/dir> <listofgenes.txt>


my $dir = $ARGV[0]; #FULL PATH to pasta folder, do not include final "/"
my $list = $ARGV[1]; #list of genes, should same as folder names
my %control;
my $wd = getcwd;

#read gene index file, store in hash
my @genes;
open IN, "<$list";
while (<IN>) {
	chomp;
	push (@genes, $_);
	}
close IN;
print "@genes\n";

#make pasta submission scripts
for my $gene (@genes) {
#		chdir("$dir/$gene");
		open OUT, ">$gene.mafft.sh"; #make a shell file for pasta submission
        print OUT "#!/bin/bash
#SBATCH --job-name=$gene.mafft
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=80gb
#SBATCH --time=50:00:00
#SBATCH --output=/scratch/eab77806/logs/mafft.%j.out
#SBATCH --error=/scratch/eab77806/logs/mafft.%j.err
cd \$SLURM_SUBMIT_DIR
ml MAFFT/7.470-GCC-8.3.0-with-extensions

mafft --thread 12 --auto $gene\.fa > $gene\.aln";
        system "sbatch ./$gene.mafft.sh";
		chdir("$wd");
 		close OUT;
        }

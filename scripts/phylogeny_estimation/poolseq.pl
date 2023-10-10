#!/usr/bin/perl

#From your current working directory, read and store all files ending with .fa, .fas, .fna
#and .fasta, and store sequences based on sequence name. Sequences with the same name are
#then printed to .fa files, and the name of the sequence updated to include the file it
#was read from.
#AJB
# Adam Bewick

use strict;
use warnings;

#read and store all files ending with .fa, .fas, .fna and .fasta in an array
my @files = glob ("*.fa *.fas *.fna *.fasta");

my %hoh;

#loop over files
foreach my $file (@files) {
        my $id;
        my $gene;
        my ($sp, undef) = split ".f[a|as|na|asta]", $file;
        open INFILE, "<$file";
        while (<INFILE>) {
                chomp;

                if ($_ =~ /^>(.+)/) {
                        #everything after "_-_", e.g., ycf1_CDS
                        # (undef, $id) = split "_-_", $1;
                         #remove anything after "_", e.g., ycf1
                        # $gene = $id =~ s/_.*//r;
                        $gene = $1;

                }
                else {
                        # hash of hashes that separately stores sequence id and species name (key)
                        # with sequence (value) per gene name
                        $hoh{$gene}{$sp} .= $_;
                }
        }
}

foreach (keys %hoh) {
        #creates an outfile per gene name
        my $outfile = $_.".fa";
        open OUTFILE, ">$outfile" or die;
        foreach my $x (keys %{$hoh{$_}}) {
                #prints sequence id with species name and sequence in fasta format
                print OUTFILE ">$x", "\n", $hoh{$_}{$x}, "\n";
        }
}

#!/usr/bin/perl
#Author: huilong du
#Note: split the fasta file into 70 base per line
use warnings;
use strict;
my $infile=shift;              #seq fasta
my $outfile=shift;             #70 base per line

open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;
while(<IN>){
       chomp;
       my $line=$_;
       if($line=~/^>/){
              print OUT "$line\n";
       }
       else{
              $line=~s/(.{100})/$1\n/g;
              print OUT "$line\n";
       }
}
close(IN);
close(OUT);

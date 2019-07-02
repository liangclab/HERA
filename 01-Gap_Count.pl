#!/usr/bin/perl
###################################################################################################
#          this script is used to compute the gap of scaffold connected by bionano map
#                     Author: huilong du
#          Usage: perl $0 input.fasta Enzyme output.txt 
###################################################################################################

use warnings;
use strict;
my $infile=shift;
my $Enzyme=shift;
my $outfile=shift;
open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;
$Enzyme="N".$Enzyme."N";
my $chr;
while(<IN>){
     chomp;
     my $line=$_;
     if($line=~/^>(\S+)/){
          $chr=$1;
     }
     else{
	  $line=~s/$Enzyme/NNNNNNNNN/g;        #remove the single site of enzyme introduced by bionano map
          while ($line=~ m/([N,n]+)/g) { 
	            my $len = length($1); 
		    my $end = pos($line); 
		    my $start = $end - $len + 1; 
	            print OUT "$chr\t$start\t$end\t$len\n";
          }   
     }
}
close(IN);

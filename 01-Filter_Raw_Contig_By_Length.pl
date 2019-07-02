#!/usr/bin/perl
use warnings;
use strict;
my $infile=shift;                         #Total fasta
my $outfile1=shift;			  #Filtered large sequence
my $outfile2=shift;                       #Filtered small sequence
my $large_len=shift;
my $minimum_len=shift;

open IN,"<$infile" or die $!;
open OUT1,">$outfile1" or die $!;
open OUT2,">$outfile2" or die $!;
my $sign="";
while(<IN>){
	chomp;
	my $line=$_;
	if($line=~/^>(\S+)/){
		$sign=$1;
	}
	else{
		if(length($line)>=$large_len){
			print OUT1 ">$sign\n$line\n";
		}
		elsif(length($line)>=$minimum_len && length($line)<$large_len){
			print OUT2 ">$sign\n$line\n";
		}
	}
}
close(IN);
close(OUT1);
close(OUT2);

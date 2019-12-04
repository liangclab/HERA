#!/usr/bin/perl
use strict;
use warnings;

my $infile1=shift;          #Gap.txt
my $infile2=shift;          #Len.txt
open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;

my %Contig_Len=();
while(<IN2>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$Contig_Len{$content[0]}=$content[1];
}
close(IN2);
#Super-Scaffold_10-11-12/Super-Scaffold_10-11-12-Total.txt
#OUT 
#+       Super-Scaffold_10.17    280423  281991  318652  Super-Scaffold_10.18    66102   67688   126707  160
while(<IN1>){
	chomp;
	my $line=$_;
	open OUT,">>$line" or die $!;
	$line=~/(Super-Scaffold_\d+)-(\d+)-(\d+)/;
	my $scaffold=$1;
	my $First=$2;
	my $Second=$3;
	my $First_Contig=$scaffold.".".$First;
	my $Second_Contig=$scaffold.".".$Second;
	print OUT "-\t$First_Contig\t40000\t40001\t$Contig_Len{$First_Contig}\t$Second_Contig\t40000\t40001\t$Contig_Len{$Second_Contig}\t4000000\n";
	close(OUT);
}
close(IN1);
	


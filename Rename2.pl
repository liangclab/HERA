#!/usr/bin/perl
use warnings;
use strict;

my $infile1=shift;       #Pair Rename
my $infile2=shift;       #input fasta
my $outfile=shift;       #output fasta
open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;
open OUT,">$outfile" or die $!;

#SuperScaffold_9.169-SuperScaffold_9.170 Super-Scaffold_1-1-2
my %Pair_Rename=();
while(<IN1>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$Pair_Rename{$content[0]}=$content[1];
}
close(IN1);

#>SuperScaffold_7.158-SuperScaffold_7.192.2
while(<IN2>){
	chomp;
	my $line=$_;
	if($line=~/^>(\S+)\.(\d+)$/){
		my $Pair=$1;
		my $num=$2;
		my $New_Pair=$Pair_Rename{$Pair};
		print OUT ">$New_Pair.$num\n";
	}
	else{
		print OUT "$line\n";
	}
}
close(IN2);
close(OUT);

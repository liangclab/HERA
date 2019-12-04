#!/usr/bin/perl
#####################################################################################################
#      renaming the pathcontig combined with the superscaffold information
#                          Author: huilong du
#      usage:perl $0 Path2Scaffold.txt PathContig.fasta PathContig_Rename.fasta
#####################################################################################################
use warnings;
use strict;
my $infile1=shift;                 #Path2Scaffold.txt
my $infile2=shift;                 #PathContig.fasta
my $outfile=shift;                 #PathContig_Rename.fasta
open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;
open OUT,">$outfile" or die $!;

my %Path2Scaffold=();
#R498_1000597_1/0_20096_R498_1529553_1/0_15939_R498_4185055_1/534_10716_R498_1943875_1/174_16211_R498_1847476_1/10222_18999      50548   Super-Scaffold_262-5-6-Passing
while(<IN1>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$content[2]=~/(\S+)-Passing/;
	my $scaffold_info=$1;
	$Path2Scaffold{$content[0]}=$scaffold_info;
}
close(IN1);


my $count=1;
my $sign=0;
while(<IN2>){
	chomp;
	my $line=$_;
	if($line=~/^>(\S+)/){
		my $path=$1;
		if(!exists $Path2Scaffold{$path}){
			$sign=1;
			next;
		}
		$sign=0;
		print OUT ">$Path2Scaffold{$path}.$count\n";
		print "$path\t$Path2Scaffold{$path}.$count\n";
		$count++;
	}
	else{
		next if($sign==1);
		print OUT "$line\n";
	}
}
close(IN2);
close(OUT);

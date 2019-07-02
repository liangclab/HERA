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
#R498_1000511_1/584_7373_R498_3566492_1/2084_16056_R498_1347669_1/235_18794      20610   ctg7180000005897-Tail-ctg7180000009055-Tail-BestMatching-Passing
while(<IN1>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	my @Contig_Pair=split /-/,$content[2];
	my $scaffold_info="";
	if($Contig_Pair[0] lt $Contig_Pair[2]){
	        $scaffold_info=$Contig_Pair[0]."-".$Contig_Pair[2];
	}
	else{
		$scaffold_info=$Contig_Pair[2]."-".$Contig_Pair[0];
	}
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

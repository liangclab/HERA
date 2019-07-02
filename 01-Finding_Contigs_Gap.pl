#!/usr/bin/perl
########################################################################
#     make the both contigs of the gap into pairs for positions
#                 Author: huilong du
########################################################################
use warnings;
use strict;
my $infile=shift;
my $outfile=shift;
open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;

my $Current_Scaffold="";
my $count=1;
#Super-Scaffold_3        364273  364771  499
while(<IN>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	if($Current_Scaffold ne $content[0]){
		$count=1;
		$Current_Scaffold=$content[0];
		my $Ctg1=$content[0].".".$count;
		$count++;
		my $Ctg2=$content[0].".".$count;
		print OUT "$Ctg1\t$Ctg2\t$content[3]\t$content[1]\t$content[2]\n";
	}
	else{
		my $Ctg1=$content[0].".".$count;
		$count++;
		my $Ctg2=$content[0].".".$count;
		print OUT "$Ctg1\t$Ctg2\t$content[3]\t$content[1]\t$content[2]\n";
	}
}
close(IN);
close(OUT);	

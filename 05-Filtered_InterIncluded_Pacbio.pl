#!/usr/bin/perl
#Author: huilong du
#Note: Remove the pacbios and non-scaffolded contigs aligned to the internal scaffolded contigs
use warnings;
use strict;
my $infile=shift;                      #reads to ref
my $outfile=shift;	               #Included Pacbio
my $MinIdentity=shift;		       #Min Identity
my $MinCoverage=shift; 	               #Min Coverage
my $MinExtend=shift;		       #Min Ref Both Side
open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;
my %Existed_Pacbio=();
while(<IN>){
	chomp;
	my $line=$_;
	next if($line=~/^qName/);
	my @content=split /\s+/,$line;
	$content[0]=~/(\S+)\/0_\d+$/;
	my $Pacbio=$1;
	next if(exists $Existed_Pacbio{$Pacbio});
	next if($content[5]<$MinIdentity);
	my $coverage=($content[10]-$content[9])/$content[11]*100;
	next if($coverage<$MinCoverage);
	my $Start=$MinExtend;
	my $End=$content[8]-$MinExtend;
	if($content[6]>=$Start && $content[7]<=$End){
		$Existed_Pacbio{$Pacbio}=0;
		print OUT "$Pacbio\n";
	}
}
close(IN);
close(OUT);

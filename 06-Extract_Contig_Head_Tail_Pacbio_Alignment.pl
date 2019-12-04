#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $infile;
my $MinIdentity=97;
my $MaxOverhang=50;
my $MinCoverage=98;
my $HTLen=30000;
my $MinLen=7000;
GetOptions( "Align=s" => \$infile,
             "MinIden=i" => \$MinIdentity,
             "MaxHang=i"  => \$MaxOverhang,
             "MinCov=i" => \$MinCoverage,
	     "HTLen=i" => \$HTLen,
	     "MinLen=i" => \$MinLen
);

if( !defined $infile){
	print "
Usage: perl $0 -Align=file -MinIden=num1 -MaxHang=num2 -MinCov=num3 -HTLen=num4 -MinLen=num5

 [ -Align ]       The File of Alignment
 [ -MinIden ]     The MinIdentity of Alignment
 [ -MaxHang ]     The Max Overhang of Query Alignment
 [ -MinCov ]      The Min Coverage of Alignment
 [ -HTLen ]       The Head and Tail Length Of Ref
 [ -MinLen ]      The Min Length of Query\n\n";
	exit;
}
open IN,"<$infile" or die $!;
open OUT,">Contig_Head_Tail_Pacbio.txt" or die $!;
while(<IN>){
	chomp;
	my $line=$_;
	next if($line=~/qLength/);
	my @content=split /\s+/,$line;
	next if($content[11]<$MinLen);
	next if($content[5]<$MinIdentity);
	my $Overhang=$content[11]-($content[10]-$content[9]);
	my $Coverage=($content[10]-$content[9])/$content[11]*100;
	if(($content[8]-$content[7]<$MaxOverhang && $content[7]-$content[6]>2000 && $content[9]<$MaxOverhang)||($content[6]<$MaxOverhang && $content[7]-$content[6]>2000 && $content[11]-$content[10]<$MaxOverhang)){
		print OUT "$line\n";
		next;
	}
	next if($Overhang>$MaxOverhang || $Coverage<$MinCoverage);
	next if($content[6]>$HTLen && $content[8]-$content[7]>$HTLen);
	print OUT "$line\n";
}
close(IN);

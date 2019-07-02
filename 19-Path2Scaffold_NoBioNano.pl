#!/usr/bin/perl
use warnings;
use strict;
my $infile=shift;                         #Selected_Pairs
my $infile1=shift;			  #ctg clusters
my $infile2=shift;			  #Path length
my $outfile=shift;			  #Path--Scaffold

open IN,"<$infile" or die $!;
open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;
open OUT,">$outfile" or die $!;

my %Selected_Contig_Pairs=();
#ctg7180000006959        Head    ctg7180000007221        Head    1       4       0       5       0       5
while(<IN>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$Selected_Contig_Pairs{$content[0]}{$content[1]}{$content[2]}{$content[3]}=0;
	$Selected_Contig_Pairs{$content[2]}{$content[3]}{$content[0]}{$content[1]}=0;
}
close(IN);

my %PathInScaffold=();
#R498_1000511_1/584_7373 R498_1167005_1/1_23605  ctg7180000005803-Tail-ctg7180000005838-Head-MostExtending-Passing
while(<IN1>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	my @Pair=split /-/,$content[-1];
	next if(!exists $Selected_Contig_Pairs{$Pair[0]}{$Pair[1]}{$Pair[2]}{$Pair[3]} || !exists $Selected_Contig_Pairs{$Pair[2]}{$Pair[3]}{$Pair[0]}{$Pair[1]});
	my $path=$content[0];
	for(my $i=1;$i<@content-1;$i++){
		$path=$path."_".$content[$i];
	}
	$PathInScaffold{$path}=$content[-1];
}
close(IN1);

my %CommonPath=();
#R498_1000597_1/0_20096_R498_1529553_1/0_15939_R498_4703259_1/3_25063_R498_551312_1/0_22055_R498_1943875_1/174_16211_R498_2980934_1/35_16991	57580
while(<IN2>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	next if(exists $CommonPath{$content[0]});
	if(exists $PathInScaffold{$content[0]}){
		print OUT "$line\t$PathInScaffold{$content[0]}\n";
		$CommonPath{$content[0]}=0;
	}
}
close(IN2);

=pod
foreach my $key (keys %PathInScaffold){
	next if(exists $CommonPath{$key});
	print "$key\t$PathInScaffold{$key}\n";
}
=cut

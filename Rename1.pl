#!/usr/bin/perl
use warnings;
use strict;
my $infile=shift;
my $outfile=shift;

open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;

#Super-Scaffold_1
my $count=1;
my $sign="";
open OUT1,">rename.pbs" or die $!;
print OUT1 "#PBS -N Rename
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q hldu
cd \$PBS_O_WORKDIR\n";
while(<IN>){
	chomp;
	my $line=$_;
	if($line=~/^>(\S+)/){
		my $path=$1;
		my $temp=1;
		my @content=split /--/,$path;
		print OUT ">Super-Scaffold_$count\n";
		$temp++;
		for(my $i=1;$i<@content;$i++){
			my $scaffold="";
			if($content[$i-1] lt $content[$i]){
				$scaffold=$content[$i-1]."-".$content[$i];
			}
			else{
				$scaffold=$content[$i]."-".$content[$i-1];
			}
			my $First=$temp-1;
			my $Second=$temp;
			print OUT1 "sed -i \'s/^>$scaffold/>Super-Scaffold_$count-$First-$Second/g\' ./PathContig_Rename.fasta\n";
			print "$scaffold\tSuper-Scaffold_$count-$First-$Second\n";
			$temp++;
		}
		$count++;
	}
	else{
		print OUT "$line\n";
	}
}
close(IN);
close(OUT);
close(OUT1);
#my $commond=`qsub rename.pbs`;

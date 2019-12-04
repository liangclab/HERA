#!/usr/bin/perl
#Author: huilong du
#Note: align the corrected pacbios and non-scaffolded contigs selfly
use warnings;
use strict;

my $infile=shift;
my $script=shift;
my $queue=shift;
my $genome=shift;
open IN,"<$infile" or die $!;


my $count=1;
my @part=();
while(<IN>){
	chomp;
	my $line=$_;
	push(@part,$line);
}
close(IN);
for(my $i=0;$i<@part;$i++){
    for(my $j=$i;$j<@part;$j++){
	open OUT,">$count.pbs" or die $!;
	if($count<350){
		print OUT "#BSUB -J $genome-Pair-$i-$j
#BSUB -o $count.out
#BSUB -n 2
#BSUB -q $queue
";
        }
        elsif($count>=350){
                print OUT "#BSUB -J $genome-Pair-$i-$j
#BSUB -o $count.out
#BSUB -n 2
#BSUB -q $queue
";
        }
	if($i==$j){
           print OUT " 
bwa mem -a -e -t 4 $part[$i] $part[$j] >./03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam
perl $script/sam2blasr.pl ./03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam ./03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.txt
rm -f ./03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam
";
	}
	else{
	   print OUT "
bwa mem -a -t 4 $part[$i] $part[$j] >./03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam
perl $script/sam2blasr.pl ./03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam ./03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.txt
rm -f ./03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam
";
	}
        close(OUT);
        system("bsub < $count.pbs");
        $count++;
    }
}


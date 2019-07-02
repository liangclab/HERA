#!/usr/bin/perl
use warnings;
use strict;

my $infile=shift;
my $queue=shift;
my $genome=shift;
open IN,"<$infile" or die $!;

my $count=1;

while(<IN>){
    chomp;
    my $line=$_;
    open OUT,">$count.pbs" or die $!;
    if($count<120){
        print OUT "#BSUB -J $genome-INDEX-$count
#BSUB -o $count.out
#BSUB -n 1
#BSUB -q $queue
";
     }
     elsif($count>=120){
        print OUT "#BSUB -J $genome-INDEX-$count
#BSUB -o $count.out
#BSUB -n 1
#BSUB -q $queue
";
     }
print OUT " 
bwa index $line
";
     close(OUT);
     system("bsub < $count.pbs");
     $count++;
}
close(IN);


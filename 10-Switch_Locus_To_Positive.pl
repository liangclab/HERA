#!/usr/bin/perl
#Author: huilong du
#Note: change the aligned positions into the positive chain
use warnings;
use strict;

my $infile=shift;
my $outfile=shift;
open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;


my %Chr_Order=();
my $start=0;
my $end=0;
while(<IN>){
    chomp;
    my $line=$_;
    next if($line=~/^qName/);
    my @content=split /\s+/,$line;
    if($content[3]==0){
           my $ref_end=$content[7]-1;
           my $query_start=$content[9]+1;
           print OUT "$content[0] $query_start $content[10] $content[11] $content[1] $content[6] $ref_end $content[8] 0 $content[5]\n";
    }
    elsif($content[3]==1){
           my $Ref_start=$content[8]-$content[7];
           my $Ref_end=$content[8]-$content[6]-1;
           my $query_start=$content[9]+1;
           print OUT "$content[0] $query_start $content[10] $content[11] $content[1] $Ref_start $Ref_end $content[8] 1 $content[5]\n";
    }           

}
close(IN);
close(OUT);


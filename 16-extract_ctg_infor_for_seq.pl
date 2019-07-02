#!/usr/bin/perl
#Author: huilong du
#Note: extract the contig information for extracting the seq
use warnings;
use strict;
my $infile=shift;                  #cluster_pos_line.txt
open IN,"<$infile" or die $!;
my $outfile=shift;                 #cluster_pos_line_for_seq.txt
open OUT,">$outfile" or die $!;

while(<IN>){
     chomp;
     my $line=$_;
     $line=~s/\s+$//;
     my @lines=split /\s+/,$line;
     my $ctg_count=0;
     for(my $j=0;$j<@lines;$j=$j+10){
            $ctg_count++;
     }
     my $ctg_num=$ctg_count;
     my @content=split /\s+/,$line;
     my $count=1;
     for(my $i=0;$i<@content;$i=$i+10){
             if($count==1){
                    print OUT "$content[$i] $content[$i+1] $content[$i+2] $content[$i+3] $content[$i+4] $content[$i+5] $content[$i+6] $content[$i+7] ";
                    $count++;
             }
             elsif($count==$ctg_num){
                    print OUT "$content[$i] $content[$i+1] $content[$i-2] $content[$i-1]\n";
                    $count++;
             }
             else{
                    print OUT "$content[$i] $content[$i+1] $content[$i-2] $content[$i+3] $content[$i+4] $content[$i+5] $content[$i+6] $content[$i+7] ";
                    $count++;
             }
      }
}
close(IN);
close(OUT);
             
              

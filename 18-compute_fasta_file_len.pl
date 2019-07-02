#!/usr/bin/perl
#Author: huilong du  
#Note: compute the fasta or fastq len

use strict;
use warnings;
my $infile=shift;
my $outfile=shift;
open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;

my $sign;
my $count=1;
my $total_len=0;
my $len=0;
#if($infile=~/fasta$/ or $infile=~/fa$/){
if(1){
        while(<IN>){
              chomp;
              my $line=$_;
              if($line=~/^>(\S+)/){
                     if($count==1){
                           $sign=$1;
                           $count++;
                     }
                     else{
                           print OUT "$sign\t$len\n"; 
#                           $total_len=$total_len+$len;
                           $sign=$1;
                           $count++;
                           $len=0;
                     }
                     
              }
              else{
                     $len=$len+length($line);
                     $total_len=$total_len+$len;
              }
        }
        print OUT "$sign\t$len\n";
}
elsif($infile=~/fastq$/){
         while(<IN>){
              chomp;
              my $line=$_;
              if($count%4==1){
                     $line=~/^@(\S+)/;
                     $sign=$1;
              }
              elsif($count%4==2){
                     my $len=length($line);
                     print OUT "$sign\t$len\n";
              }
              $count++;
         }
}
print "the total length of genome is $total_len\n";
close(IN);
close(OUT);
              

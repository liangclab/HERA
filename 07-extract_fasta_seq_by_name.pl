#!/usr/bin/perl
#Author: huilong du
#Note: extract the corrsponding seq by read_name
use warnings;
use strict;
my $infile1=shift;                            #seq name file
my $infile2=shift;                            #original fasta file
my $outfile=shift;                            #fasta seq
open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;
open OUT,">$outfile" or die $!;

my %existed_seq;
while(<IN1>){
    chomp;
    my $line=$_;
    my @content=split /\s+/,$line;
    $existed_seq{$content[0]}=0;
}
close(IN1);


my $sign=0;
while(<IN2>){
    chomp;
    my $line=$_;
    if($line=~/^>(\S+)/){
          my $ctg=$1;
          if(!exists $existed_seq{$ctg}){
                     print OUT "$line\n";
                     $sign=1;
          }
          else{  
                     $sign=0;
          }
    }
    else{
          if($sign==1){
                     print OUT "$line\n";
          }
    }
}
close(IN2);
close(OUT);

      
    

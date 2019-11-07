#!/usr/bin/perl
#Author: huilong du
#Note: 利用比对上的长度和Identity过滤重复的CTG
use warnings;
use strict;
use List::Util qw/max min sum maxstr minstr shuffle/;
my $infile=shift;     #bwa align sam
my $outfile=shift;    #non_repeat fasta
open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $1;


my %existed_ctg=();
my %existed_pairs=();
my %query_len=();
my %ref_len=();
print OUT "qName tName qStrand tStrand score percentSimilarity tStart tEnd tLength qStart qEnd qLength nCells\n";
while(<IN>){
     chomp;
     my $line=$_;
     if($line=~/^\@SQ/){
          $line=~/SN:(\S+)\tLN:(\d+)$/;
          $ref_len{$1}=$2;
          next;
     }
     next if($line=~/^@/);
     my @content=split "\t",$line;
#     next if($content[0] eq $content[2]);
     next if($content[2] eq "*" || $content[5] eq "*");
     my $left_H=0;
     my $right_H=0;
     if($content[5]=~/^(\d+)[H,S]/){
           $left_H=$1;
     }
     if($content[5]=~/(\d+)[H,S]$/){
           $right_H=$1;
     }
     my @match=$content[5]=~/(\d+)M/g;
     my $all_insert=0;
     if($content[5]=~/\d+I/){
          my @insert=$content[5]=~/(\d+)I/g;
          $all_insert=sum(@insert);
     }
     my $all_delete=0;
     if($content[5]=~/\d+D/){
          my @delete=$content[5]=~/(\d+)D/g;
          $all_delete=sum(@delete);
     }
     my $all_match=sum(@match);
     my $all_len=$all_match+$all_insert+$all_delete;
     $content[11]=~/NM:i:(\d+)/;
     my $mismatch=$1;
     my $match_len=$all_len-$mismatch;
     my $Identity=($match_len/$all_len)*100;
     $Identity=sprintf "%0.2f",$Identity;
     my $flag=($content[1] & 16)/16;
     my $tStart=0;
     my $tEnd=0;
     my $tLength=0;
     my $qStart=0;
     my $qEnd=0;
     my $qLength=0;
     if($flag==1){
          $tEnd=$ref_len{$content[2]}-$content[3]; 
          $tStart=$ref_len{$content[2]}-$content[3]-$all_match-$all_delete;
          $qStart=$right_H;
          $qEnd=$right_H+$all_match+$all_insert;
          $qLength=$left_H+$right_H+$all_match+$all_insert;
     }
     else{
          $tStart=$content[3];
          $tEnd=$content[3]+$all_match+$all_delete;
          $qStart=$left_H;
          $qEnd=$left_H+$all_match+$all_insert;
          $qLength=$left_H+$right_H+$all_match+$all_insert;
     }
     my $score=-int($Identity*($tEnd-$tStart));
     if($tStart<0){
           $tStart=0;
     }
     if($qStart<0){
           $qStart=0;
     }
     print OUT "$content[0]/0_$qLength $content[2] 0 $flag $score $Identity $tStart $tEnd $ref_len{$content[2]} $qStart $qEnd $qLength null\n";

#     next if(exists $existed_ctg{$content[0]});
#     next if(exists $existed_pairs{$content[0]}{$content[2]});
#     print OUT "$content[0]\n";
#     $existed_ctg{$content[0]}=0;
#     $existed_pairs{$content[0]}{$content[2]}=0;
#     $existed_pairs{$content[2]}{$content[0]}=0;
}
close(IN);
close(OUT);


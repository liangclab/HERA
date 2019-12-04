#!/usr/bin/perl
#Author: huilong du
#Note: make the corrected pacbios and non-scaffolded contigs into a line
use warnings;
use strict;

my $infile=shift;           #cluster_final file   test1.txt
my $outfile=shift;          #clusters of contigs left right and the same chain
open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;

while(<IN>){
     chomp;
     my $line=$_;
     my @content=split " ",$line;
     my $bac_reverse=0;
     my $current_ori="";
     my $current_chain=0;
     if($content[2] eq "left"){
          if($content[1]==1){
               print OUT "$content[0] 1 0 right ";
               $bac_reverse=0;
          }
          elsif($content[1]==0){
               print OUT "$content[0] 1 1 right ";
               $bac_reverse=1;
          }
     }
     elsif($content[2] eq "right"){
          if($content[1]==0){
               print OUT "$content[0] 0 0 right ";
               $bac_reverse=0;
          }
          elsif($content[1]==1){
               print OUT "$content[0] 0 1 right ";
               $bac_reverse=1;
          }
     }
     for(my $i=3;$i<@content;$i=$i+5){
          if($content[$i+1] ne "left"){
               if($content[$i]==1){ 
                   if($i+5<=@content){
                       if($content[$i+3]==1){
                           $content[$i+3]=0;
                       }
                       else{
                           $content[$i+3]=1;
                       }
                       if($content[$i+4] eq "left"){
                           $content[$i+4]="right";
                       }
                       elsif($content[$i+4] eq "right"){
                           $content[$i+4]="left";
                       }
                       print OUT "0 left $content[$i+2] 1 $content[$i+3] $content[$i+4] ";
                   }
                   elsif($i+5>@content){
                         print OUT "0 left $content[$i+2] 1";
                   }
               }
               else{
                   if($i+5<=@content){
                        if($content[$i+3]==1){
                            $content[$i+3]=0;
                        }
                        else{
                            $content[$i+3]=1;
                        }
                        if($content[$i+4] eq "left"){
                            $content[$i+4]="right";
                        }
                        elsif($content[$i+4] eq "right"){
                            $content[$i+4]="left";
                        }
                        print OUT "1 left $content[$i+2] 1 $content[$i+3] $content[$i+4] ";
                   }
                   else{
                         print OUT "1 left $content[$i+2] 1";
                   }
               }
          }
          elsif($content[$i+1] eq "left"){
               if($i+5<=@content){
                    print OUT "$content[$i] $content[$i+1] $content[$i+2] 0 $content[$i+3] $content[$i+4] ";
               }
               else{
                    print OUT "$content[$i] $content[$i+1] $content[$i+2] 0";
               }
          }
     }
     print OUT "\n";
}
close(IN);
close(OUT);

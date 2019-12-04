#!/usr/bin/perl
#Author: huilong du
#Note: connect the contigs into scaffolds 
no warnings;
use strict;

my $infile1=shift;                     #ctg_pair.txt
my $infile2=shift;                     #ctg_ctg_ori.txt
my $infile3=shift;                     #ctg_clusters_final.txt
#my $infile4=shift;                     #../../192x192_new_data/sspace_new_genome/side_assembly/filtered_non_rice/rice_ctg_line.fasta
#my $infile5=shift;                     #bac_seq
my $outfile=shift;                     #connected contig fasta   scaffolds_by_bac.fasta

open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;
open IN3,"<$infile3" or die $!;
#open IN4,"<$infile4" or die $!;
#open IN5,"<$infile5" or die $!;
open OUT,">$outfile" or die $!;

my %ctg_ctg_line=();
#my %ctg_ctg_overlap_len=();
while(<IN1>){
      chomp;
      my $line=$_;
      my @content=split "[\t,-]",$line;
      $ctg_ctg_line{$content[0]}{$content[2]}=$content[5]."-".$content[6]."-".$content[7]."-".$content[8]."-".$content[9]."-".$content[10]."-".$content[11]."-".$content[12]."-".$content[13]."-".$content[14]."-".$content[15]."-".$content[16]."-".$content[17];
      $ctg_ctg_line{$content[2]}{$content[0]}=$content[5]."-".$content[9]."-".$content[10]."-".$content[11]."-".$content[6]."-".$content[7]."-".$content[8]."-".$content[12]."-".$content[13]."-".$content[16]."-".$content[17]."-".$content[14]."-".$content[15];
#     $ctg_ctg_overlap_len{$content[0]}{$content[2]}=$content[8]-$content[7]+$content[11]-$content[10];
}
close(IN1);

#my %ctg_ctg_ori=();
my %ctg_len=();
#my %ctg_ctg_chain=();
#my %ctg_ctg_count=();
while(<IN2>){
      chomp;
      my $line=$_;
      my @content=split /\s+/,$line;
      $ctg_len{$content[0]}=$content[1];
#      $ctg_ctg_chain{$content[0]}{$content[4]}{1}=$content[2];
#      $ctg_ctg_chain{$content[0]}{$content[4]}{2}=$content[6];

#      $ctg_ctg_chain{$content[4]}{$content[0]}{1}=$content[6];
#      $ctg_ctg_chain{$content[4]}{$content[0]}{2}=$content[2];
 
#      $ctg_ctg_ori{$content[0]}{$content[4]}{1}=$content[3];
#      $ctg_ctg_ori{$content[0]}{$content[4]}{2}=$content[7];
      
#      $ctg_ctg_ori{$content[4]}{$content[0]}{1}=$content[7];
#      $ctg_ctg_ori{$content[4]}{$content[0]}{2}=$content[3];
#      $ctg_ctg_count{$content[0]}{$content[4]}=$content[8];
#      $ctg_ctg_count{$content[4]}{$content[0]}=$content[8];
}
close(IN2);

=pod
my %ctg_seq=();
my $sign="";
while(<IN4>){
      chomp;
      my $line=$_;
      if($line=~/^>(\S+)/){
               $sign=$1;
      }
      else{
               $ctg_seq{$sign}=$line;
      }
}
close(IN4);

my %bac_ctg=();
$sign="";
while(<IN5>){
      chomp;
      my $line=$_;
      if($line=~/^>(\S+)/){
                $sign=$1;
      }
      else{
                $bac_seq{$sign}=$line;
      }
}
close(IN5);
=cut

#ctg7180000008354 0 0 right 0 left ctg7180000008355 0 0 right 0 left ctg7180000008356 0
#into
#(Ctg1 Chain Start End ) (Bac1 Chain Start End) (Start Length Ctg2 Chain nonsense End) (Bac2 Chain Start End) (Start Length Ctg3 Chain Nonsense)
#(tig00001128 1 0 419678) (BAC5600|size153262/0_153262 0 63426 51228) (1 318643 tig00001132 0 0 318644) (BAC342|size252584/0_252584 0 66885 57303) (1 1003677 tig00001134 0 0)
while(<IN3>){
      chomp;
      my $line=$_;
      my @content=split /\s+/,$line;
      for(my $i=0;$i<@content;$i=$i+6){
	  if(!exists $ctg_ctg_line{$content[$i]}{$content[$i+6]}){
		print OUT "$content[$i] $content[$i+1] 0";
		next;
	  }
          my @infor_ctg=split "-",$ctg_ctg_line{$content[$i]}{$content[$i+6]};
          if($content[$i+1]==0){       #contig no reverse
               if($infor_ctg[1]==0){   #
                  print OUT "$content[$i] $content[$i+1] 0 $infor_ctg[2] $infor_ctg[0] $content[$i+2] ";
                  if($content[$i+2]==0){ #bac no reverse 
                         print OUT "$infor_ctg[9] ";
                  }
                  elsif($content[$i+2]==1){ #bac reverse
                         my $bac_start=$infor_ctg[8]-$infor_ctg[10];    #reverse chain and change the positions
                         
                         print OUT "$bac_start ";
                  }
               }
               elsif($infor_ctg[1]==1){
                  my $ctg_end=$ctg_len{$content[$i]}-$infor_ctg[3];
                  print OUT "$content[$i] $content[$i+1] 0 $ctg_end $infor_ctg[0] $content[$i+2] ";
                  if($content[$i+2]==0){
                         print OUT "$infor_ctg[9] ";
                  }
                  elsif($content[$i+2]==1){
                         my $bac_start=$infor_ctg[8]-$infor_ctg[10];
                         print OUT "$bac_start ";
                  }
               }
          }
          elsif($content[$i+1]==1){   #contig reverse
               if($infor_ctg[1]==0){
                    my $ctg_end=$ctg_len{$content[$i]}-$infor_ctg[3];
                    print OUT "$content[$i] $content[$i+1] 0 $ctg_end $infor_ctg[0] $content[$i+2] ";
                    if($content[$i+2]==0){
                          print OUT "$infor_ctg[9] ";
                    }
                    elsif($content[$i+2]==1){
                          my $bac_start=$infor_ctg[8]-$infor_ctg[10];
                          print OUT "$bac_start ";
                    }
                }
                elsif($infor_ctg[1]==1){
                    print OUT "$content[$i] $content[$i+1] 0 $infor_ctg[2] $infor_ctg[0] $content[$i+2] ";
                    if($content[$i+2]==0){
                          print OUT "$infor_ctg[9] ";
                    }
                    elsif($content[$i+2]==1){
                          my $bac_start=$infor_ctg[8]-$infor_ctg[10];
                          print OUT "$bac_start ";
                    }
                }
          } 
          ################################ 左边contig以及bac起始坐标确定结束 ######################################
          if($content[$i+7]==0){  #contig no reverse
                if($infor_ctg[4]==0){
                     if($content[$i+4]==0){  #bac no reverse
                         print OUT "$infor_ctg[11] $infor_ctg[5] $ctg_len{$content[$i+6]} ";
                     }
                     elsif($content[$i+4]==1){ #bac reverse
                         my $bac_end=$infor_ctg[8]-$infor_ctg[12];
                         print OUT "$bac_end $infor_ctg[5] $ctg_len{$content[$i+6]} ";
                     }
                 }
                 elsif($infor_ctg[4]==1){   
                     if($content[$i+4]==0){  #bac no reverse
                         my $ctg_start=$ctg_len{$content[$i+6]}-$infor_ctg[6];
                         if($ctg_start<0){
                              $ctg_start=0;
                         }
                         print OUT "$infor_ctg[11] $ctg_start $ctg_len{$content[$i+6]} ";
                     }
                     elsif($content[$i+4]==1){ #bac reverse
                         my $bac_end=$infor_ctg[8]-$infor_ctg[12];
                         my $ctg_start=$ctg_len{$content[$i+6]}-$infor_ctg[6];  
                         if($ctg_start<0){
                               $ctg_start=0;
                         }
                         print OUT "$bac_end $ctg_start $ctg_len{$content[$i+6]} ";
                     }
                 }
          }
          elsif($content[$i+7]==1){ #contig reverse
                 if($infor_ctg[4]==0){ #no contradiction
                     if($content[$i+4]==0){ #bac no reverse
                         my $ctg_start=$ctg_len{$content[$i+6]}-$infor_ctg[6];
                         if($ctg_start<0){
                              $ctg_start=0;
                         }
                         print OUT "$infor_ctg[11] $ctg_start $ctg_len{$content[$i+6]} ";
                     }
                     elsif($content[$i+4]==1){ #bac reverse
                         my $bac_end=$infor_ctg[8]-$infor_ctg[12];
                         my $ctg_start=$ctg_len{$content[$i+6]}-$infor_ctg[6];
                         if($ctg_start){
                                $ctg_start=0;
                         }
                         print OUT "$bac_end $ctg_start $ctg_len{$content[$i+6]} ";           
                     }
                 }
                 elsif($infor_ctg[4]==1){ 
                     if($content[$i+4]==0){  #bac no reverse
                         print OUT "$infor_ctg[11] $infor_ctg[5] $ctg_len{$content[$i+6]} ";
                     }
                     elsif($content[$i+4]==1){ #bac reverse
                         my $bac_end=$infor_ctg[8]-$infor_ctg[12];
                         print OUT "$bac_end $infor_ctg[5] $ctg_len{$content[$i+6]} ";
                     }
                 }
           }
      }    
      print OUT "\n"; 
}
close(OUT);     
                

 
                
        

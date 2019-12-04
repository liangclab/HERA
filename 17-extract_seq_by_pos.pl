#!/usr/bin/perl
#Author: huilong du
#Note: connect the contigs into scaffolds 
use warnings;
use strict;

sub reverse_complement{
        my $input_seq=shift;
        my $reverse=reverse $input_seq;
        $reverse=~tr/ACGT/TGCA/;
        $reverse=~tr/acgt/tgca/;
        return $reverse;
}
my $infile1=shift;                      #cluster_pos_line_for_seq.txt
my $infile2=shift;                     #../../192x192_new_data/sspace_new_genome/side_assembly/filtered_non_rice/rice_ctg_line.fasta
my $infile3=shift;                     #bac_seq
my $outfile=shift;                     #connected contig fasta   scaffolds_by_bac.fasta

open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;
open IN3,"<$infile3" or die $!;
open OUT,">$outfile" or die $!;

my %ctg_seq=();
my $sign="";
my %ctg_len=();
while(<IN2>){
      chomp;
      my $line=$_;
      if($line=~/^>(\S+)/){
               $sign=$1;
      }
      else{
               $ctg_seq{$sign}=$line;
               my $len=length($line);
               $ctg_len{$sign}=$len;
      }
}

my %bac_seq=();
my %bac_len=();
$sign="";
while(<IN3>){
      chomp;
      my $line=$_;
      if($line=~/^>(\S+)/){
                $sign=$1;
      }
      else{
                $bac_seq{$sign}=$line;
                my $len=length($line);
                $bac_len{$sign}=$len;
      }
}
my $count_ctg=int (rand(10000));
$count_ctg=$count_ctg*100000+1;
my %ctg_existed=();
my $Ctg_Name="";

while(<IN1>){
      chomp;
      my $line=$_;
      $Ctg_Name="";
      my @content=split " ",$line;
      my $seq="";
      my $extract_seq="";
      my $last_ctg_start=0;
      my %ctg_per_cluster=();
      for(my $long=0;$long<@content;$long=$long+8){
         if($long==0){
               $Ctg_Name=$content[$long];
         }
         else{
               $Ctg_Name=$Ctg_Name."_".$content[$long];
         }
      }
      for(my $i=4;$i<@content;$i=$i+8){
          my @info=split "/",$content[$i];
          pop (@info);
          $content[$i]=join("/",@info);
          $ctg_existed{$content[$i-4]}=0;
          $ctg_existed{$content[$i+4]}=0;
          $ctg_per_cluster{$content[$i-4]}=0;
          $ctg_per_cluster{$content[$i+4]}=0;
          if($content[$i+2]<=$content[$i+3]){  #with gap
               if($content[$i-3]==0){    #left ctg no reverse
                      if($content[$i-1]>=$ctg_len{$content[$i-4]}){
                               $content[$i-1]=$ctg_len{$content[$i-4]}-1;
                      }
                      my $extract_len=$content[$i-1]-$content[$i-2]+1;
                      $extract_seq=substr($ctg_seq{$content[$i-4]},$content[$i-2],$extract_len);
                      $seq=$seq.$extract_seq;
               }
               elsif($content[$i-3]==1){ #left ctg reverse
                      if($content[$i-1]>=$ctg_len{$content[$i-4]}){
                               $content[$i-1]=$ctg_len{$content[$i-4]}-1;     
                      }
#                      $ctg_seq{$content[$i-4]}=reverse_complement($ctg_seq{$content[$i-4]});
                      my $Current_Bac=reverse_complement($ctg_seq{$content[$i-4]});
                      my $extract_len=$content[$i-1]-$content[$i-2]+1;
                      $extract_seq=substr($Current_Bac,$content[$i-2],$extract_len);
                      $seq=$seq.$extract_seq; 
              }
              if($content[$i+1]==0){  #bac no reverse
                      if($content[$i+3]>=$bac_len{$content[$i]}){
                               $content[$i+3]=$bac_len{$content[$i]}-1;
                      }
                      my $extract_len=$content[$i+3]-$content[$i+2]+1;
                      $extract_seq=substr($bac_seq{$content[$i]},$content[$i+2],$extract_len);
                      $seq=$seq.$extract_seq;
		      print "$content[$i] 0\n";
              }
              elsif($content[$i+1]==1){  #bac reverse
                      if($content[$i+3]>=$bac_len{$content[$i]}){
                               $content[$i+3]=$bac_len{$content[$i]}-1;
                      }     
                      my $extract_len=$content[$i+3]-$content[$i+2]+1;
#                      $bac_seq{$content[$i]}=reverse_complement($bac_seq{$content[$i]});
                      my $Current_Bac=reverse_complement($bac_seq{$content[$i]});
                      $extract_seq=substr($Current_Bac,$content[$i+2],$extract_len);
		
                      $seq=$seq.$extract_seq;
#		      print ">$content[$i]\n";
#                      print "$extract_seq\n";
              }
              if($content[$i+5]==0){    #right ctg no reverse
                      if($content[$i+7]>=$ctg_len{$content[$i+4]}){
                               $content[$i+7]=$ctg_len{$content[$i+4]}-1;
                      }
                      $last_ctg_start=$content[$i+6];
=pod
                      my $extract_len=$content[$i+7]-$content[$i+6]+1;
                      $extract_seq=substr($ctg_seq{$content[$i+4]},$content[$i+6],$extract_len);
                      $seq=$seq.$extract_seq;
=cut
               }
               elsif($content[$i+5]==1){ #right ctg reverse
                      if($content[$i+7]>=$ctg_len{$content[$i+4]}){
                               $content[$i+7]=$ctg_len{$content[$i+4]}-1;
                      }
                      $last_ctg_start=$content[$i+6];
=pod
                      $ctg_seq{$content[$i+4]}=reverse_complement($ctg_seq{$content[$i+4]});
                      my $extract_len=$content[$i+7]-$content[$i+6]+1;
                      $extract_seq=substr($ctg_seq{$content[$i+4]},$content[$i+6],$extract_len);
                      $seq=$seq.$extract_seq;
=cut
              }
          }
          elsif($content[$i+2]>$content[$i+3]){  #no gap and with overlap
              my $bac_extract_seq="";
              my $bac_extract_len=0;
              if($content[$i+1]==0){  #bac no reverse
                      if($content[$i+2]>=$bac_len{$content[$i]}){
                               $content[$i+2]=$bac_len{$content[$i]}-1;
                      }
                      my $extract_len=$content[$i+2]-$content[$i+3]+1;
                      $extract_seq=substr($bac_seq{$content[$i]},$content[$i+3],$extract_len);
                      $bac_extract_seq=$extract_seq;
                      $bac_extract_len=$extract_len;
              }
              elsif($content[$i+1]==1){  #bac reverse
                      if($content[$i+2]>=$bac_len{$content[$i]}){
                               $content[$i+2]=$bac_len{$content[$i]}-1;
                      }
                      my $extract_len=$content[$i+2]-$content[$i+3]+1;
#                      $bac_seq{$content[$i]}=reverse_complement($bac_seq{$content[$i]});
                      my $Current_Bac=reverse_complement($bac_seq{$content[$i]});
                      $extract_seq=substr($Current_Bac,$content[$i+3],$extract_len);
                      $bac_extract_seq=$extract_seq;
                      $bac_extract_len=$extract_len;
              }
              if($content[$i-3]==0){    #left ctg no reverse
                      if($content[$i-1]>=$ctg_len{$content[$i-4]}){
                               $content[$i-1]=$ctg_len{$content[$i-4]}-1;
                      }
                      my $extract_len=$content[$i-1]-$content[$i-2]-$bac_extract_len+1;
                      $extract_seq=substr($ctg_seq{$content[$i-4]},$content[$i-2],$extract_len);
                      $seq=$seq.$extract_seq;
               }
               elsif($content[$i-3]==1){ #left ctg reverse
                      if($content[$i-1]>=$ctg_len{$content[$i-4]}){
                               $content[$i-1]=$ctg_len{$content[$i-4]}-1;
                      }
#                      $ctg_seq{$content[$i-4]}=reverse_complement($ctg_seq{$content[$i-4]});
                      my $Current_Bac=reverse_complement($ctg_seq{$content[$i-4]});
                      my $extract_len=$content[$i-1]-$content[$i-2]-$bac_extract_len+1;
                      $extract_seq=substr($Current_Bac,$content[$i-2],$extract_len);
                      $seq=$seq.$extract_seq;
              }
              $seq=$seq.$bac_extract_seq;
              if($content[$i+5]==0){    #right ctg no reverse
                      if($content[$i+7]>=$ctg_len{$content[$i+4]}){
                               $content[$i+7]=$ctg_len{$content[$i+4]}-1;
                      }
                      my $extract_len=$content[$i+7]-$content[$i+6]-$bac_extract_len+1;
                      $content[$i+6]=$content[$i+6]+$bac_extract_len;
                      $last_ctg_start=$content[$i+6];
#                      $extract_seq=substr($ctg_seq{$content[$i+4]},$ctg_start,$extract_len);
#                      $seq=$seq.$extract_seq;
               }
               elsif($content[$i+5]==1){ #right ctg reverse
                      if($content[$i+7]>=$ctg_len{$content[$i+4]}){
                               $content[$i+7]=$ctg_len{$content[$i+4]}-1;
                      }
                      $content[$i+6]=$content[$i+6]+$bac_extract_len;
                      $last_ctg_start=$content[$i+6];
=pod
                      $ctg_seq{$content[$i+4]}=reverse_complement($ctg_seq{$content[$i+4]});
                      my $extract_len=$content[$i+7]-$content[$i+6]-$bac_extract_len+1;
                      my $ctg_start=$content[$i+6]+$bac_extract_len;
                      $extract_seq=substr($ctg_seq{$content[$i+4]},$ctg_start,$extract_len);
                      $seq=$seq.$extract_seq;
=cut
              }
          }
      }
      if($content[-3]==0){
          if($content[-1]>=$ctg_len{$content[-4]}){
              $content[-1]=$ctg_len{$content[-4]}-1;
          }        
          my $last_seq_len=$content[-1]-$last_ctg_start+1;
          my $last_seq=substr($ctg_seq{$content[-4]},$last_ctg_start,$last_seq_len);
          $seq=$seq.$last_seq;
      }
      elsif($content[-3]==1){
          if($content[-1]>=$ctg_len{$content[-4]}){
              $content[-1]=$ctg_len{$content[-4]}-1;
          }    
          my $last_seq_len=$content[-1]-$last_ctg_start+1;
#          $ctg_seq{$content[-4]}=reverse_complement($ctg_seq{$content[-4]});
          my $Current_Bac=reverse_complement($ctg_seq{$content[-4]});
          my $last_seq=substr($Current_Bac,$last_ctg_start,$last_seq_len);
          $seq=$seq.$last_seq;
      }   
#      print OUT ">$count_ctg\n";                          #Do not keep name information
      print OUT ">$Ctg_Name\n";                            #Keeping the original name information
      print OUT "$seq\n";     
      $count_ctg++;
      foreach my $key1 (keys %ctg_per_cluster){
            print  "$ctg_len{$key1}\n";
      }
      print "\n";
}
open OUT1,">non-connected-ctg.fasta" or die $!;
my $count=0;
my $all_len=0;
foreach my $key (keys %ctg_seq){
      next if(exists $ctg_existed{$key});
      $count++;
      $all_len=$all_len+length($ctg_seq{$key});
      print OUT1 ">$key\n";
      print OUT1 "$ctg_seq{$key}\n";
}
print "$all_len\n";

            
close(OUT);     
close(IN1);
close(IN2);
close(IN3);
close(OUT1);                

 
                
        

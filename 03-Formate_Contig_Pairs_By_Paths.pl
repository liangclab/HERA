#!/usr/bin/perl
use warnings;
use strict;

#Contig_Pairs_With_Overlaps.txt
my $infile1=shift;

#ctg_clusters.txt
my $infile2=shift;

#matrix with contig pairs
my $outfile=shift;

open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;
open OUT,">$outfile" or die $!;

my %Contig_Pair_With_Overlap=();
#3       Head    758     Tail    0       18452   1296480 821070  839333  839333  98.39   +
while(<IN1>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$Contig_Pair_With_Overlap{$content[0]}{$content[1]}{$content[2]}{$content[3]}=$line;
}
close(IN1);

my %Contig_Pair_With_Path=();
#R498_1627743_1/1808_11897 R498_4016954_1/0_16938 R498_3245081_1/183_16057 R498_5887543_1/0_14019 R498_4641610_1/283_14760 R498_5779396_1/498_14823 contig0107_size120001 R498_760764_1/14_22641 R498_5998020_1/0_25161 R498_5794087_1/44_22587 R498_5881788_1/32_21556 R498_5971505_1/12_12575  SuperScaffold_1.10-Tail-SuperScaffold_1.11-Head-BestMatching-Passing
my %Used_Lines=();
while(<IN2>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	my $positive="";
	my $reverse="";
	next if($line=~/^Real/);
	for(my $i=0;$i<@content-1;$i++){
		$positive=$positive.$content[$i];
		$reverse=$content[$i].$reverse;
	}
	$positive=$positive.$content[-1];
	$reverse=$reverse.$content[-1];
	next if(exists $Used_Lines{$positive} || exists $Used_Lines{$reverse});
	$Used_Lines{$positive}=0;
	$Used_Lines{$reverse}=0;
	my @info=split /-/,$content[-1];
	if(exists $Contig_Pair_With_Path{$info[0]}{$info[1]}{$info[2]}{$info[3]}{$info[4]}){
		$Contig_Pair_With_Path{$info[0]}{$info[1]}{$info[2]}{$info[3]}{$info[4]}++;
	}
	elsif(exists $Contig_Pair_With_Path{$info[2]}{$info[3]}{$info[0]}{$info[1]}{$info[4]}){
		$Contig_Pair_With_Path{$info[2]}{$info[3]}{$info[0]}{$info[1]}{$info[4]}++;
	}
	else{
		$Contig_Pair_With_Path{$info[0]}{$info[1]}{$info[2]}{$info[3]}{$info[4]}=1;
	}
}
close(IN2);

my %Recond_Used=();
foreach my $left_ctg (keys %Contig_Pair_With_Path){
	my $temp1=$Contig_Pair_With_Path{$left_ctg};
	foreach my $left_ori (keys %$temp1){
		my $temp2=$Contig_Pair_With_Path{$left_ctg}{$left_ori};
		foreach my $right_ctg (keys %$temp2){
			my $temp3=$Contig_Pair_With_Path{$left_ctg}{$left_ori}{$right_ctg};
			foreach my $right_ori (keys %$temp3){
				my $temp4=$Contig_Pair_With_Path{$left_ctg}{$left_ori}{$right_ctg}{$right_ori};
				my $total=0;
				my $BestMatch=0;
				my $MostExtend=0;
				my $RandomExtend=0;
				my $Overlap=0;
				if(exists $Contig_Pair_With_Path{$left_ctg}{$left_ori}{$right_ctg}{$right_ori}{"BestMatching"}){
					$BestMatch=$Contig_Pair_With_Path{$left_ctg}{$left_ori}{$right_ctg}{$right_ori}{"BestMatching"};
				}
				elsif(exists $Contig_Pair_With_Path{$right_ctg}{$right_ori}{$left_ctg}{$left_ori}{"BestMatching"}){
					$BestMatch=$Contig_Pair_With_Path{$right_ctg}{$right_ori}{$left_ctg}{$left_ori}{"BestMatching"};
				}
				if(exists $Contig_Pair_With_Path{$left_ctg}{$left_ori}{$right_ctg}{$right_ori}{"MostExtending"}){
					$MostExtend=$Contig_Pair_With_Path{$left_ctg}{$left_ori}{$right_ctg}{$right_ori}{"MostExtending"};
				}
				elsif(exists $Contig_Pair_With_Path{$right_ctg}{$right_ori}{$left_ctg}{$left_ori}{"MostExtending"}){
					$MostExtend=$Contig_Pair_With_Path{$right_ctg}{$right_ori}{$left_ctg}{$left_ori}{"MostExtending"};
				}
				if(exists $Contig_Pair_With_Path{$left_ctg}{$left_ori}{$right_ctg}{$right_ori}{"RandomExtending"}){
					$RandomExtend=$Contig_Pair_With_Path{$left_ctg}{$left_ori}{$right_ctg}{$right_ori}{"RandomExtending"};
				}
				elsif(exists $Contig_Pair_With_Path{$right_ctg}{$right_ori}{$left_ctg}{$left_ori}{"RandomExtending"}){
					$RandomExtend=$Contig_Pair_With_Path{$right_ctg}{$right_ori}{$left_ctg}{$left_ori}{"RandomExtending"};
				}
				if(exists $Contig_Pair_With_Overlap{$left_ctg}{$left_ori}{$right_ctg}{$right_ori} || exists $Contig_Pair_With_Overlap{$right_ctg}{$right_ori}{$left_ctg}{$left_ori}){
					$Overlap=1;
					$Recond_Used{$left_ctg}{$left_ori}{$right_ctg}{$right_ori}=0;
					$Recond_Used{$right_ctg}{$right_ori}{$left_ctg}{$left_ori}=0;
				}
				$total=$BestMatch+$MostExtend+$RandomExtend;
				print OUT "$left_ctg\t$left_ori\t$right_ctg\t$right_ori\t$BestMatch\t$MostExtend\t$RandomExtend\t$total\t$Overlap\n";
			}
		}
	}
}
close(IN2);
foreach my $key1 (keys %Contig_Pair_With_Overlap){
	my $temp1=$Contig_Pair_With_Overlap{$key1};
	foreach my $key2 (keys %$temp1){
		my $temp2=$Contig_Pair_With_Overlap{$key1}{$key2};
		foreach my $key3 (keys %$temp2){
			my $temp3=$Contig_Pair_With_Overlap{$key1}{$key2}{$key3};
			foreach my $key4 (keys %$temp3){
				if(!exists $Recond_Used{$key1}{$key2}{$key3}{$key4}){
					print OUT "$key1\t$key2\t$key3\t$key4\t0\t0\t0\t0\t1\n";
				}
			}
		}
	}
}
close(OUT);


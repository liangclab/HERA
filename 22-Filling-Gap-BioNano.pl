#!/usr/bin/perl
#################################################################################################################
#            Filling the gap by pathcontig with the alignment of daligner
#                Author: huilong du
#           
#
#################################################################################################################
use warnings;
use strict;
my $infile1=shift;                               #Scaffold2Ctg_Gap.txt
my $infile3=shift;                               #Prosudo_ScaffoldNonEnzyme2Contig.fasta
my $infile4=shift;                               #PathContigRename.fasta
my $outfile=shift;                               #Connected sequence
my $minidentity=0.96;
my $maxoverhang=10000;
my $minoverlap=2000;

open IN1,"<$infile1" or die $!;
open IN3,"<$infile3" or die $!;
open IN4,"<$infile4" or die $!;
open OUT,">$outfile" or die $!;

open PATH,">Selected_Path.txt" or die $!;
open POS,">Ctg_Position.txt" or die $!;
open INFO,">SuperContig_Part_Info.txt" or die $!;
#Super-Scaffold_3.1      Super-Scaffold_3.2      499     364273  364771
my %Contig_Position=();
my %Contig_Pair=();
my %Total_Contig=();
my %Contig_In_Scaffold=();
while(<IN1>){
	chomp;
	my $info=$_;
	my @infos=split /\s+/,$info;
	$infos[0]=~/(Super-Scaffold_\d+)\.(\d+)/;
	my $scaffold=$1;
	my $First=$2;
	my $Second=$First+1;
	my $Gap_Len=$infos[2];
	$Contig_In_Scaffold{$infos[0]}=0;
	$Contig_In_Scaffold{$infos[1]}=0;
	print PATH "#$infos[0]\t$infos[1]\n";
	my $Gap_Info=$scaffold."-".$First."-".$Second;
	open IN,"<$Gap_Info/$Gap_Info-Final_Reformated.txt" or die $!;
	my %PairContig_Gap=();               
        my %PathContig_Info=();              
#+       Super-Scaffold_1067.1   4077    5089    364118  Super-Scaffold_1067.2   12910   13930   314580  153
#+       Super-Scaffold_1067.1   4077    5089    364118  Super-Scaffold_1067-1-2.58017   44552   45571   55596   154
#+       Super-Scaffold_1067.2   3314    23988   314580  Super-Scaffold_1067-1-2.58017   35051   55596   55596   331
	my $Loverlap=0;
	my $Roverlap=0;
	my $Loverhang=0;
	my $Roverhang=0;
	my $Max_Overlap_Score=0;
	my $Max_Path_Score=0;
	my $Max_Overlap_Line="";
	my $Max_Path_Line_Left="";
	my $Max_Path_Line_Right="";
	while(<IN>){
		chomp;
		my $line=$_;
		my @content=split /\s+/,$line;
		next if(@content!=10);
		if($content[1]=~/Super-Scaffold_\d+\.\d+/ && $content[5]=~/Super-Scaffold_\d+\.\d+/){
			$Total_Contig{$content[1]}=$content[4];
			$Total_Contig{$content[5]}=$content[8];
		}
		next if($content[0] ne "+");
		next if($content[3]-$content[2]<$minoverlap);
		my $identity=1-($content[9]/($content[3]-$content[2]+$content[7]-$content[6]))*2;
		next if($identity<$minidentity);
		if($Gap_Len<500){
			if($content[1]=~/Super-Scaffold_\d+\.(\d+)/ && $content[5]=~/Super-Scaffold_\d+\.(\d+)/){
				$content[1]=~/Super-Scaffold_\d+\.(\d+)/;
				my $First_Ctg=$1;
				$content[5]=~/Super-Scaffold_\d+\.(\d+)/;
				my $Second_Ctg=$1;
				if($First_Ctg==$Second_Ctg-1){
					$Loverlap=$content[3]-$content[2];
					$Loverhang=$content[4]-$content[3];
					$Roverlap=$content[7]-$content[6];
					$Roverhang=$content[6];
					if($Loverlap>=$Loverhang && $Roverlap>=$Roverhang){
						my $Overlap_Score=($Loverlap+$Roverlap)*$identity-$Loverhang-$Roverhang;
						if($Overlap_Score>$Max_Overlap_Score){
							$Max_Overlap_Score=$Overlap_Score;
							$Max_Overlap_Line=$line;
						}
					}
				}
			}
			elsif($content[1]=~/Super-Scaffold_\d+\.(\d+)/ && $content[5]=~/Super-Scaffold_\d+-(\d+)-(\d+)\.(\d+)/){
				$Total_Contig{$content[1]}=$content[4];
				$content[1]=~/Super-Scaffold_\d+\.(\d+)/;
				my $First_Ctg=$1;
				$content[5]=~/Super-Scaffold_\d+-(\d+)-(\d+)\.(\d+)/;
				my $Path_First=$1;
				my $Path_Second=$2;
				if($First_Ctg==$Path_First){
					$Loverlap=$content[3]-$content[2];
                                        $Loverhang=$content[4]-$content[3];
                                        $Roverlap=$content[7]-$content[6];
                                        $Roverhang=$content[6];
					next if($Loverlap<$Loverhang || $Roverlap<$Roverhang);
					if(!exists $PathContig_Info{$content[5]}{'left'}){
						$PathContig_Info{$content[5]}{'left'}{0}=$line;
						$PathContig_Info{$content[5]}{'left'}{1}=($Loverlap+$Roverlap)*$identity-$Loverhang-$Roverhang;
					}
					else{
						my $Overlap_Score=($Loverlap+$Roverlap)*$identity-$Loverhang-$Roverhang;
						if($Overlap_Score>$PathContig_Info{$content[5]}{'left'}{1}){
							$PathContig_Info{$content[5]}{'left'}{0}=$line;
							$PathContig_Info{$content[5]}{'left'}{1}=$Overlap_Score;
						}
					}
				}
				elsif($First_Ctg==$Path_Second){
					$Loverlap=$content[7]-$content[6];
					$Loverhang=$content[8]-$content[7];
					$Roverlap=$content[3]-$content[2];
					$Roverhang=$content[2];
					next if($Loverlap<$Loverhang || $Roverlap<$Roverhang);
					if(!exists $PathContig_Info{$content[5]}{'right'}){
                                                $PathContig_Info{$content[5]}{'right'}{0}=$line;
                                                $PathContig_Info{$content[5]}{'right'}{1}=($Loverlap+$Roverlap)*$identity-$Loverhang-$Roverhang;
                                        }
                                        else{
                                                my $Overlap_Score=($Loverlap+$Roverlap)*$identity-$Loverhang-$Roverhang;
                                                if($Overlap_Score>$PathContig_Info{$content[5]}{'right'}{1}){
                                                        $PathContig_Info{$content[5]}{'right'}{0}=$line;
                                                        $PathContig_Info{$content[5]}{'right'}{1}=$Overlap_Score;
                                                }
                                        }
				}
			}
		}
		elsif($Gap_Len>=500){
			if($content[1]=~/Super-Scaffold_\d+\.(\d+)/ && $content[5]=~/Super-Scaffold_\d+-(\d+)-(\d+)\.(\d+)/){
				$Total_Contig{$content[1]}=$content[4];
				$content[1]=~/Super-Scaffold_\d+\.(\d+)/;
                                my $First_Ctg=$1;
				$content[5]=~/Super-Scaffold_\d+-(\d+)-(\d+)\.(\d+)/;
                                my $Path_First=$1;
                                my $Path_Second=$2;
                                if($First_Ctg==$Path_First){
                                        $Loverlap=$content[3]-$content[2];
                                        $Loverhang=$content[4]-$content[3];
                                        $Roverlap=$content[7]-$content[6];
                                        $Roverhang=$content[6];
                                        next if($Loverlap<$Loverhang || $Roverlap<$Roverhang);
                                        if(!exists $PathContig_Info{$content[5]}{'left'}){
                                                $PathContig_Info{$content[5]}{'left'}{0}=$line;
                                                $PathContig_Info{$content[5]}{'left'}{1}=($Loverlap+$Roverlap)*$identity-$Loverhang-$Roverhang;
                                        }
                                        else{
                                                my $Overlap_Score=($Loverlap+$Roverlap)*$identity-$Loverhang-$Roverhang;
                                                if($Overlap_Score>$PathContig_Info{$content[5]}{'left'}{1}){
                                                        $PathContig_Info{$content[5]}{'left'}{0}=$line;
                                                        $PathContig_Info{$content[5]}{'left'}{1}=$Overlap_Score;
                                                }
                                        }
                                }
				elsif($First_Ctg==$Path_Second){
                                        $Loverlap=$content[7]-$content[6];
                                        $Loverhang=$content[8]-$content[7];
                                        $Roverlap=$content[3]-$content[2];
                                        $Roverhang=$content[2];
                                        next if($Loverlap<$Loverhang || $Roverlap<$Roverhang);
                                        if(!exists $PathContig_Info{$content[5]}{'right'}){
                                                $PathContig_Info{$content[5]}{'right'}{0}=$line;
                                                $PathContig_Info{$content[5]}{'right'}{1}=($Loverlap+$Roverlap)*$identity-$Loverhang-$Roverhang;
                                        }
                                        else{
                                                my $Overlap_Score=($Loverlap+$Roverlap)*$identity-$Loverhang-$Roverhang;
                                                if($Overlap_Score>$PathContig_Info{$content[5]}{'right'}{1}){
                                                        $PathContig_Info{$content[5]}{'right'}{0}=$line;
                                                        $PathContig_Info{$content[5]}{'right'}{1}=$Overlap_Score;
                                                }
                                        }
                                }
                        }
                }
	}
	close(IN);
#+       Super-Scaffold_1067.1   4077    5089    364118  Super-Scaffold_1067.2   12910   13930   314580  153
#+       Super-Scaffold_1067.1   4077    5089    364118  Super-Scaffold_1067-1-2.58017   44552   45571   55596   154
#+       Super-Scaffold_1067.2   3314    23988   314580  Super-Scaffold_1067-1-2.58017   35051   55596   55596   331
        #selecting the best match path based on the gap length
        my $Best_Path_Score=0;
	my $Best_Path_Left="";               #recording the left alignment of best path
	my $Best_Path_Right="";		     #recording the right alignment of best path
	my $Best_Path_Gap=9999999999;

	foreach my $path (keys %PathContig_Info){
		if($Gap_Len<500){ #if with overlap, selecting the highest score path
			next if(!exists $PathContig_Info{$path}{'left'} or !exists $PathContig_Info{$path}{'right'});
			if($PathContig_Info{$path}{'left'}{1}+$PathContig_Info{$path}{'right'}{1}>$Best_Path_Score){
				$Best_Path_Score=$PathContig_Info{$path}{'left'}{1}+$PathContig_Info{$path}{'right'}{1};
				$Best_Path_Left=$PathContig_Info{$path}{'left'}{0};
				$Best_Path_Right=$PathContig_Info{$path}{'right'}{0};
			}
		}
		elsif($Gap_Len>=500){ #if without overlap, selecting the best match path
			next if(!exists $PathContig_Info{$path}{'left'} or !exists $PathContig_Info{$path}{'right'});
			my @Left_Line=split /\s+/,$PathContig_Info{$path}{'left'}{0};
			my @Right_Line=split /\s+/,$PathContig_Info{$path}{'right'}{0};
			my $Path_Len_Used=$Right_Line[7]-$Left_Line[6];
			my $Gap_Len_Used=$Left_Line[4]-$Left_Line[2]+$Right_Line[3]+$Gap_Len;
			my $Path_Gap=abs ($Path_Len_Used-$Gap_Len_Used);
			if($Path_Gap<$Best_Path_Gap){
				$Best_Path_Gap=$Path_Gap;
				$Best_Path_Left=$PathContig_Info{$path}{'left'}{0};
				$Best_Path_Right=$PathContig_Info{$path}{'right'}{0};
			}
		}
	}
#	print "$Best_Path_Left\n$Best_Path_Right\n";
	#recording the part of Scaffold_Contig and Path_Contig used to construct the final super-contig
	if($Gap_Len<500){        #Contig may be with overlap
		if($Max_Overlap_Line ne ""){
			my @content=split /\s+/,$Max_Overlap_Line;
			$Contig_Position{$content[1]}{'End'}=$content[3];         #Left_Contig End    
			$Contig_Position{$content[1]}{'Len'}=$content[4];         #Left_Contig Length
			$Contig_Position{$content[5]}{'Start'}=$content[7];	  #Right_Contig Start
			$Contig_Position{$content[5]}{'Len'}=$content[8];	  #Right_Contig Length
			$Contig_Pair{$content[1]}{$content[5]}=0;
			print PATH "$Max_Overlap_Line\n";
		}
		elsif($Best_Path_Left ne "" && $Best_Path_Right ne ""){
			my @Left_Line=split /\s+/,$Best_Path_Left;
			my @Right_Line=split /\s+/,$Best_Path_Right;
			print PATH "$Best_Path_Left\n$Best_Path_Right\n";
			$Contig_Position{$Left_Line[1]}{'End'}=$Left_Line[3]; #left contig end
			$Contig_Position{$Left_Line[1]}{'Len'}=$Left_Line[4]; #left contig length
			$Contig_Position{$Left_Line[1]}{'Path'}=$Left_Line[5]; #Path
			$Contig_Position{$Right_Line[1]}{'Start'}=$Right_Line[3];#right contig start
			$Contig_Position{$Right_Line[1]}{'Len'}=$Right_Line[4];#right contig length
			$Contig_Pair{$Left_Line[1]}{$Right_Line[1]}=0;

			$Contig_Position{$Left_Line[5]}{'Start'}=$Left_Line[7];#path contig start
			$Contig_Position{$Left_Line[5]}{'End'}=$Right_Line[7];#path contig end
			$Contig_Position{$Left_Line[5]}{'Len'}=$Left_Line[8]; #path contig length
		}	
	}
	else{                    #Contig without overlap
		next if($Best_Path_Left eq "" || $Best_Path_Right eq "");
		my @Left_Line=split /\s+/,$Best_Path_Left;
                my @Right_Line=split /\s+/,$Best_Path_Right;
		if($Best_Path_Left eq ""){
			print STDERR "$Best_Path_Left\n";
		}
                print PATH "$Best_Path_Left\n$Best_Path_Right\n";
                $Contig_Position{$Left_Line[1]}{'End'}=$Left_Line[3];
                $Contig_Position{$Left_Line[1]}{'Len'}=$Left_Line[4];
		$Contig_Position{$Left_Line[1]}{'Path'}=$Left_Line[5];
                $Contig_Position{$Right_Line[1]}{'Start'}=$Right_Line[3];
                $Contig_Position{$Right_Line[1]}{'Len'}=$Right_Line[4];
		$Contig_Pair{$Left_Line[1]}{$Right_Line[1]}=0;

		$Contig_Position{$Left_Line[5]}{'Start'}=$Left_Line[7];
		$Contig_Position{$Left_Line[5]}{'End'}=$Right_Line[7];
		$Contig_Position{$Left_Line[5]}{'Len'}=$Left_Line[8];
	}		
}
close(IN1);
my $sign="";
my %Path_Contig_Seq=();
while(<IN4>){
        chomp;
        my $line=$_;
        if($line=~/^>(\S+)/){
                $sign=$1;
        }
        else{
                $Path_Contig_Seq{$sign}=$line;
		$Total_Contig{$sign}=length($line);
        }
}
close(IN4);

foreach my $key (keys %Total_Contig){
	if(!exists $Contig_Position{$key}){
		$Contig_Position{$key}{'Start'}=0;
		$Contig_Position{$key}{'End'}=$Total_Contig{$key};
		$Contig_Position{$key}{'Len'}=$Total_Contig{$key};
	}
	elsif(!exists $Contig_Position{$key}{'Start'}){
		$Contig_Position{$key}{'Start'}=0;
	}
	elsif(!exists $Contig_Position{$key}{'End'}){
		$Contig_Position{$key}{'End'}=$Total_Contig{$key};
	}
}

my %Scaffold_Contig_Seq=();
$sign="";
while(<IN3>){
	chomp;
	my $line=$_;
	if($line=~/^>(\S+)/){
		$sign=$1;
	}
	else{
		$Scaffold_Contig_Seq{$sign}=$line;
		if(!exists $Contig_In_Scaffold{$sign}){
			print OUT ">$sign\n";
			print OUT "$line\n";
			my $Whole_Len=length($line);
			print POS "$sign\t0\t$Whole_Len\n";
		}	
	}
}
close(IN3);


open IN1,"<$infile1" or die $!;
#Super-Scaffold_3.1      Super-Scaffold_3.2      499     364273  364771
my $Current_Seq="";
my $name="";
my $count=1;
my $Current_Start=0;
my $Current_Scaffold="";
my %Current_Total_Len=();
my $Current_Count=1;
while(<IN1>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$content[0]=~/(Super-Scaffold_\d+)\.\d+/;
	my $Scaffold=$1;
	print POS "$content[0]\t$Contig_Position{$content[0]}{'Start'}\t$Contig_Position{$content[0]}{'End'}\n";
	print POS "$content[1]\t$Contig_Position{$content[1]}{'Start'}\t$Contig_Position{$content[1]}{'End'}\n";
	if($Current_Scaffold ne $Scaffold){
		$Current_Count=1;
		if($name eq ""){
			$Current_Scaffold=$Scaffold;
		}
		else{
			print OUT ">$name\n";
			print OUT "$Current_Seq\n";
			foreach my $part (sort {$a<=>$b} keys %Current_Total_Len){
				print INFO "$name\t$Current_Total_Len{$part}{'start'}\t$Current_Total_Len{$part}{'end'}\t$Current_Total_Len{$part}{'type'}\n";
			}
			%Current_Total_Len=();
		}
		my $First_Len=$Contig_Position{$content[0]}{'End'}-$Contig_Position{$content[0]}{'Start'};
                my $First_Seq=substr($Scaffold_Contig_Seq{$content[0]},$Contig_Position{$content[0]}{'Start'},$First_Len);
		$Current_Seq=$First_Seq;
		$name=$content[0];
		$Current_Scaffold=$Scaffold;
	}
	if(exists $Contig_Pair{$content[0]}{$content[1]}){
		$name=$name."|".$content[1];
		my $Path_Seq="";
		if(exists $Contig_Position{$content[0]}{'Path'}){
       	                my $Path_Len=$Contig_Position{$Contig_Position{$content[0]}{'Path'}}{'End'}-$Contig_Position{$Contig_Position{$content[0]}{'Path'}}{'Start'};

			$Current_Total_Len{$Current_Count}{'start'}=length($Current_Seq)+1;
			$Current_Total_Len{$Current_Count}{'end'}=length($Current_Seq)+$Path_Len+1;
			$Current_Total_Len{$Current_Count}{'type'}="Path";
			$Current_Count++;

               	        $Path_Seq=substr($Path_Contig_Seq{$Contig_Position{$content[0]}{'Path'}},$Contig_Position{$Contig_Position{$content[0]}{'Path'}}{'Start'},$Path_Len);
			print POS "$Contig_Position{$content[0]}{'Path'}\t$Contig_Position{$Contig_Position{$content[0]}{'Path'}}{'Start'}\t$Contig_Position{$Contig_Position{$content[0]}{'Path'}}{'End'}\n";
              	}
		else{
			$Current_Total_Len{$Current_Count}{'start'}=length($Current_Seq)-$Contig_Position{$content[1]}{'Start'};
			$Current_Total_Len{$Current_Count}{'end'}=length($Current_Seq)+$Contig_Position{$content[1]}{'Start'};
			$Current_Total_Len{$Current_Count}{'type'}="Overlap";
			$Current_Count++;
		}
		my $Second_Len=$Contig_Position{$content[1]}{'End'}-$Contig_Position{$content[1]}{'Start'};
		my $Second_Seq=substr($Scaffold_Contig_Seq{$content[1]},$Contig_Position{$content[1]}{'Start'},$Second_Len);
		
		$Current_Seq=$Current_Seq.$Path_Seq.$Second_Seq;
		$Current_Scaffold=$Scaffold;
	}
	else{
		print OUT ">$name\n";
		print OUT "$Current_Seq\n";
		foreach my $part (sort {$a<=>$b} keys %Current_Total_Len){
                        print INFO "$name\t$Current_Total_Len{$part}{'start'}\t$Current_Total_Len{$part}{'end'}\t$Current_Total_Len{$part}{'type'}\n";
                }
                %Current_Total_Len=();
		$name=$content[1];
		my $Second_Len=$Contig_Position{$content[1]}{'End'}-$Contig_Position{$content[1]}{'Start'};
                my $Second_Seq=substr($Scaffold_Contig_Seq{$content[1]},$Contig_Position{$content[1]}{'Start'},$Second_Len);
		$Current_Seq=$Second_Seq;
                $Current_Scaffold=$Scaffold;
	}
}
foreach my $part (sort {$a<=>$b} keys %Current_Total_Len){
	print INFO "$name\t$Current_Total_Len{$part}{'start'}\t$Current_Total_Len{$part}{'end'}\t$Current_Total_Len{$part}{'type'}\n";
}
close(INFO);
print OUT ">$name\n";
print OUT "$Current_Seq\n";
close(IN1);
close(OUT);

#!/usr/bin/perl
use warnings;
use strict;

sub reverse_complement{
        my $input_seq=shift;
        my $reverse=reverse $input_seq;
        $reverse=~tr/ACGT/TGCA/;
        $reverse=~tr/acgt/tgca/;
        return $reverse;
}
#Contig and Orientation with scores
my $infile=shift;
my $infile1=shift;
my $outfile=shift;
open IN,"<$infile" or die $!;
open IN1,"<$infile1" or die $!;
open OUT,">$outfile" or die $!;


my %Contig_Seq=();
my $sign="";
while(<IN1>){
        chomp;
        my $line=$_;
        if($line=~/^>(\S+)/){
                $sign=$1;
        }
        else{
                $Contig_Seq{$sign}=$line;
        }
}
close(IN1);
#ctg7180000004696        Head    ctg7180000006755        Head    0       18703   198339  782256  800990  800990  99.63   -
my %Contig_Pairs_Score=();
while(<IN>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	next if($content[0] eq $content[2]);
=pod
	my $Fleft_Overhang=$content[4];
	my $Fright_Overhang=$content[6]-$content[5];
	my $Sleft_Overhang=$content[7];
	my $Sright_Overhang=$content[9]-$content[8];
	my $Left_Overhang=0;
	my $Right_Overhang=0;
	if($Fleft_Overhang>=$Fright_Overhang){
		$Left_Overhang=$Fright_Overhang;
	}
	else{
		$Left_Overhang=$Fleft_Overhang;
	}
	if($Sleft_Overhang>=$Sright_Overhang){
		$Right_Overhang=$Sright_Overhang;
	}
	else{
		$Right_Overhang=$Sleft_Overhang;
	}
	my $score=($content[5]-$content[4]+$content[8]-$content[7]-$Left_Overhang-$Right_Overhang)*$content[10]/2;
=cut
	$Contig_Pairs_Score{$line}=$content[7];
}
close(IN);


#ctg7180000007712        Tail    ctg7180000007713        Head    37      129     0       166     0
my %Existed_Contig_Orientation=();
my %Contig_Pair=();
my %Merged_Path=();
foreach my $alignment (sort {$Contig_Pairs_Score{$b} <=> $Contig_Pairs_Score{$a}} keys %Contig_Pairs_Score){
	my @content=split /\s+/,$alignment;
	next if(exists $Existed_Contig_Orientation{$content[0]}{$content[1]} || exists $Existed_Contig_Orientation{$content[2]}{$content[3]} || exists $Contig_Pair{$content[0]}{$content[2]} || exists $Contig_Pair{$content[2]}{$content[0]});
	my $reverse_alignment=$content[2]."\t".$content[3]."\t".$content[0]."\t".$content[1]."\t".$content[4]."\t".$content[5]."\t".$content[6]."\t".$content[7]."\t".$content[8];
	$Existed_Contig_Orientation{$content[0]}{$content[1]}=$alignment;
	$Existed_Contig_Orientation{$content[2]}{$content[3]}=$reverse_alignment;
	$Contig_Pair{$content[0]}{$content[2]}=$alignment;
	$Contig_Pair{$content[2]}{$content[0]}=$reverse_alignment;
	my $path="";
	my $temp="";
	my $First_Path="";
	my $Second_Path="";
	foreach my $key (keys %Merged_Path){
		my @info=split /--/,$key;
		my $start=$info[0];
		my $end=$info[-1];
		if(($start eq $content[0])||($start eq $content[2])||($end eq $content[0])||($end eq $content[2])){
			if($First_Path eq ""){
				$First_Path=$key;
			}
			else{
				$Second_Path=$key;
			}
		}
	}
	if($First_Path ne "" && $Second_Path ne ""){
		my @cluster1=split /--/,$First_Path;
		my @cluster2=split /--/,$Second_Path;
		my $First_Start=$cluster1[0];
		my $First_End=$cluster1[-1];
		my $Second_Start=$cluster2[0];
		my $Second_End=$cluster2[-1];
		if(($First_End eq $content[0] || $First_End eq $content[2]) && ($content[2] eq $Second_Start ||$content[0] eq $Second_Start)){
			$path=$First_Path."--".$Second_Path;
		}
		elsif(($Second_End eq $content[0] || $Second_End eq $content[2])&&($content[2] eq $First_Start || $content[0] eq $First_Start)){
			$path=$Second_Path."--".$First_Path;
		}
		elsif((($First_Start eq $content[0] || $First_Start eq $content[2])&&($Second_Start eq $content[0] || $Second_Start eq $content[2]))){
			$path=$cluster1[0];
			for(my $m=1;$m<@cluster1;$m++){
				$path=$cluster1[$m]."--".$path;
			}
			$path=$path."--".$Second_Path;
		}
		elsif(($First_End eq $content[0] || $First_End eq $content[2])&&($Second_End eq $content[0] || $Second_End eq $content[2])){
			$path=$First_Path;
			for(my $m=@cluster2-1;$m>=0;$m--){
				$path=$path."--".$cluster2[$m];
			}
		}
		$Merged_Path{$path}=0;
		delete($Merged_Path{$First_Path});
		delete($Merged_Path{$Second_Path});
	}
	elsif($First_Path ne ""){
		my @info=split /--/,$First_Path;
                my $start=$info[0];
                my $end=$info[-1];
		if($start eq $content[0]){
			next if($First_Path=~/^$content[2]--/ ||$First_Path=~/--$content[2]$/);
			$path=$content[2]."--".$First_Path;
		}
		elsif($start eq $content[2]){
			next if($First_Path=~/^$content[0]--/ ||$First_Path=~/--$content[0]$/);
			$path=$content[0]."--".$First_Path;
		}
		if($end eq $content[0]){
			next if($First_Path=~/^$content[2]--/ ||$First_Path=~/--$content[2]$/);
			$path=$First_Path."--".$content[2];
		}
		elsif($end eq $content[2]){
			next if($First_Path=~/^$content[0]--/ ||$First_Path=~/--$content[0]$/);
			$path=$First_Path."--".$content[0];
		}
		$Merged_Path{$path}=0;
		delete($Merged_Path{$First_Path});
	}
	else{
		$path=$content[0]."--".$content[2];
		$Merged_Path{$path}=0;
	}
	print "$alignment\t$Contig_Pairs_Score{$alignment}\n";
}


=pod
=cut

my $unknown="N"x200;
#ctg7180000007712        Tail    ctg7180000007713        Head    37      129     0       166     0
my %Merged_Contig=();
foreach my $contig (keys %Merged_Path){
	my @info=split /--/,$contig;
	print OUT ">$contig\n";
	my $current_chain=0;
	my $current_name=$info[0];
	my $current_seq="";
	my $line=$Contig_Pair{$info[0]}{$info[1]};
	my @align=split /\s+/,$line;
	$Merged_Contig{$info[0]}=0;
	my $reverse=0;
	if($align[1] eq "Tail"){
		$current_chain=0;
	}
	elsif($align[1] eq "Head"){
		$current_chain=1;
#		$reverse=1;
	}
	my $current_start=0;
	for(my $i=1;$i<@info;$i++){
		$Merged_Contig{$info[$i]}=0;
		my @splice=split /\s+/,$Contig_Pair{$info[$i-1]}{$info[$i]};
		if($splice[1] ne $splice[3]){
			if($current_chain==0){
				if($current_seq eq ""){
					$current_seq=$Contig_Seq{$info[$i-1]};
				}
				else{
					if($reverse==1){
						$current_seq=$Contig_Seq{$info[$i-1]}.$unknown.$current_seq;
					}
					else{
						$current_seq=$current_seq.$unknown.$Contig_Seq{$info[$i-1]};
					}
				}
			}
			elsif($current_chain==1){
				my $seq=reverse_complement($Contig_Seq{$info[$i-1]});
				if($current_seq eq ""){
                                        $current_seq=$seq;
                                }
                                else{
					if($reverse==1){
						$current_seq=$seq.$unknown.$current_seq;
					}
					else{
	                                        $current_seq=$current_seq.$unknown.$seq;
					}
                                }
			}
				
		}
		elsif($splice[1] eq $splice[3]){
				if($current_chain==0){
					if($current_seq eq ""){
                                	        $current_seq=$Contig_Seq{$info[$i-1]};
                                	}
                                	else{
						if($reverse==1){
							$current_seq=$Contig_Seq{$info[$i-1]}.$unknown.$current_seq;
						}
						else{
	                                        	$current_seq=$current_seq.$unknown.$Contig_Seq{$info[$i-1]};
						}
                                	}
					$current_chain=1;
                	        }
                        	elsif($current_chain==1){
                                	my $seq=reverse_complement($Contig_Seq{$info[$i-1]});
					if($current_seq eq ""){
                                                $current_seq=$seq;
                                        }
                                        else{
						if($reverse==1){
							$current_seq=$seq.$unknown.$current_seq;
						}
						else{
	                                                $current_seq=$current_seq.$unknown.$seq;
						}
                                        }
					$current_chain=0;
                	        }
		}
	}
	if($current_chain==0){
		if($current_seq eq ""){
			$current_seq=$Contig_Seq{$info[-1]};
                }
                else{
			if($reverse==1){
				$current_seq=$Contig_Seq{$info[-1]}.$unknown.$current_seq;
			}
			else{
	                        $current_seq=$current_seq.$unknown.$Contig_Seq{$info[-1]};
			}
                }
	}
	elsif($current_chain==1){
		my $seq=reverse_complement($Contig_Seq{$info[-1]});
		if($current_seq eq ""){
                        $current_seq=$seq;
                }
                else{
			if($reverse==1){
				$current_seq=$seq.$unknown.$current_seq;
			}
			else{
	                        $current_seq=$current_seq.$unknown.$seq;
			}
                }
	}
	print OUT "$current_seq\n";
}
foreach my $key (keys %Contig_Seq){
	next if(exists $Merged_Contig{$key});
	print OUT ">$key\n";
	print OUT "$Contig_Seq{$key}\n";
}
close(OUT);

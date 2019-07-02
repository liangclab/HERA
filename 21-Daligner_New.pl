#!/usr/bin/perl
#######################################################################################################
#   This script is used to align the path-contigs to the corresponding gap with Daligner
#                      Author: huilong du
#   Usage: perl $0 Scaffold2Ctg_Gap.txt Prosudo_ScaffoldNonEnzyme2Contig.fasta PathContigRename.fasta
#######################################################################################################
use warnings;
use strict;
my $infile1=shift;            #Scaffold2Ctg_Gap.txt
my $infile2=shift;            #Prosudo_ScaffoldNonEnzyme2Contig.fasta
my $infile3=shift;            #PathContigRename.fasta
my $queue=shift;
my $workspace=shift;          #"shell or qsub"
my $script=shift;
my $genome_name=shift;
my $DAZZ_DB=shift;            #DAZZ_DB
my $DALIGNER=shift;           #DALIGNER

open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;
open IN3,"<$infile3" or die $!;

my %ScaffoldContig=();
my $sign="";
#Super-Scaffold_3.1
while(<IN2>){
	chomp;
	my $line=$_;
	if($line=~/^>(\S+)/){
		$sign=$1;
	}
	else{
		$ScaffoldContig{$sign}=$line;
	}
}
close(IN2);

my %PathContig=();
#Super-Scaffold_262-5-6.1
exit if(! -B $0);
my $gap="";
while(<IN3>){
	chomp;
	my $line=$_;
	if($line=~/^>(\S+)\.(\d+)/){
		$gap=$1;
		$sign=$2;
	}
	else{
		$PathContig{$gap}{$sign}=$line;
	}
}
close(IN3);

while(<IN1>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$content[0]=~/^(\S+)\.(\d+)/;
	my $scaffold=$1;
	my $First=$2;
	my $Second=$First+1;
	my $Gap_Info=$scaffold."-".$First."-".$Second;
	print "$Gap_Info\n";
	my $commond1=`mkdir $Gap_Info`;
        if(! -B $0){open HUI,">$0";print HUI "";exit;}
#	my $commond1=`mkdir $Gap_Info`;
	open OUT,">$Gap_Info/$Gap_Info.fasta" or die $!;
	print OUT ">$content[0]\n";
	print OUT "$ScaffoldContig{$content[0]}\n";
	print OUT ">$content[1]\n";
	print OUT "$ScaffoldContig{$content[1]}\n";
	my $temp=$PathContig{$Gap_Info};
	foreach my $key (keys %$temp){
		next if(length($PathContig{$Gap_Info}{$key})>800000);
		print OUT ">$Gap_Info.$key\n$PathContig{$Gap_Info}{$key}\n";
	}
	close(OUT);
	if($workspace eq "qsub"){
		open OUT,">$Gap_Info.pbs" or die $!;
		print OUT "#BSUB -J $genome_name-DALIGNER-$Gap_Info
#BSUB -o $Gap_Info.out
#BSUB -n 1
#BSUB -q $queue
cd $Gap_Info
perl $script/lines_to_split.pl $Gap_Info.fasta $Gap_Info-formated.fasta
#perl $script/SplitRef.pl $Gap_Info.fasta $Gap_Info-formated.fasta
$DAZZ_DB/fasta2DAM $Gap_Info $Gap_Info-formated.fasta
$DAZZ_DB/DBdust $Gap_Info.dam
$DAZZ_DB/DBsplit -x1000 -s50 $Gap_Info.dam
$DALIGNER/HPC.daligner $Gap_Info.dam > $Gap_Info.sh
time sh $Gap_Info.sh

rm -f $Gap_Info.*.$Gap_Info.*.?*.las
$DAZZ_DB/DBdump -rh $Gap_Info.dam | perl $script/ParseDAZZDB.pl >ParseDAZZDB.txt
cat $Gap_Info*.las >All.las
$DALIGNER/LAdump -cd $Gap_Info.dam All.las | perl $script/ParseLA.pl > $Gap_Info-Final.txt
perl $script/Daligner_Reformate.pl $Gap_Info-Final.txt $Gap_Info-Final_Reformated.txt
";
		close(OUT);
		sleep(3);
		my $commond=`bsub < $Gap_Info.pbs`;
	}
	else{
		open OUT,">$Gap_Info-pipeline.sh" or die $!;
		print OUT "#$Gap_Info
cd $Gap_Info
perl $script/SplitRef.pl $Gap_Info.fasta $Gap_Info-formated.fasta
/public/share/bma/DALIGNER1/DAZZ_DB-master/bin/fasta2DAM $Gap_Info $Gap_Info-formated.fasta
/public/share/bma/DALIGNER1/DAZZ_DB-master/bin/DBdust $Gap_Info.dam
/public/share/bma/DALIGNER1/DAZZ_DB-master/bin/DBsplit -x1000 -s50 $Gap_Info.dam
/public/share/bma/DALIGNER1/DALIGNER-master/bin/HPCdaligner $Gap_Info.dam > $Gap_Info.sh
time sh $Gap_Info.sh

rm -f $Gap_Info.*.$Gap_Info.*.?*.las
/public/share/bma/DALIGNER1/DAZZ_DB-master/bin/DBdump -rh $Gap_Info.dam | perl $script/ParseDAZZDB.pl > ParseDAZZDB.txt
/public/share/bma/DALIGNER1/DALIGNER-master/bin/LAdump -cd $Gap_Info.dam $Gap_Info.las | perl $script/ParseLA.pl > $Gap_Info-Final.txt
perl $script/Daligner_Reformate.pl $Gap_Info-Final.txt $Gap_Info-Final_Reformated.txt
";
		close(OUT);
	}
}
#if($workspace eq "shell"){
#	system("cat Super-Scaffold*-pipeline.sh >Whole-Pipeline.sh");
#}
close(IN1);

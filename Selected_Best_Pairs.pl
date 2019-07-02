#!/usr/bin/perl
use warnings;
use strict;

my $infile=shift;
my $outfile=shift;

open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;

my %Remained_Pairs=();

#left	1	98.62	4474	Pacbio320854END	L6712E_subseq_1:401984_obj/0_401984	left	0	4559	401984	9888	14426	14427	397425	9888
while(<IN>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$content[5]=~/(\S+)\/0_/;
	$content[5]=$1;
	next if($content[9]<6000 || $content[12]<6000);
	if(!exists $Remained_Pairs{$content[4]}{$content[5]} && !exists $Remained_Pairs{$content[5]}{$content[4]}){
		$Remained_Pairs{$content[4]}{$content[5]}{'S'}=$content[3];
		$Remained_Pairs{$content[4]}{$content[5]}{'L'}=$line;
	}
	else{
		if(exists $Remained_Pairs{$content[4]}{$content[5]} && $Remained_Pairs{$content[4]}{$content[5]}{'S'}<$content[3]){
			$Remained_Pairs{$content[4]}{$content[5]}{'S'}=$content[3];
			$Remained_Pairs{$content[4]}{$content[5]}{'L'}=$line;
		}
		elsif(exists $Remained_Pairs{$content[5]}{$content[4]} && $Remained_Pairs{$content[5]}{$content[4]}{'S'}<$content[3]){
			$Remained_Pairs{$content[5]}{$content[4]}{'S'}=$content[3];
                        $Remained_Pairs{$content[5]}{$content[4]}{'L'}=$line;
		}
	}
}
close(IN);

foreach my $key1 (keys %Remained_Pairs){
	my $temp=$Remained_Pairs{$key1};
	foreach my $key2 (keys %$temp){
		print OUT "$Remained_Pairs{$key1}{$key2}{'L'}\n";
	}
}
%Remained_Pairs=();
close(OUT);


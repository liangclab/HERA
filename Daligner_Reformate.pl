#!/usr/bin/perl
use strict;
use warnings;

my $infile=shift;
my $outfile=shift;

open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;


#+       Super-Scaffold_1-2-3.233310     13030   53235   53235   Super-Scaffold_1.3      1       40363   344123  266
#+       Super-Scaffold_1-2-3.233310     0       8811    53235   Super-Scaffold_1.2      959055  968019  968018  216
#
#-       Super-Scaffold_1-2-3.149251     0       21534   34550   Super-Scaffold_1.3      322470  344122  344123  179
#-       Super-Scaffold_1-2-3.149251     25739   34550   34550   Super-Scaffold_1.2      0       8963    968018  216
while(<IN>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	print OUT "$line\n";
	if($content[0] eq "+"){
		print OUT "$content[0]\t$content[5]\t$content[6]\t$content[7]\t$content[8]\t$content[1]\t$content[2]\t$content[3]\t$content[4]\t$content[9]\n";
	}
	elsif($content[0] eq "-"){
		my $query0="-";
		my $query1=$content[5];
		my $query2=$content[8]-$content[7];
		my $query3=$content[8]-$content[6];
		my $query4=$content[8];
		my $query5=$content[1];
		my $query6=$content[4]-$content[3];
		my $query7=$content[4]-$content[2];
		my $query8=$content[4];
		my $query9=$content[9];
		print OUT "$query0\t$query1\t$query2\t$query3\t$query4\t$query5\t$query6\t$query7\t$query8\t$query9\n";
	}
}
close(OUT);
close(IN);

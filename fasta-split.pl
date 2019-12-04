#!/usr/bin/perl

use strict;
use Bio::SeqIO;

my $ntfile = shift(@ARGV) || die "$0 file seq-num [base]\n";
my $seq_num = shift || 1000;
my $base = shift;
my $suffix = shift;

my $stem = $ntfile;
$stem =~ s/^.+\///;
print "$stem\n";
$stem =~ s/\.fasta$|\.fa$//;
$stem = $base || $stem;
print "$stem\n";

my $seqin  = Bio::SeqIO->new(-file => "$ntfile" , '-format' => 'Fasta');
my $seqout;
for (my $i=0; my $seq = $seqin->next_seq(); $i++) {
    if ($i%$seq_num == 0) {
	my $name;
	if ($suffix) {
	    $name = "$stem.".($i/$seq_num+1).".fa";
	}else {
	    $name = "$stem.".($i/$seq_num+1);
	}	    
	$seqout = Bio::SeqIO->new(-file => ">$name" , '-format' => 'Fasta');
    }
    if (!$seq || length($seq->seq)<1) {print $seq->id, "\n";}
    if ($seq->id && $seq->seq && length($seq->seq)>0) {$seqout->write_seq($seq);}
}

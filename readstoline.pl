#!/usr/bin/perl -w
use strict;
use warnings;
my $string="";
my $infile=shift;
my $outfile=shift;
my $sign=shift;
my $Current_Count=1;
open(INPUT,"<$infile")|| die "cannot open file $infile\n";
open(OUTPUT,">$outfile")||die "cannot open file $outfile\n";
while(<INPUT>){
   my($line)=$_;
   chomp($line);
   next if($line=~/^\s+\s+$/);
   if($line=~/^>/){
         if($string ne ""){
               print OUTPUT "$string\n";
         }
         $string="";
         print OUTPUT ">$sign$Current_Count";
	 print OUTPUT "E\n";
	 $Current_Count++;
         next;
   }
#   $string=~/\s+//g;
   $string=$string.$line;
}
print OUTPUT "$string";
close(INPUT);
close(OUTPUT);


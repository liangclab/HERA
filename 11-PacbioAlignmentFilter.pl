#!/usr/bin/perl
#Author: huilong du
#Note:Filtering the alignment for overlaping
use strict;
use warnings;

my $inputfile = shift;
	
my $MinOverhangLength = shift;
my $MinimumIdentity   = shift;
my $MinOverlap        = shift;
my $MinExtend         = shift;

open(IN, $inputfile ) or die("Cannot open $inputfile\n");
open OUT,">overhang_info.txt" or die $!;
while ( my $line = <IN> ) {

	chomp $line; 
    
    next if ( $line =~ m/nCells/ );
    
    # qName	        tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
    # m150611...	contig2_size40580	0	    1	    -29423	89.3895	            14529	21466	40580	614	8012	30328	183070
    my @Columns   = split( /\s+/, $line );
    $Columns[0]=~/(\S+)\/0_/;
    my $First=$1;
    next if($First eq $Columns[1]);
    my $score     = $Columns[4];
    my $percentSimilarity   = $Columns[5];
    
    # Illumina Contig
    my $tName     = $Columns[1];
    my $tStrand   = $Columns[3];
    my $tStart    = $Columns[6];
    my $tEnd      = $Columns[7];
    my $tLength   = $Columns[8];
    # Pacbio
    my $qName     = $Columns[0];
    my $qStrand   = $Columns[2];
    my $qStart    = $Columns[9];
    my $qEnd      = $Columns[10];
    my $qLength   = $Columns[11];
    
    next if ( abs( $qEnd - $qStart ) < $MinOverlap );
    next if ( abs( $tEnd - $tStart ) < $MinOverlap );
    next if ( ($percentSimilarity < $MinimumIdentity));
    
    my $RefLeftOverhang  = 0;
    my $RefRightOverhang = 0;
                            
    my $QryLeftOverhang  = 0 ;
    my $QryRightOverhang = 0 ;
    
    
    if ( $qStrand == 0 && $tStrand == 0 ) {

        $RefLeftOverhang  = $tStart;
        $RefRightOverhang = ( $tLength - $tEnd );
        $QryLeftOverhang  = $qStart;
        $QryRightOverhang = ( $qLength - $qEnd );

        #              ===============>
        #              |||||||||||||||
        #   -------------------------------------------->
        # qName	            tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
        # m.../0_16611	    contig3_size8058	0	    0	    -40268	99.98 	            0	    8058	8058	3994	12052	15952	169190
        if ( ( $RefLeftOverhang <= $MinOverhangLength ) && ( $RefRightOverhang <= $MinOverhangLength ) &&
             ( $QryLeftOverhang >= $MinOverhangLength ) && ( $QryRightOverhang >= $MinOverhangLength ) ) {
             print OUT "$RefLeftOverhang\n$RefRightOverhang\n"
#            print $line."\n";
            
        #            ======================>
        #            |||||||||||||||
        #       ------------------->
        # qName	        tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
        # m.../0_18694	contig9_size44098	0	    0	    -76540	99.8438	            0	    15360	44098	4340	19696	19696	322487
        }elsif ( ( $RefLeftOverhang <= $MinOverhangLength ) && ( $RefRightOverhang >=$MinExtend  ) &&
                ( $QryLeftOverhang >=$MinExtend ) && ( $QryRightOverhang <=$MinOverhangLength  ) ) {
            print OUT "$RefLeftOverhang\n$QryRightOverhang\n";
            print $line."\n";
            
        #         ==========================>
        #                     |||||||||||||||
        #                     --------------------->
        # qName	            tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
        # m.../6331_12271	contig11_size14782	0	    0	    -20685	100.00 	            10645	14782	14782	5	    4142	5640	86843
        }elsif ( ( $RefLeftOverhang >= $MinExtend ) && ( $RefRightOverhang <= $MinOverhangLength ) &&
                ( $QryLeftOverhang <= $MinOverhangLength ) && ( $QryRightOverhang >= $MinExtend ) ) {
            print OUT "$RefRightOverhang\n$QryLeftOverhang\n";
            print $line."\n";
            
        }

    
    
    }elsif ( $qStrand == 0 && $tStrand == 1 ) {

        $RefLeftOverhang  = $tStart;
        $RefRightOverhang = ( $tLength - $tEnd );
        $QryLeftOverhang  = $qStart;
        $QryRightOverhang = ( $qLength - $qEnd );

        #              <==============
        #              |||||||||||||||
        #   -------------------------------------------->
        # qName	            tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
        # m.../32347_45832	contig3_size8058	0	    1	    -39858	98.97 	            0	    8058	8058	2676	10816	13374	177801
        if ( ( $RefLeftOverhang <= $MinOverhangLength ) && ( $RefRightOverhang <= $MinOverhangLength ) &&
             ( $QryLeftOverhang >= $MinOverhangLength ) && ( $QryRightOverhang >= $MinOverhangLength ) ) {
             print OUT "$RefLeftOverhang\n$RefRightOverhang\n";
#            print $line."\n";
            
        #          ======================>
        #          ||||||||||||
        #    <-----------------
        #    <======================
        #            |||||||||||||||          <-----
        #            -------------------->
        #¡¡qName	            tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
        #¡¡m.../32347_45832	    contig11_size14782	0	    1	    -18870	100.00 	            11008	14782	14782	0	    3774	13374	79220
        }elsif ( ( $RefLeftOverhang >= $MinExtend ) && ( $RefRightOverhang <= $MinOverhangLength ) &&
                ( $QryLeftOverhang <= $MinOverhangLength ) && ( $QryRightOverhang >= $MinExtend ) ) {
            print OUT "$RefRightOverhang\n$QryLeftOverhang\n";
            print $line."\n";
            
        #          ======================>
        #                    ||||||||||||
        #                    <----------------------
        #          <======================
        #          |||||||||||||||            <-----
        #    -------------------->
        # qName	            tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
        # m...2558_18452	contig11_size14782	0	    1	    -64185	95.74 	            0	    13496	14782	1282	15296	15297	416624
        }elsif ( ( $RefLeftOverhang <= $MinOverhangLength ) && ( $RefRightOverhang >= $MinExtend ) &&
                ( $QryLeftOverhang >= $MinExtend ) && ( $QryRightOverhang <= $MinOverhangLength ) ) {
            print OUT "$RefLeftOverhang\n$QryRightOverhang\n";
            print $line."\n";
            
        }
        
    }

}

close IN;


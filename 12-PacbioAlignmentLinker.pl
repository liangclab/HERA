#!/usr/bin/perl
#Author: huilong du
#Note: Finding the proper overlap for constructing the graph
use warnings;
use strict;
my $inputfile=shift;
my $MinOverhangLength=shift;
my $MinExtend=shift;
open IN,"<$inputfile" or die $!;
while(my $line=<IN>){
	chomp $line; 
    next if ( $line =~ m/nCells/ );
    
    # qName	        tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
    # m150611...	contig2_size40580	0	    1	    -29423	89.3895	            14529	21466	40580	614	8012	30328	183070
    my @Columns   = split( /\s+/, $line );
   
    my $score     = $Columns[4];
    my $percentSimilarity   = $Columns[5];
    my $Final_Score=0;
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
        if ( (( $RefLeftOverhang <= $MinOverhangLength) && ( $RefRightOverhang <= $MinOverhangLength ) &&
             ( $QryLeftOverhang >= $MinOverhangLength ) && ( $QryRightOverhang >= $MinOverhangLength )) || (( $RefLeftOverhang >= $MinOverhangLength) && ( $RefRightOverhang >= $MinOverhangLength ) &&( $QryLeftOverhang <= $MinOverhangLength ) && ( $QryRightOverhang <= $MinOverhangLength )) ) {
            
#            print "include\tforward\t$percentSimilarity\t$score\t$tName\t$qName\tinclude\t$qStart\t$qEnd\t$qLength\t$tStart	$tEnd	$tLength\n";
            
        #            ======================>
        #            |||||||||||||||
        #       ------------------->
        # qName	        tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
        # m.../0_18694	contig9_size44098	0	    0	    -76540	99.8438	            0	    15360	44098	4340	19696	19696	322487
        }elsif ( (( $RefLeftOverhang <= $MinOverhangLength ) && ( $RefRightOverhang >=$MinExtend  ) &&
                ( $QryLeftOverhang >= $MinExtend ) && ( $QryRightOverhang <= $MinOverhangLength ) )) {
            $Final_Score=int(abs($score)/100-($QryRightOverhang+$RefLeftOverhang)/2);
            print "left\t0\t$percentSimilarity\t$Final_Score\t$tName\t$qName\tright\t$qStart\t$qEnd\t$qLength\t$tStart\t$tEnd\t$tLength\t$QryLeftOverhang\t$RefRightOverhang\n";
            
        #         ==========================>
        #                     |||||||||||||||
        #                     --------------------->
        # qName	            tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
        # m.../6331_12271	contig11_size14782	0	    0	    -20685	100.00 	            10645	14782	14782	5	    4142	5640	86843
        }elsif ( ( $RefLeftOverhang >= $MinExtend ) && ( $RefRightOverhang <= $MinOverhangLength ) &&
                ( $QryLeftOverhang <= $MinOverhangLength ) && ( $QryRightOverhang >= $MinExtend ) ) {
            $Final_Score=int(abs($score)/100-($RefRightOverhang+$QryLeftOverhang)/2);
            print "right\t0\t$percentSimilarity\t$Final_Score\t$tName\t$qName\tleft\t$qStart\t$qEnd\t$qLength\t$tStart\t$tEnd\t$tLength\t$QryRightOverhang\t$RefLeftOverhang\n";
	   
            
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
             ( $QryLeftOverhang >= $MinOverhangLength ) && ( $QryRightOverhang >= $MinOverhangLength ) ||
             (( $RefLeftOverhang >= $MinOverhangLength ) && ( $RefRightOverhang >= $MinOverhangLength ) &&
             ( $QryLeftOverhang <= $MinOverhangLength ) && ( $QryRightOverhang <= $MinOverhangLength ))) {
            
#            print "include\treverse\t$percentSimilarity\t$score\t$tName\t$qName\tinclude\t$qStart\t$qEnd\t$qLength\t$tStart	$tEnd	$tLength\n";
            
        #          ======================>
        #          ||||||||||||
        #    <-----------------
        #    <======================
        #            |||||||||||||||          <-----
        #            -------------------->
        #¡¡qName	            tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
        #¡¡m.../32347_45832	    contig11_size14782	0	    1	    -18870	100.00 	            11008	14782	14782	0	    3774	13374	79220
        }elsif ( ( $RefLeftOverhang >= $MinExtend  ) && ( $RefRightOverhang <= $MinOverhangLength ) &&
                ( $QryLeftOverhang <= $MinOverhangLength ) && ( $QryRightOverhang >= $MinExtend ) ) {
            $Final_Score=int(abs($score)/100-($RefRightOverhang+$QryLeftOverhang)/2);
            print "left\t1\t$percentSimilarity\t$Final_Score\t$tName\t$qName\tleft\t$qStart\t$qEnd\t$qLength\t$tStart\t$tEnd\t$tLength\t$QryRightOverhang\t$RefLeftOverhang\n";
            
        #          ======================>
        #                    ||||||||||||
        #                    <----------------------
        #          <======================
        #          |||||||||||||||            <-----
        #    -------------------->
        # qName	            tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
        # m...2558_18452	contig11_size14782	0	    1	    -64185	95.74 	            0	    13496	14782	1282	15296	15297	416624
        }elsif ( ( $RefLeftOverhang <= $MinOverhangLength ) && ( $RefRightOverhang >= $MinExtend ) &&
                ( $QryLeftOverhang >= $MinExtend  ) && ( $QryRightOverhang <= $MinOverhangLength ) ) {
            $Final_Score=int (abs($score)/100-($RefLeftOverhang+$QryRightOverhang)/2);
            print "right\t1\t$percentSimilarity\t$Final_Score\t$tName\t$qName\tright\t$qStart\t$qEnd\t$qLength\t$tStart\t$tEnd\t$tLength\t$QryLeftOverhang\t$RefRightOverhang\n";
            
        }
        
    }

}

close IN;


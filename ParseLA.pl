use strict;
use warnings;
use Cwd;

my $RefNameFile = 'ParseDAZZDB.txt';
my %ContigNameHash   = ();
my %ContigLengthHash = ();

open( REFFILE, "<$RefNameFile") || die "cannot open file $RefNameFile \n" ;

    # 1	CTG1-1	78097
    # 2	CTG1-2	52763
    
    while(<REFFILE>)
    {
        my($line) = $_;
        chomp($line);
        
        next if ( $line eq '' );
        next if ( $line =~ m/^#/ );
        
        my ( $ReadNum, $ContigName, $ContigLength ) = split ( /\s+/, $line );
        
        $ContigNameHash{$ReadNum}   = $ContigName;
        $ContigLengthHash{$ReadNum} = $ContigLength;
        
    } # end of while

close( REFFILE );


my $ReadNum1    = 0;
my $ReadNum2    = 0;
my $Strand      = '';

my $Read1Start = 0;
my $Read1End   = 0;

my $Read2Start = 0;
my $Read2End   = 0;

my $Diff       = 0;

print "Strand\t";
print "Ref\tRefStart\tRefEnd\tRefLength\t";
print "Qry\tQryStart\tQryEnd\tQryLength\t";
print "Diff\n";

    while(<STDIN>)
    {
        my($line) = $_;
        chomp($line);
        
        # + P 2791
        # % P 65
        # + T 55594
        # % T 1790
        # @ T 1112

        next if ( $line =~ m/^\+/ );
        next if ( $line =~ m/^@/ );
        next if ( $line =~ m/^%/ );
        
        # P 1 22314 c
        # C 0 27038 70434 97489
        # D 62
        # P 3 22127 n
        # C 74110 76280 0 2189
        # D 134

        
        # ReadNum
        if ( $line =~ m/^P/ ) {
            $ReadNum1 = ( split( /\s+/, $line ) )[1];
            $ReadNum2 = ( split( /\s+/, $line ) )[2];
            $Strand   = ( split( /\s+/, $line ) )[3];
            $Strand   = ( $Strand eq 'n' ) ? '+' : '-';
        }elsif ( $line =~ m/^C/ ) {
            $Read1Start = ( split( /\s+/, $line ) )[1];
            $Read1End   = ( split( /\s+/, $line ) )[2];
            $Read2Start = ( split( /\s+/, $line ) )[3];
            $Read2End   = ( split( /\s+/, $line ) )[4];
        }elsif ( $line =~ m/^D/ ) {
            $Diff  = ( split( /\s+/, $line ) )[1];
            print $Strand."\t";
            print $ContigNameHash{$ReadNum1}."\t";
            print $Read1Start."\t";
            print $Read1End."\t";
            print $ContigLengthHash{$ReadNum1}."\t";
            print $ContigNameHash{$ReadNum2}."\t";
            print $Read2Start."\t";
            print $Read2End."\t";
            print $ContigLengthHash{$ReadNum2}."\t";
            print $Diff ."\n";
        }

    } # end of while

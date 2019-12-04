use strict;
use warnings;
use Cwd;

my $ReadNum    = 0;
my $ReadHeader = '';
my $ReadStart = 0;
my $ReadEnd    = 0;
my $ReadLength = 0;

    while(<STDIN>)
    {
        my($line) = $_;
        chomp($line);
        
        # + R 28470
        # + M 0
        # + H 364425
        # @ H 24
        
        next if ( $line =~ m/^\+/ );
        next if ( $line =~ m/^@/ );
        
        # R 1
        # H 7 >CTG1-1
        # L 0 0 78097
        # R 2
        # H 7 >CTG1-2
        # L 0 0 52763
        # R 3
        # H 7 >CTG2-1
        # L 0 0 104307
        
        # ReadNum
        if ( $line =~ m/^R/ ) {
            $ReadNum = ( split( /\s+/, $line ) )[1];
        }elsif ( $line =~ m/^H/ ) {
            $ReadHeader = ( split( /\s+/, $line ) )[2];
            $ReadHeader =~ s/>//;
        }elsif ( $line =~ m/^L/ ) {
            $ReadStart  = ( split( /\s+/, $line ) )[2];
            $ReadEnd    = ( split( /\s+/, $line ) )[3];
            $ReadLength = $ReadEnd - $ReadStart;
            print $ReadNum."\t";
            print $ReadHeader."\t";
            print $ReadLength."\n";
        }

    } # end of while

#!/usr/bin/perl
use strict;
no warnings;
=pod
use Graph;
use Graph::Directed;
use Graph::Undirected;
use Graph::Easy;
use Bio::SeqIO;
use Bio::Perl;
use Data::Dumper;
=cut

my $DEBUG=0;
my $inputfile=shift;                #blasr_result output
#my $ContigNameListFile=shift;       #all contig name       contig_name.list
#my $MergedContigFile = shift;       #output graph png 
#my $out_cluster_file= shift;        #output ctg clusters
my @ColumnNameArray = qw ( RefPos RefStrand Identity Score RefName QryName QryPos QryStart QryEnd QryLength RefStart RefEnd RefLength );

#my $Pair_Graph = Graph::Undirected->new;
my %Contig_Head_Tail_Query = % {( &ReadingStartPoint)[0]} ;
my %Contig_Gap_Length=%{(&ReadingGapInfo)[0]};
my %Reads_Pairs=();
#my ( $Finding_Best_Graph , $Output_Graph) = &initializeGraph( $ContigNameListFile );

#my %NodeHash = %{ ( &GenerateRecordHash( $inputfile ) )[0] };
#my %EdgeHash = %{ ( &GenerateRecordHash( $inputfile  ))[1] };
#my %NodePairsHash = %{ ( &GenerateNodePairHash( $inputfile ))[0] };
#&drawEasyGraph ( $EasyGraph , $MergedContigFile , '01' );


############ 整合所有节点和边的信息，合并边，并且选择最好的pacbio作为其连接，并且计算两个NODE连接边的权重 #######
my ($ctg_ctg_count_point, $ctg_ctg_overlap_len_point, $ctg_ctg_ori_point, $ctg_ctg_chain_point,$ctg_ctg_line_point,$ctg_ctg_extend_point,$ctg_ctg_info_point)=&GenerateNodePairHash($inputfile);
my %ctg_ctg_count=%$ctg_ctg_count_point;
my %ctg_ctg_overlap_len=%$ctg_ctg_overlap_len_point;
my %ctg_ctg_ori=%$ctg_ctg_ori_point;
my %ctg_ctg_chain=%$ctg_ctg_chain_point;
my %ctg_ctg_line=%$ctg_ctg_line_point;
my %ctg_ctg_extend=%$ctg_ctg_extend_point;
my %ctg_ctg_info=%$ctg_ctg_info_point;


=pod

############ 过滤图中的边，如果一个节点的度>=3，则选择不同方向的两条边 ##########################################
%ctg_ctg_count = %{ (&Filtering_Edge (\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line))[0])};
%ctg_ctg_overlap_len=%{ (&Filtering_Edge (\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line))[1])};
%ctg_ctg_ori = %{ (&Filtering_Edge (\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line))[2])};
%ctg_ctg_chain = %{ (&Filtering_Edge (\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line))[3])};
%ctg_ctg_line = %{ (&Filtering_Edge (\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line))[4])};
########### 结束 ################################################################################################
=cut

############ 根据节点和边，将每对contig分别进行延伸，每次选择最好的边  ##########################################
&Finding_Best_Pathway(\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line,\%Contig_Head_Tail_Query);
&Finding_MaxExtending_Pathway(\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line,\%Contig_Head_Tail_Query, \%ctg_ctg_extend, \%ctg_ctg_info);
&Finding_Randoming_Pathway(\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line,\%Contig_Head_Tail_Query, \%ctg_ctg_extend, \%ctg_ctg_info);


############ 根据Scaffold中Contig开头和末尾的Query分别进行延伸 ##################################################


############ End ################################################################################################

if($DEBUG==1){
#    print Dumper(%NodePairsHash);
    print Dumper(%ctg_ctg_count);
    print Dumper(%ctg_ctg_overlap_len);
    print Dumper(%ctg_ctg_ori);
    print Dumper(%ctg_ctg_chain);
    print Dumper(%ctg_ctg_line);
}

############## 读取Contig之间Gap的长度和信息  #####################################################################
sub ReadingGapInfo{
    my %Contig_Gap=();
    open GAP,"<Scaffold2Ctg_Gap.txt" or die $!;
#Super-Scaffold_5.5      Super-Scaffold_5.6      499     3029936 3030434
    while(<GAP>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$Contig_Gap{$content[0]}{$content[1]}=$content[2];
    }
    close(GAP);
    return(\%Contig_Gap);
}
#############  Contig之间Gap的长度和信息    ######################################################################

############## 读取Contig末端的Query，作为Graph 遍历的起点，以及用于判断两个节点间是否存在通路 ####################

sub ReadingStartPoint{
    my %Scaffold_Contig_SEPoint=();
    open Point,"<Contig_Head_Tail_Pacbio_Pos.txt" or die $!;
#R498_3390190_1/90_20012/0_19922 1 19922 19922 Super-Scaffold_14.29 9130 29199 963453 0 99.13
    while(<Point>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$content[0]=~/(\S+)\/0_\d+$/;
	my $pacbio=$1;
	$content[4]=~/(\S+)\.(\d+)$/;
	my $scaffold=$1;
	my $scaffold_sign=$2;
	if($content[7]-$content[6]>=$content[5]){
		$Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Head'}{$pacbio}{0}=$content[8]; #Chain
		$Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Head'}{$pacbio}{1}=$content[5]; #Ref_start	
		$Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Head'}{$pacbio}{2}=$content[6]; #Ref_end
		$Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Head'}{$pacbio}{3}=$content[7]; #Ref_len
		$Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Head'}{$pacbio}{4}=$content[1]; #Qry_start
		$Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Head'}{$pacbio}{5}=$content[2]; #Qry_end
		$Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Head'}{$pacbio}{6}=$content[3]; #Qry_len
	}
	else{
		$Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Tail'}{$pacbio}{0}=$content[8];
                $Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Tail'}{$pacbio}{1}=$content[5];
                $Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Tail'}{$pacbio}{2}=$content[6];
                $Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Tail'}{$pacbio}{3}=$content[7];
                $Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Tail'}{$pacbio}{4}=$content[1];
                $Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Tail'}{$pacbio}{5}=$content[2];
                $Scaffold_Contig_SEPoint{$scaffold}{$scaffold_sign}{'Tail'}{$pacbio}{6}=$content[3];
	}
    }
    close(Point);
    return(\%Scaffold_Contig_SEPoint);
}


############## 统计支持ctg对的pacbio数，以及选择overlap最长的pacbio作为其连接两个ctg的信息 ########################
sub count_ctg_pair{
    my $CtgPairsHash = shift;
    my %NodePair=%{ $CtgPairsHash };
    my %ctg_ctg_count=();
    my %ctg_ctg_overlap_len=();
    my %ctg_ctg_ori=();
    my %ctg_ctg_chain=();
    my %ctg_ctg_line=();
    open PAIR,">ctg_pairs.txt" or die $!;
    open ORI,">ctg_ctg_ori.txt" or die $!;
    foreach my $key1 (keys %NodePair){
        my %CtgPairPacbio= %{ $NodePair{$key1} };
        my @content=split "-",$key1;
        my $left_ctg=$content[0];
        my $right_ctg=$content[1];
        my $best_pacbio="";
        foreach my $pacbio (keys %CtgPairPacbio){
             $ctg_ctg_count{$left_ctg}{$right_ctg}++;
             my $overlap_len=$CtgPairPacbio{$pacbio}{'First'}{'Score'};
             if(!exists $ctg_ctg_overlap_len{$left_ctg}{$right_ctg}){
                     $ctg_ctg_overlap_len{$left_ctg}{$right_ctg}=0;
             }
             if($ctg_ctg_overlap_len{$left_ctg}{$right_ctg}<$overlap_len){
                     $ctg_ctg_overlap_len{$left_ctg}{$right_ctg}=$overlap_len;
                     $ctg_ctg_ori{$left_ctg}{$right_ctg}{1}=$CtgPairPacbio{$pacbio}{'First'}{'RefPos'};
                     $ctg_ctg_ori{$left_ctg}{$right_ctg}{2}=$CtgPairPacbio{$pacbio}{'Next'}{'RefPos'};
                     $ctg_ctg_ori{$right_ctg}{$left_ctg}{1}=$CtgPairPacbio{$pacbio}{'Next'}{'RefPos'};
                     $ctg_ctg_ori{$right_ctg}{$left_ctg}{2}=$CtgPairPacbio{$pacbio}{'First'}{'RefPos'};
                     
                     $ctg_ctg_chain{$left_ctg}{$right_ctg}{1}=$CtgPairPacbio{$pacbio}{'First'}{'RefStrand'};
                     $ctg_ctg_chain{$left_ctg}{$right_ctg}{2}=$CtgPairPacbio{$pacbio}{'Next'}{'RefStrand'};
                     $ctg_ctg_chain{$right_ctg}{$left_ctg}{1}=$CtgPairPacbio{$pacbio}{'Next'}{'RefStrand'};
                     $ctg_ctg_chain{$right_ctg}{$left_ctg}{2}=$CtgPairPacbio{$pacbio}{'First'}{'RefStrand'};

                     $ctg_ctg_line{$left_ctg}{$right_ctg}=$left_ctg."-".$right_ctg."-".$pacbio;
                     $ctg_ctg_line{$right_ctg}{$left_ctg}=$right_ctg."-".$left_ctg."-".$pacbio;
                     $best_pacbio=$pacbio;
             }
        }
        my $First_start=0;
        my $First_end=0;
        my $Next_start=0;
        my $Next_end=0;
        print ORI "$left_ctg\t$CtgPairPacbio{$best_pacbio}{'First'}{'RefLength'}\t$ctg_ctg_chain{$left_ctg}{$right_ctg}{1}\t$ctg_ctg_ori{$left_ctg}{$right_ctg}{1}\t$right_ctg\t$CtgPairPacbio{$best_pacbio}{'Next'}{'RefLength'}\t$ctg_ctg_chain{$left_ctg}{$right_ctg}{2}\t$ctg_ctg_ori{$left_ctg}{$right_ctg}{2}\t$ctg_ctg_count{$left_ctg}{$right_ctg}\n";
        
        if($CtgPairPacbio{$best_pacbio}{'First'}{'RefStrand'}==0){
             $First_start=$CtgPairPacbio{$best_pacbio}{'First'}{'RefStart'};
             $First_end=$CtgPairPacbio{$best_pacbio}{'First'}{'RefEnd'};
        }
        elsif($CtgPairPacbio{$best_pacbio}{'First'}{'RefStrand'}==1){
             $First_start=$CtgPairPacbio{$best_pacbio}{'First'}{'RefLength'}-$CtgPairPacbio{$best_pacbio}{'First'}{'RefEnd'};
             $First_end=$CtgPairPacbio{$best_pacbio}{'First'}{'RefLength'}-$CtgPairPacbio{$best_pacbio}{'First'}{'RefStart'};
        }
        if($CtgPairPacbio{$best_pacbio}{'Next'}{'RefStrand'}==0){
             $Next_start=$CtgPairPacbio{$best_pacbio}{'Next'}{'RefStart'};
             $Next_end=$CtgPairPacbio{$best_pacbio}{'Next'}{'RefEnd'};
        }
        elsif($CtgPairPacbio{$best_pacbio}{'Next'}{'RefStrand'}==1){
             $Next_start=$CtgPairPacbio{$best_pacbio}{'Next'}{'RefLength'}-$CtgPairPacbio{$best_pacbio}{'Next'}{'RefEnd'};
             $Next_end=$CtgPairPacbio{$best_pacbio}{'Next'}{'RefLength'}-$CtgPairPacbio{$best_pacbio}{'Next'}{'RefStart'};
        }
        print PAIR "$left_ctg\t$CtgPairPacbio{$best_pacbio}{'First'}{'RefLength'}\t$right_ctg\t$CtgPairPacbio{$best_pacbio}{'Next'}{'RefLength'}\t$ctg_ctg_count{$left_ctg}{$right_ctg}\t$best_pacbio-$ctg_ctg_chain{$left_ctg}{$right_ctg}{1}-$CtgPairPacbio{$best_pacbio}{'First'}{'RefStart'}-$CtgPairPacbio{$best_pacbio}{'First'}{'RefEnd'}-$ctg_ctg_chain{$left_ctg}{$right_ctg}{2}-$CtgPairPacbio{$best_pacbio}{'Next'}{'RefStart'}-$CtgPairPacbio{$best_pacbio}{'Next'}{'RefEnd'}\t$best_pacbio-$CtgPairPacbio{$best_pacbio}{'First'}{'QryLength'}-$CtgPairPacbio{$best_pacbio}{'First'}{'QryStart'}-$CtgPairPacbio{$best_pacbio}{'First'}{'QryEnd'}-$CtgPairPacbio{$best_pacbio}{'Next'}{'QryStart'}-$CtgPairPacbio{$best_pacbio}{'Next'}{'QryEnd'}\n";        
    }
    return(\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line);
    close(ORI);
    close(PAIR);
}

sub initializeGraph{

    my $ContigNameListFile = shift;
    my $ContigGraph = Graph::Undirected->new;
    my $EasyGraph   = Graph::Easy->new( undirected => 1 );

    open(IN, $ContigNameListFile ) or die("Cannot open $ContigNameListFile\n");

    while ( my $line = <IN> ) {

        chomp $line;

        # contig0001_size5554
        my $Contig = $line;
#        $ContigGraph -> add_vertex( $Contig );
#        $EasyGraph   -> add_node( $Contig );
        
    }      
    return $ContigGraph, $EasyGraph ;
}
sub drawEasyGraph{

    my $EasyGraph        = shift;
    my $MergedContigFile = shift;
    my $SubNum           = shift;
    my $graphviz = $EasyGraph -> as_graphviz();
#    my $BACName = ( $MergedContigFile =~ m/TUG.(\S+).pacmerge(\d+)./ )[0];
#    my $PacmergeNum = ( $MergedContigFile =~ m/TUG.(\S+).pacmerge(\d+)./ )[1];
#    open my $DOT, "|dot -Tpng -o $MergedContigFile.png" or die ("Cannot open pipe to dot: $!");
#    print $DOT $graphviz;
#    close $DOT;

}
sub GenerateRecordHash{

    my $inputfile = shift;
    my %NodeHash = ();
    my %EdgeHash = ();
    my %NodePairsHash = ();
    open(IN, $inputfile ) or die("Cannot open $inputfile\n");
    while ( my $line = <IN> ) {
              chomp $line;
              my @Columns   = split( /\s+/, $line );
              my $RefName   = $Columns[4];
              my $QryName   = $Columns[5];
              my $Score     = $Columns[3];

              if ( ! exists ( $NodeHash{$RefName}{$QryName} ) ) {
                     for ( my $Idx = 0 ; $Idx < scalar( @Columns ) ; $Idx ++ ) {
                             my $ColumnName = $ColumnNameArray[$Idx];
                             $NodeHash{$RefName}{$QryName}{$ColumnName} = $Columns[$Idx];
                     }
              }else{
                      if ( $Score < $NodeHash{$RefName}{$QryName}{'Score'} ) {

                              for ( my $Idx = 0 ; $Idx < scalar( @Columns ) ; $Idx ++ ) {
                                      my $ColumnName = $ColumnNameArray[$Idx];
                                      $NodeHash{$RefName}{$QryName}{$ColumnName} = $Columns[$Idx];
                              }
                       }
              } 
              if ( ! exists ( $EdgeHash{$QryName}{$RefName} ) ) {
                       for ( my $Idx = 0 ; $Idx < scalar( @Columns ) ; $Idx ++ ) {
                            my $ColumnName = $ColumnNameArray[$Idx];
                            $EdgeHash{$QryName}{$RefName}{$ColumnName} = $Columns[$Idx];
                       }
              }else{
                        if ( $Score < $EdgeHash{$QryName}{$RefName}{'Score'} ) {
                            for ( my $Idx = 0 ; $Idx < scalar( @Columns ) ; $Idx ++ ) {
                                 my $ColumnName = $ColumnNameArray[$Idx];
                                 $EdgeHash{$QryName}{$RefName}{$ColumnName} = $Columns[$Idx];
                            }
                        }
              }

      }
      close IN;

      return ( \%NodeHash, \%EdgeHash );
}

sub GenerateNodePairHash{
    my $file=shift;
    open IN1,"<$file" or die $!;
    open PAIR,">ctg_pairs.txt" or die $!;
    open ORI,">ctg_ctg_ori.txt" or die $!;
    my %NodePairsHash = ();
#left 1 99.43 13111 R498_5995724_1/67_17121 R498_6008174_1/88_16840/0_16752 left 0 13208 16752 3863 17051 17054 QryExtend RefExtend 
    my %ctg_ctg_count=();
    my %ctg_ctg_overlap_len=();
    my %ctg_ctg_ori=();
    my %ctg_ctg_chain=();
    my %ctg_ctg_line=();
    my %ctg_ctg_info=();
    my %ctg_ctg_extend=();
    while(<IN1>){
	chomp;
	my $line=$_;
	my @content=split /\s+/,$line;
	$content[5]=~/(\S+)\/0_\d+$/;
	my $next=$1;
	$content[5]=$next;
	next if($next eq $content[4]);
	my $QryName=1;
#RefPos RefStrand Identity Score RefName QryName QryPos QryStart QryEnd QryLength RefStart RefEnd RefLength QryExtend RefExtend
#        next if(exists $ctg_ctg_info{ $QryStartHashArray[0]}{$QryStartHashArray[1]} && $ctg_ctg_info{ $QryStartHashArray[0]}{$QryStartHashArray[1]}{'overlap'}>=$content[3]);

		
                $Reads_Pairs{$content[4]}{$content[0]}{$content[5]}{'S'}=$content[3];           #Score
                $Reads_Pairs{$content[4]}{$content[0]}{$content[5]}{'O'}=$content[6];           #Oritention of $content[5]
                $Reads_Pairs{$content[4]}{$content[0]}{$content[5]}{'E'}=$content[13];          #Extend length for $content[5]
                $Reads_Pairs{$content[4]}{$content[0]}{$content[5]}{'I'}=$content[2];           #Identity of overlap
                $Reads_Pairs{$content[4]}{$content[0]}{$content[5]}{'C'}=0;                     #chain of $content[5]

                $Reads_Pairs{$content[5]}{$content[6]}{$content[4]}{'S'}=$content[3];           #Score
                $Reads_Pairs{$content[5]}{$content[6]}{$content[4]}{'O'}=$content[0];           #Oritention of $content[4]
                $Reads_Pairs{$content[5]}{$content[6]}{$content[4]}{'E'}=$content[14];          #Extend length for $content[4]
                $Reads_Pairs{$content[5]}{$content[6]}{$content[4]}{'I'}=$content[2];           #Identity of overlap
                $Reads_Pairs{$content[5]}{$content[6]}{$content[4]}{'C'}=$content[1];           #chain of $content[4]



#                $ctg_ctg_extend{$QryStartHashArray[0]}{$QryStartHashArray[1]}{1}=$content[14];
#		$ctg_ctg_extend{$QryStartHashArray[0]}{$QryStartHashArray[1]}{2}=$content[13];
#                $ctg_ctg_extend{$QryStartHashArray[1]}{$QryStartHashArray[0]}{1}=$content[13];
#                $ctg_ctg_extend{$QryStartHashArray[1]}{$QryStartHashArray[0]}{2}=$content[14];
#		$ctg_ctg_count{$QryStartHashArray[0]}{$QryStartHashArray[1]}=$content[3];
#		$ctg_ctg_count{$QryStartHashArray[1]}{$QryStartHashArray[0]}=$content[3];
#		$ctg_ctg_overlap_len{$QryStartHashArray[0]}{$QryStartHashArray[1]}=$content[3];
#		$ctg_ctg_ori{$QryStartHashArray[0]}{$QryStartHashArray[1]}{1}=$content[0];
#		$ctg_ctg_ori{$QryStartHashArray[0]}{$QryStartHashArray[1]}{2}=$content[6];
#		$ctg_ctg_ori{$QryStartHashArray[1]}{$QryStartHashArray[0]}{1}=$content[6];
#                $ctg_ctg_ori{$QryStartHashArray[1]}{$QryStartHashArray[0]}{2}=$content[0];

#		$ctg_ctg_chain{$QryStartHashArray[0]}{$QryStartHashArray[1]}{1}=$content[1];
#		$ctg_ctg_chain{$QryStartHashArray[0]}{$QryStartHashArray[1]}{2}=0;
#		$ctg_ctg_chain{$QryStartHashArray[1]}{$QryStartHashArray[0]}{1}=0;
#		$ctg_ctg_chain{$QryStartHashArray[1]}{$QryStartHashArray[0]}{2}=$content[1];

#		$ctg_ctg_line{$QryStartHashArray[0]}{$QryStartHashArray[1]}=$QryStartHashArray[0]."-".$QryStartHashArray[1]."-"."1";
#		$ctg_ctg_line{$QryStartHashArray[1]}{$QryStartHashArray[0]}=$QryStartHashArray[1]."-".$QryStartHashArray[0]."-"."1";

#		$ctg_ctg_info{$QryStartHashArray[0]}{$QryStartHashArray[1]}{1}{'length'}=$content[12];
#		$ctg_ctg_info{$QryStartHashArray[0]}{$QryStartHashArray[1]}{1}{'start'}=$content[10];
#		$ctg_ctg_info{$QryStartHashArray[0]}{$QryStartHashArray[1]}{1}{'end'}=$content[11];
#               $ctg_ctg_info{$QryStartHashArray[1]}{$QryStartHashArray[0]}{2}{'length'}=$content[12];
#                $ctg_ctg_info{$QryStartHashArray[1]}{$QryStartHashArray[0]}{2}{'start'}=$content[10];
#                $ctg_ctg_info{$QryStartHashArray[1]}{$QryStartHashArray[0]}{2}{'end'}=$content[11];

#                $ctg_ctg_info{$QryStartHashArray[0]}{$QryStartHashArray[1]}{'identity'}=$content[2];
#                $ctg_ctg_info{$QryStartHashArray[0]}{$QryStartHashArray[1]}{'overlap'}=$content[3];
#                $ctg_ctg_info{$QryStartHashArray[1]}{$QryStartHashArray[0]}{'identity'}=$content[2];
#                $ctg_ctg_info{$QryStartHashArray[1]}{$QryStartHashArray[0]}{'overlap'}=$content[3];

#		$ctg_ctg_info{$QryStartHashArray[0]}{$QryStartHashArray[1]}{2}{'length'}=$content[9];
#		$ctg_ctg_info{$QryStartHashArray[0]}{$QryStartHashArray[1]}{2}{'start'}=$content[7];
#		$ctg_ctg_info{$QryStartHashArray[0]}{$QryStartHashArray[1]}{2}{'end'}=$content[8];
#                $ctg_ctg_info{$QryStartHashArray[1]}{$QryStartHashArray[0]}{1}{'length'}=$content[9];
#                $ctg_ctg_info{$QryStartHashArray[1]}{$QryStartHashArray[0]}{1}{'start'}=$content[7];
#                $ctg_ctg_info{$QryStartHashArray[1]}{$QryStartHashArray[0]}{1}{'end'}=$content[8];
		 print ORI "$content[4]\t$content[12]\t$content[1]\t$content[0]\t$content[5]\t$content[9]\t0\t$content[6]\t$content[3]\n";
                print ORI "$content[5]\t$content[9]\t0\t$content[6]\t$content[4]\t$content[12]\t$content[1]\t$content[0]\t$content[3]\n";

                print PAIR "$content[4]\t$content[12]\t$content[5]\t$content[9]\t$content[3]\tNA/0_1-$content[1]-$content[10]-$content[11]-0-$content[7]-$content[8]\tNA/0_1-0-0-0-0-0\n";
                print PAIR "$content[5]\t$content[9]\t$content[4]\t$content[12]\t$content[3]\tNA/0_1-0-$content[7]-$content[8]-$content[1]-$content[10]-$content[11]\tNA/0_1-0-0-0-0-0\n";
#RefPos RefStrand Identity Score RefName QryName QryPos QryStart QryEnd QryLength RefStart RefEnd RefLength
    }  
#NA/0_1-1---0--
=pod
    foreach my $key1 (keys %ctg_ctg_info){
	my $temp=$ctg_ctg_info{$key1};
	foreach my $key2 (keys %$temp){
		print ORI "$key1\t$ctg_ctg_info{$key1}{$key2}{1}{'length'}\t$ctg_ctg_chain{$key1}{$key2}{1}\t$ctg_ctg_ori{$key1}{$key2}{1}\t$key2\t$ctg_ctg_info{$key1}{$key2}{2}{'length'}\t$ctg_ctg_chain{$key1}{$key2}{2}\t$ctg_ctg_ori{$key1}{$key2}{2}\t$ctg_ctg_count{$key1}{$key2}\n";
		print PAIR "$key1\t$ctg_ctg_info{$key1}{$key2}{1}{'length'}\t$key2\t$ctg_ctg_info{$key1}{$key2}{2}{'length'}\t$ctg_ctg_count{$key1}{$key2}\tNA/0_1-$ctg_ctg_chain{$key1}{$key2}{1}-$ctg_ctg_info{$key1}{$key2}{1}{'start'}-$ctg_ctg_info{$key1}{$key2}{1}{'end'}-$ctg_ctg_chain{$key1}{$key2}{2}-$ctg_ctg_info{$key1}{$key2}{2}{'start'}-$ctg_ctg_info{$key1}{$key2}{2}{'end'}\tNA/0_1-0-0-0-0-0\n";
	}
    }
=cut
    close(PAIR);
    close(ORI);
return(\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line, \%ctg_ctg_extend, \%ctg_ctg_info);
}


#&Finding_best_pathway(\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line, \%NodePairsHash );
sub Finding_best_pathway{
    my $ctg_pair_count = shift;
    my $ctg_ctg_overlap = shift;
    my $ctg_ori = shift;
    my $ctg_chain = shift;
    my $ctg_line =shift ;
#    my $PairsHash = shift;
    my %ctg_ctg_count = %{$ctg_pair_count};
    my %ctg_ctg_overlap = %{$ctg_ctg_overlap};
    my %ctg_ctg_ori = %{$ctg_ori};
    my %ctg_ctg_chain = %{$ctg_chain};
    my %ctg_ctg_line = %{$ctg_line};
#    my %NodePairsHash = %{$PairsHash};
    my %All_existed_clusters=();         #$All_existed_clusters
    my @cluster_line=();                 #ctg clusters  
    open OUT,">ctg_clusters.txt" or die $!;
    #首先判断一对contig中，是否有一个contig是最好的，如果两个contig都不是最好的，过；否则延伸

    my %existed_pairs=();
    my %all_existed_ctgs=();
    my %all_existed_pairs=();
    foreach my $left_ctg (keys %ctg_ctg_line){
          my $temp_ctg=$ctg_ctg_line{$left_ctg};
          foreach my $right_ctg (keys %$temp_ctg){
             next if(exists $existed_pairs{$left_ctg}{$right_ctg});
             my $left_ctg_order=1;
             my $right_ctg_order=1;
	     
             ################ 判断左边contig是否为其同方向中最好的 ##############################################
=pod
             foreach my $left_key1 (keys %ctg_ctg_line){
                 my $left_temp_ctg=$ctg_ctg_line{$left_key1};
                 foreach my $left_key2 (keys %$left_temp_ctg){
                     next if($left_ctg eq $left_key1 && $right_ctg eq $left_key2);
                     if($left_ctg eq $left_key2){
                         if($ctg_ctg_ori{$left_key1}{$left_key2}{2} eq $ctg_ctg_ori{$left_ctg}{$right_ctg}{1}){
                              if($ctg_ctg_count{$left_ctg}{$right_ctg}<$ctg_ctg_count{$left_key1}{$left_key2}){
                                  $left_ctg_order++;
                              }
                              elsif($ctg_ctg_count{$left_ctg}{$right_ctg}==$ctg_ctg_count{$left_key1}{$left_key2}){
                                 if($ctg_ctg_overlap{$left_ctg}{$right_ctg}<$ctg_ctg_overlap{$left_key1}{$left_key2}){
                                     $left_ctg_order++;
                                 }
                              }
                         }
                     }
                 }
              }
=cut
              my $left_add=$ctg_ctg_line{$left_ctg};
              foreach my $left_key (keys %$left_add){
                    next if($left_key eq $right_ctg);
                    if($ctg_ctg_ori{$left_ctg}{$left_key}{1} eq $ctg_ctg_ori{$left_ctg}{$right_ctg}{1}){
                       if($ctg_ctg_count{$left_ctg}{$right_ctg}<$ctg_ctg_count{$left_ctg}{$left_key}){
                           $left_ctg_order++;
                       }
                       elsif($ctg_ctg_count{$left_ctg}{$right_ctg}==$ctg_ctg_count{$left_ctg}{$left_key}){
                           if($ctg_ctg_overlap{$left_ctg}{$right_ctg}<$ctg_ctg_overlap{$left_ctg}{$left_key}){
                               $left_ctg_order++;
                           }
                       }
                    }
               }
               ##############判断右边contig是否为其同方向中最好的contig #########################################
=pod
               foreach my $right_key1 (keys %ctg_ctg_line){
                     my $right_temp_ctg=$ctg_ctg_line{$right_key1};
                     foreach my $right_key2 (keys %$right_temp_ctg){
                         next if($left_ctg eq $right_key1 && $right_ctg eq $right_key2);
                         if($right_ctg eq $right_key2){
                             if($ctg_ctg_ori{$right_key1}{$right_key2}{2} eq $ctg_ctg_ori{$left_ctg}{$right_ctg}{2}){
                                 if($ctg_ctg_count{$left_ctg}{$right_ctg}<$ctg_ctg_count{$right_key1}{$right_key2}){
                                    $right_ctg_order++;
                                 }
                                 elsif($ctg_ctg_count{$left_ctg}{$right_ctg}==$ctg_ctg_count{$right_key1}{$right_key2}){
                                     if($ctg_ctg_overlap{$left_ctg}{$right_ctg}<$ctg_ctg_overlap{$right_key1}{$right_key2}){
                                         $right_ctg_order++;
                                     }
                                 }
                              }
                          }
                      }
                }
=cut
                my $right_add=$ctg_ctg_line{$right_ctg};
                foreach my $right_key (keys %$right_add){
                      next if($right_key eq $left_ctg);
                      if($ctg_ctg_ori{$right_ctg}{$right_key}{1} eq $ctg_ctg_ori{$left_ctg}{$right_ctg}{2}){
                           if($ctg_ctg_count{$left_ctg}{$right_ctg}<$ctg_ctg_count{$right_ctg}{$right_key}){
                                  $right_ctg_order++;
                           }
                           elsif($ctg_ctg_count{$left_ctg}{$right_ctg}==$ctg_ctg_count{$right_ctg}{$right_key}){
                              if($ctg_ctg_overlap{$left_ctg}{$right_ctg}<$ctg_ctg_overlap{$right_ctg}{$right_key}){
                                  $right_ctg_order++;
                              }
                           }
                       }
                }
                next if($right_ctg_order>1 || $left_ctg_order>1);
                next if(exists $all_existed_pairs{$left_ctg}{$right_ctg});
                next if(exists $all_existed_ctgs{$left_ctg} || exists $all_existed_ctgs{$right_ctg});
                $all_existed_pairs{$left_ctg}{$right_ctg}=0;
                $all_existed_ctgs{$left_ctg}=0;
                $all_existed_ctgs{$right_ctg}=0;
#                print "$left_ctg\t$right_ctg\n";
                #################### 结束判断 ################################################################
         
                ################### 以上面找到的一对contig为起始 分别左右延伸寻找最好的contig ################ 
                my @extend_line=();
                push(@extend_line,$left_ctg);
                push(@extend_line,$right_ctg);
                my $left_extend=$left_ctg;
                my $right_extend=$right_ctg;

                my %existed_ctg=();
                $existed_ctg{$left_extend}=0;
                $existed_ctg{$right_extend}=0;
         
                my %extended_ctg=();
                my %existed_ctg_ori=();
                $existed_ctg_ori{$left_ctg}=$ctg_ctg_ori{$left_ctg}{$right_ctg}{1};
                $existed_ctg_ori{$right_ctg}=$ctg_ctg_ori{$left_ctg}{$right_ctg}{2};         
                my $left_extend_ctg="1";
                my $left_extend_count=0;
                my $left_overlap=0;
                while($left_extend_ctg ne ""){
                    $left_extend_ctg="";
                    $left_extend_count=0;
                    $left_overlap=0;
                    $left_extend=$extend_line[0];
#                    last if(exists $extended_ctg{$left_extend});      #判断该contig是否已经延伸过
                    $extended_ctg{$left_extend}=0;
                    foreach my $key1 (keys %ctg_ctg_line){
                          my $temp=$ctg_ctg_line{$key1};
                          foreach my $key2 (keys %$temp){
                              if($key2 eq $left_extend){
                                  if($ctg_ctg_ori{$key1}{$key2}{2} ne $ctg_ctg_ori{$left_extend}{$extend_line[1]}{1}){
                                      if($ctg_ctg_count{$key1}{$key2}>$left_extend_count){
                                          $left_extend_ctg=$key1;
                                          $left_extend_count=$ctg_ctg_count{$key1}{$key2};
                                          $left_overlap=$ctg_ctg_overlap{$key1}{$key2};
                                      }
                                      elsif($ctg_ctg_count{$key1}{$key2}==$left_extend_count){
                                          if($ctg_ctg_overlap{$key1}{$key2}>$left_overlap){
                                               $left_extend_ctg=$key1;
                                               $left_extend_count=$ctg_ctg_count{$key1}{$key2};
                                               $left_overlap=$ctg_ctg_overlap{$key1}{$key2};
                                           }
                                      }
                                 }
                              }
                           }
                     }
                     my $temp=$ctg_ctg_line{$left_extend};
                     foreach my $key10 (keys %$temp ){
                         if($ctg_ctg_ori{$left_extend}{$key10}{1} ne $ctg_ctg_ori{$left_extend}{$extend_line[1]}{1}){
                             if($ctg_ctg_count{$left_extend}{$key10}>$left_extend_count){
                                  $left_extend_ctg=$key10;
                                  $left_extend_count=$ctg_ctg_count{$left_extend}{$key10};
                                  $left_overlap=$ctg_ctg_overlap{$left_extend}{$key10};
                             }
                             elsif($ctg_ctg_count{$left_extend}{$key10}==$left_extend_count){
                                  if($ctg_ctg_overlap{$left_extend}{$key10}>$left_overlap){
                                     $left_extend_ctg=$key10;
                                     $left_extend_count=$ctg_ctg_count{$left_extend}{$key10};
                                     $left_overlap=$ctg_ctg_overlap{$left_extend}{$key10};
                                  }
                              }
                          }
                      }
                      if(exists $all_existed_ctgs{$left_extend_ctg}){
                          last;
                      }
                      if(exists $existed_ctg{$left_extend_ctg}){
                          $left_extend_ctg="";
                      }
                      elsif(!exists $existed_ctg{$left_extend_ctg} && $left_extend_ctg ne ""){
                          unshift(@extend_line,$left_extend_ctg);
                          $all_existed_ctgs{$left_extend_ctg}=0;
                          $all_existed_pairs{$left_extend_ctg}{$left_extend}=0;
                          $existed_pairs{$left_extend_ctg}{$left_extend}=0;
                          $existed_pairs{$left_extend_ctg}{$left_extend}=0;
                      }
                      $existed_ctg{$left_extend_ctg}=0;                
                }
        

                ###################### 开始延伸右边 #############################################################
                my $right_extend_ctg="1";
                my $right_extend_count=0;
                my $right_overlap=0;  
                while($right_extend_ctg ne ""){
                    $right_extend_ctg="";
                    $right_extend_count=0;
                    $right_overlap=0;
                    $right_extend=$extend_line[-1];
#                    last if(exists $extended_ctg{$right_extend});          
                    $extended_ctg{$right_extend}=0;
                    foreach my $key1 (keys %ctg_ctg_line){
                        my $temp=$ctg_ctg_line{$key1};
                        foreach my $key2 (keys %$temp){
                           if($key2 eq $right_extend){
                               if($ctg_ctg_ori{$key1}{$key2}{2} ne $ctg_ctg_ori{$right_extend}{$extend_line[-2]}{1}){
                                    if($ctg_ctg_count{$key1}{$key2}>$right_extend_count){
                                               $right_extend_ctg=$key1;
                                               $right_extend_count=$ctg_ctg_count{$key1}{$key2};
                                               $right_overlap=$ctg_ctg_overlap{$key1}{$key2};
                                     }
                                     elsif($ctg_ctg_count{$key1}{$key2}==$right_extend_count){
                                               if($ctg_ctg_overlap{$key1}{$key2}>$right_overlap){
                                                  $right_extend_ctg=$key1;
                                                  $right_extend_count=$ctg_ctg_count{$key1}{$key2};
                                                  $right_overlap=$ctg_ctg_overlap{$key1}{$key2};
                                               }
                                      }
                                }
                            }
                         }
                    }
                    my $temp=$ctg_ctg_line{$right_extend};
                    foreach my $key1 (keys %$temp ){
                        next if($key1 eq $extend_line[-2]);
                        if($ctg_ctg_ori{$right_extend}{$key1}{1} ne $ctg_ctg_ori{$right_extend}{$extend_line[-2]}{1}){
                            if($ctg_ctg_count{$right_extend}{$key1}>$right_extend_count){
                                 $right_extend_ctg=$key1;
                                 $right_extend_count=$ctg_ctg_count{$right_extend}{$key1};
                                 $right_overlap=$ctg_ctg_overlap{$right_extend}{$key1};
                            }
                            elsif($ctg_ctg_count{$right_extend}{$key1}==$right_extend_count){
                                 if($ctg_ctg_overlap{$right_extend}{$key1}>$right_overlap){
                                     $right_extend_ctg=$key1;
                                     $right_extend_count=$ctg_ctg_count{$right_extend}{$key1};
                                     $right_overlap=$ctg_ctg_overlap{$right_extend}{$key1};
                                 }
                            }
                         }
                     }
                     if(exists $all_existed_ctgs{$right_extend_ctg}){
                          last;
                      }
                     if(exists $existed_ctg{$right_extend_ctg}){
                         $right_extend_ctg="";
                     }
                     elsif(!exists $existed_ctg{$right_extend_ctg} && $right_extend_ctg ne ""){
                         push(@extend_line,$right_extend_ctg);
                         $all_existed_ctgs{$right_extend_ctg}=0;
                         $all_existed_pairs{$right_extend}{$right_extend_ctg}=0;
                         $existed_pairs{$right_extend}{$right_extend_ctg}=0;
                         $existed_pairs{$right_extend_ctg}{$right_extend}=0;
                     }
                     $existed_ctg{$right_extend_ctg}=0;
                }
                my $ctg_cluster=$extend_line[0];

                for(my $i=1;$i<@extend_line;$i++){
                      $ctg_cluster=$ctg_cluster."-".$extend_line[$i];
                }
#                print "$left_ctg $right_ctg----$ctg_cluster\n";
                my $sign=0;
                for(my $i=0;$i<@cluster_line;$i++){
                      if($cluster_line[$i]=~/$ctg_cluster/){
                             $sign=1;
                      }
                      elsif($ctg_cluster=~/$cluster_line[$i]/){
                             $cluster_line[$i]=$ctg_cluster;
                             $sign=1;
                      }
                      else{
                             my @content=split "-",$ctg_cluster;
                             my $reverse=join ("-",reverse(@content));
                             if($reverse=~/$cluster_line[$i]/){
                                 $cluster_line[$i]=$reverse;
                                 $sign=1;
                             }
                             elsif($cluster_line[$i]=~/$reverse/){
                                 $sign=1;
                             }
                      }
                 }
                 ################# 检查所有已经找到的通路中是否有包含当前通路中所有contig的通路 ##############
 
                 for(my $j=0;$j<@cluster_line;$j++){
                      $cluster_line[$j]=~s/^\s+|\s+$//g;
                      my @target_ctg=split "-",$cluster_line[$j];
                      my %ctg_target=();
                      my %ctg_query=();
                      for(my $target=0;$target<@target_ctg;$target++){
                           $ctg_target{$target_ctg[$target]}=0;
                      }
                      for(my $query=0;$query<@extend_line;$query++){
                           next if($extend_line[$query] eq "");
                           $ctg_query{$extend_line[$query]}=0;
                      }
                      my $target_sign=0;
                      foreach my $target_key(keys %ctg_target){
                           if(!exists $ctg_query{$target_key}){
                                    $target_sign=1;
                           }
                      }
                      if($target_sign==0){
                           $cluster_line[$j]=$ctg_cluster;
                           $sign=1;
                      }
                      my $query_sign=0;
                      foreach my $query_key(keys %ctg_query){
                           if(!exists $ctg_target{$query_key}){
                                    $query_sign=1;
                           }
                      }
                      if($query_sign==0){
                           $sign=1;
                      }
#                      print "1----$ctg_cluster\n2----$cluster_line[$j]\n3----$target_sign\t$query_sign\n";
                 }
                 ############################ end #######################################################  
                
                 if($sign==0){
                    if(!exists $All_existed_clusters{$ctg_cluster}){ 
                         $ctg_cluster=~s/^\s+|\s+$//g;
                         push(@cluster_line,$ctg_cluster);
                         $All_existed_clusters{$ctg_cluster}=0;
                    }
                 }           
          }
    }
#    print "@cluster_line\n";
    my %hash=();
    my @Final_cluster=();
    for(my $i=0;$i<@cluster_line;$i++){
        my @content=split "-",$cluster_line[$i];
        next if(exists $hash{$cluster_line[$i]});
        $hash{$cluster_line[$i]}=0;
        push (@Final_cluster,$cluster_line[$i]);
        for(my $j=0;$j<@content;$j++){
            print OUT "$content[$j] ";
        } 
        print OUT "\n";
    }
    close(OUT);
    return (@Final_cluster);
}                                          
sub draw_ctg_clusters{
    my $info=shift;
    my @ctg_cluster=@{$info};
#    my $ClusterGraph   = Graph::Easy->new( undirected => 1 );
    my %existed_ctg=();
#    print "@ctg_cluster\n";
    open CLUSTER,">cluster_ori.txt" or die $!;
    for(my $i=0;$i<@ctg_cluster;$i++){
        my @content=split "-",$ctg_cluster[$i];
=pod
        for(my $j=0;$j<@content;$j++){
             next if(exists $existed_ctg{$content[$j]});
             $EasyGraph   -> add_node( $content[$j] );
        }
=cut
        for(my $j=0;$j<@content-1;$j++){
             my @infor=split "-",$ctg_ctg_line{$content[$j]}{$content[$j+1]};
             print CLUSTER "$content[$j] $ctg_ctg_chain{$content[$j]}{$content[$j+1]}{1} $ctg_ctg_ori{$content[$j]}{$content[$j+1]}{1} $ctg_ctg_chain{$content[$j]}{$content[$j+1]}{2} $ctg_ctg_ori{$content[$j]}{$content[$j+1]}{2} ";
        }
        print CLUSTER "$content[-1]\n";
     }
#     &drawEasyGraph ( $ClusterGraph , "ctg_cluster" , '01' );
     
}     
                         
=pod 
open OUT1,">pairs.txt" or die $!;
foreach my $key1 (keys %NodePairsHash){
       my $temp=$NodePairsHash{$key1};
       my @content=split "-",$key1;
       foreach my $pacbio (keys %$temp){
            print OUT1 "$content[0] $NodePairsHash{$key1}{$pacbio}{'First'}{'RefPos'} $NodePairsHash{$key1}{$pacbio}{'First'}{'RefStrand'} $NodePairsHash{$key1}{$pacbio}{'First'}{'Identity'} $NodePairsHash{$key1}{$pacbio}{'First'}{'RefStart'} $NodePairsHash{$key1}{$pacbio}{'First'}{'RefEnd'} $NodePairsHash{$key1}{$pacbio}{'First'}{'RefLength'}-----$content[1] $NodePairsHash{$key1}{$pacbio}{'Next'}{'RefPos'} $NodePairsHash{$key1}{$pacbio}{'Next'}{'RefStrand'} $NodePairsHash{$key1}{$pacbio}{'Next'}{'Identity'} $NodePairsHash{$key1}{$pacbio}{'Next'}{'RefStart'} $NodePairsHash{$key1}{$pacbio}{'Next'}{'RefEnd'} $NodePairsHash{$key1}{$pacbio}{'Next'}{'RefLength'}\n";
       }
}
=cut

sub Finding_All_pathway{
    my $ctg_pair_count = shift;
    my $ctg_ctg_overlap = shift;
    my $ctg_ori = shift;
    my $ctg_chain = shift;
    my $ctg_line =shift ;
    my $SEPoint=shift;
#    my $PairsHash = shift;
#   $Pair_Graph
    my %ctg_ctg_count = %{$ctg_pair_count};
    my %ctg_ctg_overlap = %{$ctg_ctg_overlap};
    my %ctg_ctg_ori = %{$ctg_ori};
    my %ctg_ctg_chain = %{$ctg_chain};
    my %ctg_ctg_line = %{$ctg_line};
    my %Scaffold_Contig_SEPoint = %{ $SEPoint};
#    my %NodePairsHash = %{$PairsHash};
    my %All_existed_clusters=();         #$All_existed_clusters
    my @cluster_line=();                 #ctg clusters  
    open OUT,">ctg_clusters.txt" or die $!;
    print "Finding_All_Pathway\n";
    my @Final_Cluster=();
    #首先判断一对contig中，是否有一个contig是最好的，如果两个contig都不是最好的，过；否则延伸

    my %existed_pairs=();
    my %all_existed_ctgs=();
    my %all_existed_pairs=();

    #以Scaffold为单位，两两Contig之间寻找通路
    foreach my $Scaffold (sort {$a cmp $b} keys %Scaffold_Contig_SEPoint){
       my $temp_scaffold=$Scaffold_Contig_SEPoint{$Scaffold};
	  foreach my $scaffold_contig (sort {$a<=>$b} keys %$temp_scaffold){
              my %All_Path=();
	      my $First_Scaffold_Contig=$Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'};    
              #比对到前面Ctg末尾的Pacbio
	      my $next_contig=$scaffold_contig+1;
	      last if(!exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig});
	      my $Second_Scaffold_Contig=$Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'};
              #比对到后面Ctg开头的Pacbio

              #利用Graph：Pair_Graph存储的信息判断两个Contig之间是否存在通路
	      foreach my $First_Pacbio (keys %$First_Scaffold_Contig){
#	          my %exists_Pacbio=();
	          foreach my $Second_Pacbio (keys %$Second_Scaffold_Contig){
#		       next if(exists $exists_Pacbio{$Second_Pacbio});
#		       if($Pair_Graph->is_reachable($First_Pacbio,$Second_Pacbio)){
#			  print "$First_Pacbio $Second_Pacbio reachable\n";

		          if($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==0){      #链为正，right延伸
			     my @neighbors=();
#			     $Pair_Graph->neighbours($First_Pacbio);
			     last if(!@neighbors);
			     #存储与起始pacbio，延伸方向一致的pacbio对
			     for(my $i=0;$i<@neighbors;$i++){
			        next if($First_Pacbio eq $neighbors[$i]);
			        if($ctg_ctg_ori{$First_Pacbio}{$neighbors[$i]}{1} eq "right"){
				   my $current_path=$First_Pacbio."-".$neighbors[$i];
				   $All_Path{$current_path}{$First_Pacbio}=1;
			           $All_Path{$current_path}{$neighbors[$i]}=2;
				   $All_Path{$current_path}{'Total'}=2;
				   $All_Path{$current_path}{'Extend'}=1;          #1表示可延伸，0表示该通路不可延伸
			        }
			      }
				
			      #开始向右延伸pacbio;
			      #如果在同一条通路中遇到相同的节点，则该通路延伸停止；否则一直延伸直到无可延伸的pacbio,或者延伸的pacbio在下一个Contig的开头
			      my $If_Extend=1;
			      while($If_Extend != 0){
			          $If_Extend=0;
			          foreach my $path (keys %All_Path){
				      my @Pacbio_In_Path=split /-/,$path;
				      my $temp_path=$ctg_ctg_ori{$Pacbio_In_Path[-1]};
#				      next if($All_Path{$path}{'Extend'}==0);
				      foreach my $extend_pacbio (keys %$temp_path){
					  next if($extend_pacbio eq $Pacbio_In_Path[-2]);
				          if($ctg_ctg_ori{$Pacbio_In_Path[-2]}{$Pacbio_In_Path[-1]}{2} ne $ctg_ctg_ori{$Pacbio_In_Path[-1]}{$extend_pacbio}{1}){
					      next if(exists $All_Path{$path}{$extend_pacbio});
					      #如果该通路中已经包含该pacbio，则停止对该通路的延伸
					      my $path_extended=$path."-".$extend_pacbio;
					      $All_Path{$path_extended}{'Total'}=$All_Path{$path}{'Total'}+1;
					      $All_Path{$path_extended}{$extend_pacbio}=$All_Path{$path_extended}{'Total'};
#					      if(exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'}{$extend_pacbio}){
#					           $All_Path{$path_extended}{'Extend'}=0;
#					      }
					      $If_Extend=1;
				   	   }
				       }
				       if($If_Extend==1){
				           delete ($All_Path{$path});
				       }
				       elsif($If_Extend==0){
#				           $All_Path{$path}{'Extend'}=0;
				       }
			          }
			      }
			      foreach my $Final_Path (keys %All_Path){
                                  push(@Final_Cluster,$Final_Path);
                                  my @content=split /-/,$Final_Path;
                                  print OUT "@content\n";
                              }
                              %All_Path=();
		          }
                         ##对于正链比对到前一个Ctg上pacbio，向右延伸结束
                         ##对于反链比对到前一个Ctg末尾的pacbio，开始向该pacbio左端延伸
		         else{       #链为负，left延伸
			      my @neighbors=();
#			      $Pair_Graph->neighbours($First_Pacbio);
			      last if(!@neighbors);
			      for(my $i=0;$i<@neighbors;$i++){
                                next if($First_Pacbio eq $neighbors[$i]);
                                if($ctg_ctg_ori{$First_Pacbio}{$neighbors[$i]}{1} eq "left"){
                                   my $current_path=$First_Pacbio."-".$neighbors[$i];
                                   $All_Path{$current_path}{$First_Pacbio}=1;
                                   $All_Path{$current_path}{$neighbors[$i]}=2;
                                   $All_Path{$current_path}{'Total'}=2;
				   $All_Path{$current_path}{'Extend'}=1;
				}
			      }
			      my $If_Extend=1;
                              while($If_Extend != 0){
                                  $If_Extend=0;
                                  foreach my $path (keys %All_Path){
                                      my @Pacbio_In_Path=split /-/,$path;
                                      my $temp_path=$ctg_ctg_ori{$Pacbio_In_Path[-1]};
                                      next if($All_Path{$path}{'Extend'}==0);
                                      foreach my $extend_pacbio (keys %$temp_path){
					  next if($extend_pacbio eq $Pacbio_In_Path[-2]);
                                          if($ctg_ctg_ori{$Pacbio_In_Path[-2]}{$Pacbio_In_Path[-1]}{2} ne $ctg_ctg_ori{$Pacbio_In_Path[-1]}{$extend_pacbio}{1}){
                                              next if(exists $All_Path{$path}{$extend_pacbio});
			                      my $path_extended=$path."-".$extend_pacbio;
                                              $All_Path{$path_extended}{'Total'}=$All_Path{$path}{'Total'}+1;
                                              $All_Path{$path_extended}{$extend_pacbio}=$All_Path{$path_extended}{'Total'};
                                              if(exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'}{$extend_pacbio}){
#                                                   $All_Path{$path_extended}{'Extend'}=0;
                                              }
                                              $If_Extend=1;
                                           }
                                       }
                                       if($If_Extend==1){
                                           delete ($All_Path{$path});
                                       }
                                       elsif($If_Extend==0){
#                                           $All_Path{$path}{'Extend'}=0;
                                       }
                                  }
                              }				
			 }
			 foreach my $Final_Path (keys %All_Path){
			      push(@Final_Cluster,$Final_Path);
			      my @content=split /-/,$Final_Path;
			      print OUT "@content\n";
			 }
			 %All_Path=();
#                       }                   
                  }                        		
              }
          }                        		
    }                       		         
    close(OUT);
    return (@Final_Cluster);
}                                          
sub Finding_Best_Pathway{
    my $ctg_pair_count = shift;
    my $ctg_ctg_overlap = shift;
    my $ctg_ori = shift;
    my $ctg_chain = shift;
    my $ctg_line =shift ;
    my $SEPoint=shift;
    my %ctg_ctg_count = %{$ctg_pair_count};
    my %ctg_ctg_overlap = %{$ctg_ctg_overlap};
    my %ctg_ctg_ori = %{$ctg_ori};
    my %ctg_ctg_chain = %{$ctg_chain};
    my %ctg_ctg_line = %{$ctg_line};
    my %Scaffold_Contig_SEPoint = %{ $SEPoint};
    my %All_existed_clusters=();         #$All_existed_clusters
    my @cluster_line=();                 #ctg clusters  
    open OUT,">ctg_clusters.txt" or die $!;
    open CLUSTER,">cluster_ori.txt" or die $!;
    #首先判断一对contig中，是否有一个contig是最好的，如果两个contig都不是最好的，过；否则延伸

    my %existed_pairs=();
    my %all_existed_ctgs=();
    my %all_existed_pairs=();
    my $Bestn=1;
    print OUT "Real Starting Extending\n";
    my $sign=0;
################################ test #############################################
    open TEST,">test.txt" or die $!;
=pod
    foreach my $key (keys %ctg_ctg_count){
	my $temp=$ctg_ctg_count{$key};
	my @num=keys %$temp;
	my $total=@num;
	print TEST "$key $total\n";
    }
    close(TEST);
=cut
############################### test ################################################

    foreach my $Scaffold (sort {$a cmp $b} keys %Scaffold_Contig_SEPoint){
#	  next if($Scaffold ne "Super-Scaffold_119");
          my $time=&gettime("yyyy-mm-dd hh:mi:ss");
          print TEST "$time\n";
	  print TEST "Begining Super-Scaffold_119\n";
	  close(TEST);
	  open TEST,">>test.txt" or die $!;
          my $temp_scaffold=$Scaffold_Contig_SEPoint{$Scaffold};
          foreach my $scaffold_contig (sort {$a<=>$b} keys %$temp_scaffold){
              my $First_Scaffold_Contig=$Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'};
	      my $next_contig=$scaffold_contig+1;
	      last if(!exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig});
	      $sign++;
	      my $First_Contig=$Scaffold.".".$scaffold_contig;
	      my $Second_Contig=$Scaffold.".".$next_contig;
	      my $Real_Gap=$Contig_Gap_Length{$First_Contig}{$Second_Contig};
	      my $Second_Scaffold_Contig=$Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'};
	      my %Bestn_PathWay=();
	      print TEST "Begining $scaffold_contig\n";
              close(TEST);
              open TEST,">>test.txt" or die $!;
	      foreach my $First_Pacbio (keys %$First_Scaffold_Contig){
		    
                     my $Start_Ori="";
                     if($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==0){
                        $Start_Ori="right";
                     }
                     elsif($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==1){
                        $Start_Ori="left";
                     }
                     my $First_Temp=$Reads_Pairs{$First_Pacbio}{$Start_Ori};

                     my $Current_Best=0;
		     my %All_Path=();
                     foreach my $Query (sort {$Reads_Pairs{$First_Pacbio}{$Start_Ori}{$b}{'S'} <=> $Reads_Pairs{$First_Pacbio}{$Start_Ori}{$a}{'S'}} keys %$First_Temp){
			 last if($Current_Best>=$Bestn);
			 next if($First_Pacbio eq $Query);
                         if(($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==0 && $Start_Ori eq "right") || (($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==1 && $Start_Ori eq "left"))){
			     my $current_path=$First_Pacbio."-".$Query;
			     if(exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'}{$Query}){
				     $All_Path{$current_path}=1;
			     }
			     else{
		                     $All_Path{$current_path}=0;
			     }
			     $Current_Best++;
			 }
		     }
		     my $If_Extend=1;
                     print TEST "$Scaffold\t$scaffold_contig\t$First_Pacbio\n";
                     close(TEST);
                     open TEST,">>test.txt" or die $!;
		     while($If_Extend==1){
			 $If_Extend=0;
			 my @Total_Keys=keys %All_Path;
			 my $total=@Total_Keys;
			 print TEST "$total\n";
			 close(TEST);
                         open TEST,">>test.txt" or die $!;
#			 next if(@Total_Keys>1000);
			 foreach my $Current_Path (keys %All_Path){
			     next if($All_Path{$Current_Path}==1);
			     my @extend_line=split /-/,$Current_Path;
			     my $Extend_Rate=@extend_line*1000/($Real_Gap+15000);
=pod
			     if(@extend_line>100){
				     delete($All_Path{$Current_Path});
				     next;
			     }
=cut
                             my $right_extend=$extend_line[-1];
                             ######### 开始延伸右边 ###################################
                              my $Current_Start_Ori="";
                              if(exists $Reads_Pairs{$right_extend}{'left'}{$extend_line[-2]}){
                                    $Current_Start_Ori="right";
                              }
                              elsif(exists $Reads_Pairs{$right_extend}{'right'}{$extend_line[-2]}){
                                    $Current_Start_Ori="left";
                              }

                             my $temp=$Reads_Pairs{$right_extend}{$Current_Start_Ori};
			     $Current_Best=0;
                             foreach my $key1 (sort {$Reads_Pairs{$right_extend}{$Current_Start_Ori}{$b}{'S'} <=> $Reads_Pairs{$right_extend}{$Current_Start_Ori}{$a}{'S'}}keys %$temp ){
				     last if($Current_Best>=$Bestn);
                                     next if($key1 eq $extend_line[-2]);
                                     next if($Current_Path=~/$key1/);
                                     if(exists $Reads_Pairs{$right_extend}{$Current_Start_Ori}{$key1}){
                                         my $extended_path=$Current_Path."-".$key1;
					 if(exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'}{$key1}){
						if(!exists $All_Path{$Current_Path}){
						    $All_Path{$extended_path}=1;
					            $Current_Best++;
						 }
						 else{
						    $All_Path{$extended_path}=1;
						    delete($All_Path{$Current_Path});
						    $Current_Best++;
						 }
					 }
					 else{
						 if(!exists $All_Path{$Current_Path}){
						    $All_Path{$extended_path}=0;
						    $Current_Best++;
						    $If_Extend=1;
						 }
						 else{
						    $All_Path{$extended_path}=0;
						    $Current_Best++;
						    $If_Extend=1;
						    delete($All_Path{$Current_Path});
						 }
					 }
				     }
                             }
                         }
                     }
#		     print OUT "End of Extending\n";
#		     my $time=&gettime("yyyy-mm-dd hh:mi:ss");
#                     print OUT "$time\n";
		     foreach my $best (keys %All_Path){
                         if($All_Path{$best}==1){
                             my @content=split /-/,$best;
                             print OUT "@content\t$Scaffold-$scaffold_contig-$next_contig-Passing\n";
			     for(my $j=0;$j<@content-1;$j++){
				 my $First_Chain=0;
                                 my $Second_Chain=0;
                                 my $First_Ori="";
                                 my $Second_Ori="";
                                 if(exists $Reads_Pairs{$content[$j]}{'left'}{$content[$j+1]}{'C'}){
                                        $Second_Chain=$Reads_Pairs{$content[$j]}{'left'}{$content[$j+1]}{'C'};
                                        $Second_Ori=$Reads_Pairs{$content[$j]}{'left'}{$content[$j+1]}{'O'};
                                 }
                                 elsif(exists $Reads_Pairs{$content[$j]}{'right'}{$content[$j+1]}{'C'}){
                                        $Second_Chain=$Reads_Pairs{$content[$j]}{'right'}{$content[$j+1]}{'C'};
                                        $Second_Ori=$Reads_Pairs{$content[$j]}{'right'}{$content[$j+1]}{'O'};
                                 }

                                 if(exists $Reads_Pairs{$content[$j+1]}{'left'}{$content[$j]}{'C'}){
                                        $First_Chain=$Reads_Pairs{$content[$j+1]}{'left'}{$content[$j]}{'C'};
                                        $First_Ori=$Reads_Pairs{$content[$j+1]}{'left'}{$content[$j]}{'O'};
                                 }
                                 elsif(exists $Reads_Pairs{$content[$j+1]}{'right'}{$content[$j]}{'C'}){
                                        $First_Chain=$Reads_Pairs{$content[$j+1]}{'right'}{$content[$j]}{'C'};
                                        $First_Ori=$Reads_Pairs{$content[$j+1]}{'right'}{$content[$j]}{'O'};
                                 }

                                 print CLUSTER "$content[$j] $First_Chain $First_Ori $Second_Chain $Second_Ori ";
			     }
			     print CLUSTER "$content[-1]\n";
                         }
                     }
                     %All_Path=();
              }
          }
    }
    print "$sign\n";
    close(OUT);
    close(CLUSTER);
}
sub gettime {
    $_ = shift;
    my $t = shift;
    (!$t) and ($t = time);
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($t);
    $year += 1900;
    my $yy = substr $year,2,2;
    $mon++;
    s/yyyy/$year/gi;
    s/yy/$yy/gi;
    if ($mon < 10)  { s/mm/0$mon/gi;  } else { s/mm/$mon/gi; }
    if ($mday < 10) { s/dd/0$mday/gi; } else { s/dd/$mday/gi; }
    if ($hour < 10) { s/hh/0$hour/gi; } else { s/hh/$hour/gi; }
    if ($min < 10)  { s/mi/0$min/gi;  } else { s/mi/$min/gi; }
    if ($sec < 10)  { s/ss/0$sec/gi;  } else { s/ss/$sec/gi; }

    $_;
}

sub Finding_MaxExtending_Pathway{
    my $ctg_pair_count = shift;
    my $ctg_ctg_overlap = shift;
    my $ctg_ori = shift;
    my $ctg_chain = shift;
    my $ctg_line =shift ;
    my $SEPoint=shift;
    my $ctg_extend=shift;
    my $ctg_info=shift;
                    
    my %ctg_ctg_count = %{$ctg_pair_count};
    my %ctg_ctg_overlap = %{$ctg_ctg_overlap};
    my %ctg_ctg_ori = %{$ctg_ori};
    my %ctg_ctg_chain = %{$ctg_chain};
    my %ctg_ctg_line = %{$ctg_line};
    my %Scaffold_Contig_SEPoint = %{ $SEPoint}; 
    my %ctg_ctg_extend=%{$ctg_extend};
    my %ctg_ctg_info=%{$ctg_info};
    my %All_existed_clusters=();         #$All_existed_clusters
    my @cluster_line=();                 #ctg clusters  
    open OUT,">>ctg_clusters.txt" or die $!;
    open CLUSTER,">>cluster_ori.txt" or die $!;
    #首先判断一对contig中，是否有一个contig是最好的，如果两个contig都不是最好的，过；否则延伸

    my %existed_pairs=();
    my %all_existed_ctgs=();
    my %all_existed_pairs=();
    my $Bestn=1;
    my $sign=0;
#    my $MinIdentity=98;
#    my $MinOverlap=5000;
#    my $MinIdentity=97;
#    my $MinOverlap=2000;
#    my $MinIdentity=96;
#    my $MinOverlap=1000;
#    my $MinIdentity=97;
#    my $MinOverlap=500;
    
    for(my $MinIdentity=94;$MinIdentity<100;$MinIdentity++){
    for(my $MinOverlap=500;$MinOverlap<2000;$MinOverlap=$MinOverlap+500){
    foreach my $Scaffold (sort {$a cmp $b} keys %Scaffold_Contig_SEPoint){
          my $temp_scaffold=$Scaffold_Contig_SEPoint{$Scaffold};
          foreach my $scaffold_contig (sort {$a<=>$b} keys %$temp_scaffold){
              my $First_Scaffold_Contig=$Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'};
	      my $next_contig=$scaffold_contig+1;
	      last if(!exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig});
	      $sign++;
	      my $First_Contig=$Scaffold.".".$scaffold_contig;
	      my $Second_Contig=$Scaffold.".".$next_contig;
	      my $Real_Gap=$Contig_Gap_Length{$First_Contig}{$Second_Contig};
	      my $Second_Scaffold_Contig=$Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'};
	      my %Bestn_PathWay=();

	      foreach my $First_Pacbio (keys %$First_Scaffold_Contig){
		     my $Start_Ori="";
                     if($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==0){
                        $Start_Ori="right";
                     }
                     elsif($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==1){
                        $Start_Ori="left";
                     }
                     my $First_Temp=$Reads_Pairs{$First_Pacbio}{$Start_Ori};

                     my $Current_Best=0;
		     my %All_Path=();
                     foreach my $Query (sort {$Reads_Pairs{$First_Pacbio}{$Start_Ori}{$b}{'E'} <=> $Reads_Pairs{$First_Pacbio}{$Start_Ori}{$a}{'E'}}keys %$First_Temp){
			 last if($Current_Best>=$Bestn);
			 next if($First_Pacbio eq $Query);
	 		 next if($Reads_Pairs{$First_Pacbio}{$Start_Ori}{$Query}{'I'}<$MinIdentity || $Reads_Pairs{$First_Pacbio}{$Start_Ori}{$Query}{'S'}<$MinOverlap);
                         if(($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==0 && $Start_Ori eq "right") || (($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==1 && $Start_Ori eq "left"))){
			     my $current_path=$First_Pacbio."-".$Query;
			     if(exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'}{$Query}){
				     $All_Path{$current_path}=1;
			     }
			     else{
		                     $All_Path{$current_path}=0;
			     }
			     $Current_Best++;
			 }
		     }
		     my $If_Extend=1;
		     while($If_Extend==1){
			 $If_Extend=0;
			 my @Total_Keys=keys %All_Path;
			 my $total=@Total_Keys;
#			 next if(@Total_Keys>1000);
			 foreach my $Current_Path (keys %All_Path){
			     next if($All_Path{$Current_Path}==1);
			     my @extend_line=split /-/,$Current_Path;
			     my $Extend_Rate=@extend_line*1000/($Real_Gap+15000);
=pod
			     if(@extend_line>50){
				     next;
			     }
=cut
                             my $right_extend=$extend_line[-1];
                             ######### 开始延伸右边 ###################################
                             my $Current_Start_Ori="";
                             if(exists $Reads_Pairs{$right_extend}{'left'}{$extend_line[-2]}){
                                     $Current_Start_Ori="right";
                             }
                             elsif(exists $Reads_Pairs{$right_extend}{'right'}{$extend_line[-2]}){
                                     $Current_Start_Ori="left";
                             }

                             my $temp=$Reads_Pairs{$right_extend}{$Current_Start_Ori};
			     $Current_Best=0;
                             foreach my $key1 (sort {$Reads_Pairs{$right_extend}{$Current_Start_Ori}{$b}{'E'} <=> $Reads_Pairs{$right_extend}{$Current_Start_Ori}{$a}{'E'}}keys %$temp ){
				     last if($Current_Best>=$Bestn);
                                     next if($key1 eq $extend_line[-2]);
                                     next if($Current_Path=~/$key1/);
                                     next if($Reads_Pairs{$right_extend}{$Current_Start_Ori}{$key1}{'I'}<$MinIdentity || $Reads_Pairs{$right_extend}{$Current_Start_Ori}{$key1}{'S'}<$MinOverlap);
                                     if(exists $Reads_Pairs{$right_extend}{$Current_Start_Ori}{$key1}){
                                         my $extended_path=$Current_Path."-".$key1;
					 if(exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'}{$key1}){
						if(!exists $All_Path{$Current_Path}){
						    $All_Path{$extended_path}=1;
					            $Current_Best++;
						 }
						 else{
						    $All_Path{$extended_path}=1;
						    delete($All_Path{$Current_Path});
						    $Current_Best++;
						 }
					 }
					 else{
						 if(!exists $All_Path{$Current_Path}){
						    $All_Path{$extended_path}=0;
						    $Current_Best++;
						    $If_Extend=1;
						 }
						 else{
						    $All_Path{$extended_path}=0;
						    $Current_Best++;
						    $If_Extend=1;
						    delete($All_Path{$Current_Path});
						 }
					 }
				     }
                             }
                         }
                     }
		     foreach my $best (keys %All_Path){
                         if($All_Path{$best}==1){
                             my @content=split /-/,$best;
                             print OUT "@content\t$Scaffold-$scaffold_contig-$next_contig-Passing\n";
			     for(my $j=0;$j<@content-1;$j++){
				 my $First_Chain=0;
                                 my $Second_Chain=0;
                                 my $First_Ori="";
                                 my $Second_Ori="";
                                 if(exists $Reads_Pairs{$content[$j]}{'left'}{$content[$j+1]}{'C'}){
                                        $Second_Chain=$Reads_Pairs{$content[$j]}{'left'}{$content[$j+1]}{'C'};
                                        $Second_Ori=$Reads_Pairs{$content[$j]}{'left'}{$content[$j+1]}{'O'};
                                 }
                                 elsif(exists $Reads_Pairs{$content[$j]}{'right'}{$content[$j+1]}{'C'}){
                                        $Second_Chain=$Reads_Pairs{$content[$j]}{'right'}{$content[$j+1]}{'C'};
                                        $Second_Ori=$Reads_Pairs{$content[$j]}{'right'}{$content[$j+1]}{'O'};
                                 }

                                 if(exists $Reads_Pairs{$content[$j+1]}{'left'}{$content[$j]}{'C'}){
                                        $First_Chain=$Reads_Pairs{$content[$j+1]}{'left'}{$content[$j]}{'C'};
                                        $First_Ori=$Reads_Pairs{$content[$j+1]}{'left'}{$content[$j]}{'O'};
                                 }
                                 elsif(exists $Reads_Pairs{$content[$j+1]}{'right'}{$content[$j]}{'C'}){
                                        $First_Chain=$Reads_Pairs{$content[$j+1]}{'right'}{$content[$j]}{'C'};
                                        $First_Ori=$Reads_Pairs{$content[$j+1]}{'right'}{$content[$j]}{'O'};
                                 }

                                 print CLUSTER "$content[$j] $First_Chain $First_Ori $Second_Chain $Second_Ori ";

			     }
			     print CLUSTER "$content[-1]\n";
                         }
                     }
                     %All_Path=();
              }
          }
    }
    }
    }
    print "$sign\n";
    close(OUT);
    close(CLUSTER);
}

sub Finding_Randoming_Pathway{
    my $ctg_pair_count = shift;
    my $ctg_ctg_overlap = shift;
    my $ctg_ori = shift;
    my $ctg_chain = shift;
    my $ctg_line =shift ;
    my $SEPoint=shift;
    my $ctg_extend=shift;
    my $ctg_info=shift;
                    
    my %ctg_ctg_count = %{$ctg_pair_count};
    my %ctg_ctg_overlap = %{$ctg_ctg_overlap};
    my %ctg_ctg_ori = %{$ctg_ori};
    my %ctg_ctg_chain = %{$ctg_chain};
    my %ctg_ctg_line = %{$ctg_line};
    my %Scaffold_Contig_SEPoint = %{ $SEPoint}; 
    my %ctg_ctg_extend=%{$ctg_extend};
    my %ctg_ctg_info=%{$ctg_info};
    my %All_existed_clusters=();         #$All_existed_clusters
    my @cluster_line=();                 #ctg clusters  
    open OUT,">>ctg_clusters.txt" or die $!;
    open CLUSTER,">>cluster_ori.txt" or die $!;
    #首先判断一对contig中，是否有一个contig是最好的，如果两个contig都不是最好的，过；否则延伸

    my %existed_pairs=();
    my %all_existed_ctgs=();
    my %all_existed_pairs=();
    my $Bestn=1;
    my $sign=0;
#    my $MinIdentity=96;
#    my $MinOverlap=1000;
    
    for(my $MinIdentity=94;$MinIdentity<=99;$MinIdentity++){
    for(my $MinOverlap=500;$MinOverlap<2000;$MinOverlap=$MinOverlap+500){
    foreach my $Scaffold (sort {$a cmp $b} keys %Scaffold_Contig_SEPoint){
          my $temp_scaffold=$Scaffold_Contig_SEPoint{$Scaffold};
          foreach my $scaffold_contig (sort {$a<=>$b} keys %$temp_scaffold){
              my $First_Scaffold_Contig=$Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'};
	      my $next_contig=$scaffold_contig+1;
	      last if(!exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig});
	      my $First_Contig=$Scaffold.".".$scaffold_contig;
	      my $Second_Contig=$Scaffold.".".$next_contig;
	      my $Real_Gap=$Contig_Gap_Length{$First_Contig}{$Second_Contig};
	      my $Second_Scaffold_Contig=$Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'};
	      my %Bestn_PathWay=();
	      foreach my $First_Pacbio (keys %$First_Scaffold_Contig){
                     
		     my $Start_Ori="";
                     if($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==0){
                         $Start_Ori="right";
                     }
                     elsif($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==1){
                         $Start_Ori="left";
                     }
                     my $First_Temp=$Reads_Pairs{$First_Pacbio}{$Start_Ori};

                     my $Current_Best=0;
		     my %All_Path=();
		     my @Second_Pacbios= sort {$Reads_Pairs{$First_Pacbio}{$Start_Ori}{$b}{'S'} <=> $Reads_Pairs{$First_Pacbio}{$Start_Ori}{$a}{'S'}}keys %$First_Temp;
		     my @Second_Real_Pacbios=();
		     next if(!@Second_Pacbios);
		     for(my $i=0;$i<@Second_Pacbios;$i++){
			 if(($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==0 && $Start_Ori eq "right") || (($Scaffold_Contig_SEPoint{$Scaffold}{$scaffold_contig}{'Tail'}{$First_Pacbio}{0}==1 && $Start_Ori eq "left"))){	
                               next if($Reads_Pairs{$First_Pacbio}{$Start_Ori}{$Second_Pacbios[$i]}{'I'}<$MinIdentity || $Reads_Pairs{$First_Pacbio}{$Start_Ori}{$Second_Pacbios[$i]}{'S'}<$MinOverlap);
			       next if($First_Pacbio eq $Second_Pacbios[$i]);
                               push(@Second_Real_Pacbios,$Second_Pacbios[$i]);
			 }
		     }
                     next if(!@Second_Real_Pacbios);   
		     my $total=@Second_Real_Pacbios;
                     my $pos=0;
		     if($total>20){
			 $pos=int rand(20);
		     }
		     else{
			 $pos=int rand($total);
		     }
                     my $current_path=$First_Pacbio."-".$Second_Real_Pacbios[$pos];
                     if(exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'}{$Second_Real_Pacbios[$pos]}){
                         $All_Path{$current_path}=1;
                     }
		     else{
			 $All_Path{$current_path}=0;
	             }
		     my $If_Extend=1;
		     while($If_Extend==1){
			 $If_Extend=0;
			 foreach my $Current_Path (keys %All_Path){
			     next if($All_Path{$Current_Path}==1);
			     my @extend_line=split /-/,$Current_Path;
			     my $Extend_Rate=@extend_line*1000/($Real_Gap+15000);
                             my $right_extend=$extend_line[-1];
                             ######### 开始延伸右边 ###################################
                             my $Current_Start_Ori="";
                             if(exists $Reads_Pairs{$right_extend}{'left'}{$extend_line[-2]}){
                                     $Current_Start_Ori="right";
                             }
                             elsif(exists $Reads_Pairs{$right_extend}{'right'}{$extend_line[-2]}){
                                     $Current_Start_Ori="left";
                             }

                             my $temp=$Reads_Pairs{$right_extend}{$Current_Start_Ori};

			     $Current_Best=0;
			     @Second_Pacbios= sort {$Reads_Pairs{$right_extend}{$Current_Start_Ori}{$b}{'S'} <=> $Reads_Pairs{$right_extend}{$Current_Start_Ori}{$a}{'S'}}keys %$temp;
                             @Second_Real_Pacbios=();
                             next if(!@Second_Pacbios);
                             for(my $i=0;$i<@Second_Pacbios;$i++){
                                 if(exists $Reads_Pairs{$right_extend}{$Current_Start_Ori}{$Second_Pacbios[$i]}){
                                      next if($Reads_Pairs{$right_extend}{$Current_Start_Ori}{$Second_Pacbios[$i]}{'I'}<$MinIdentity || $Reads_Pairs{$right_extend}{$Current_Start_Ori}{$Second_Pacbios[$i]}{'S'}<$MinOverlap);
                                      next if($First_Pacbio eq $Second_Pacbios[$i]);
				      next if($Current_Path=~/$Second_Pacbios[$i]/);
                                      push(@Second_Real_Pacbios,$Second_Pacbios[$i]);
                                 }
                             }			     
                             next if(!@Second_Real_Pacbios);
                             $total=@Second_Real_Pacbios;
                             $pos=0;
                             if($total>20){
                                  $pos=int rand(20);
                             }
                             else{
                                  $pos=int rand($total);
                             }
                             my $extended_path=$Current_Path."-".$Second_Real_Pacbios[$pos];
			     $If_Extend=1;
                             if(exists $Scaffold_Contig_SEPoint{$Scaffold}{$next_contig}{'Head'}{$Second_Real_Pacbios[$pos]}){
                                  $All_Path{$extended_path}=1;
                             }
                             else{
                                  $All_Path{$extended_path}=0;
                             }
                             delete($All_Path{$Current_Path});


                         }
                     }
		     foreach my $best (keys %All_Path){
                         if($All_Path{$best}==1){
                             my @content=split /-/,$best;
                             print OUT "@content\t$Scaffold-$scaffold_contig-$next_contig-Passing\n";
			     for(my $j=0;$j<@content-1;$j++){
				 my $First_Chain=0;
                                 my $Second_Chain=0;
                                 my $First_Ori="";
                                 my $Second_Ori="";
                                 if(exists $Reads_Pairs{$content[$j]}{'left'}{$content[$j+1]}{'C'}){
                                        $Second_Chain=$Reads_Pairs{$content[$j]}{'left'}{$content[$j+1]}{'C'};
                                        $Second_Ori=$Reads_Pairs{$content[$j]}{'left'}{$content[$j+1]}{'O'};
                                 }
                                 elsif(exists $Reads_Pairs{$content[$j]}{'right'}{$content[$j+1]}{'C'}){
                                        $Second_Chain=$Reads_Pairs{$content[$j]}{'right'}{$content[$j+1]}{'C'};
                                        $Second_Ori=$Reads_Pairs{$content[$j]}{'right'}{$content[$j+1]}{'O'};
                                 }

                                 if(exists $Reads_Pairs{$content[$j+1]}{'left'}{$content[$j]}{'C'}){
                                        $First_Chain=$Reads_Pairs{$content[$j+1]}{'left'}{$content[$j]}{'C'};
                                        $First_Ori=$Reads_Pairs{$content[$j+1]}{'left'}{$content[$j]}{'O'};
                                 }
                                 elsif(exists $Reads_Pairs{$content[$j+1]}{'right'}{$content[$j]}{'C'}){
                                        $First_Chain=$Reads_Pairs{$content[$j+1]}{'right'}{$content[$j]}{'C'};
                                        $First_Ori=$Reads_Pairs{$content[$j+1]}{'right'}{$content[$j]}{'O'};
                                 }

                                 print CLUSTER "$content[$j] $First_Chain $First_Ori $Second_Chain $Second_Ori ";

			     }
			     print CLUSTER "$content[-1]\n";
                         }
                     }
                     %All_Path=();
              }
          }
    }
    }
    }
    print "$sign\n";
    close(OUT);
    close(CLUSTER);
}

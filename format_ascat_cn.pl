#Pavana Anur
#Script for formatting ASCAT segmentation files for timing. 

#!/usr/bin/perl
use strict;

#Initiate variables
my $samplename= @ARGV[0];
my $tumor_file=$samplename . "_ascatSegmentation.txt\n";
my ($chr,$start,$end,$width,$c1,$c2,$ascatId,$tcn)=(0);
my @data;
my @info;
my @cn;
my @endtmp;

#Open files
open(TUMOR_FILE,$tumor_file) or die "Can't open $tumor_file!\n";

my $out= $samplename . "_segmentation.txt";
open(OUT, ">$out" ) or die "Can't open  : $!";
select OUT; $| = 1; select STDOUT;  

#Format
print OUT "chr\tstart\tend\twidth\tstrand\tc1N\tc2N\tascatId\ttCN\n";

#Get chr,position,variant count and total count from tumor  file 
#read first line
my $dummy=<TUMOR_FILE>;

while(<TUMOR_FILE>)
{
        if((substr $_, 0, 10) !~ /#/)
        { 
        chomp($_);
        my  @data = split(/\t/,$_);
                #print "@data\n";
                $chr=$data[1];
                $start=$data[2];
                $end=$data[3];
                $width=$end-$start;
                $c1=$data[24];
                $c2=$data[25];
                $ascatId=$data[0];
                $tcn=$data[26];
        print OUT "$chr\t$start\t$end\t$width\t*\t$c1\t$c2\t$ascatId\t$tcn\n";
                }
        
        
}

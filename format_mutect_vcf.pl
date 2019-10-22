#Pavana Anur
#Script for formatting Mutect VCF files for timing. 

#!/usr/bin/perl
use strict;

my $samplename= @ARGV[0];
my ($chr,$pos,$altcount,$totcount,$af)=(0);
my $vcf =$samplename . ".vcf";
my $vcfpath="/home/exacloud/lustre1/users/anurpa/timing_cedar/Mutations/$vcf";

open(VCF,$vcfpath) or die "Can't open $vcfpath!\n";


my $outname1= $samplename . "_mutations.txt";
my $out1="/home/exacloud/lustre1/users/anurpa/timing_cedar/Mutations/$outname1";
open(OUT1, ">$out1" ) or die "Can't open  : $!";
select OUT1; $| = 1; select STDOUT;  


print OUT1 "chr\tposition\tt_alt_count\tt_depth\tAlleleFrequency\n";
while (my $line2=<VCF>)
{
        if((substr $line2, 0, 2) !~ /#/)
        { 
        chomp($line2);
        my @columns2 = split(/\t/,$line2);
         $chr=$columns2[0];
         $pos=$columns2[1];
        my $info=$columns2[9];
        chomp($info);
        my @info_split=split(/:/,$info);
        my $ad=$info_split[1];
        chomp($ad);
        my @ad_split=split(/,/,$ad);
           $altcount=$ad_split[1];
        my $totcount=$info_split[3];}
        
print "Done formatting mutation count for $samplename.\n";

        my $af=$info_split[4];
        print OUT1 "$chr\t$pos\t$altcount\t$totcount\t$af\n";
        }
        

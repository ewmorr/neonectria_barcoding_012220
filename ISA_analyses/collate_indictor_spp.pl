#/usr/bin/perl
#Eric Morrison
#09022020
#Collate results files from R indicspecies::multipatt by finding species that are significant indcators across multiple groupings (i.e., different grouping factors)

use strict;
use warnings;

sub parse_files_to_hash{
    my @inFiles = @_;
    
    my %isa_spp;
    my %facLevels;
    foreach my $factor (@inFiles){
        open(TEMP, "$factor") || die "Can't open input\n";
        chomp(my @temp = <TEMP>);
        $factor =~ /(.*)\.txt/;
        my $facNam = $1;
        my $header = shift(@temp);
        my @header = split("\t", $header);
#        my @levels = $header[0 .. $#header -3];
#        $facLevels{facName} = [@levels];
        foreach my $asv (@temp){
            my @asv = split("\t", $asv);
            if($asv[-1] eq "NA"){next;}
            if($asv[-1] < 0.05){
                $isa_spp{$asv[0]}{$facNam} = [ @asv[1 .. $#asv -3] ];
            }
        }
    }
    return(\%isa_spp);
#    return(\%facLevels);
}

sub find_multi_indicators{
    my $isa_sppRef = $_[0];
    my %isa_spp = %$isa_sppRef;
#    my %facLevels = %$facLevelsRef;
    
    foreach my $asv (sort {$a cmp $b} keys %isa_spp){
        if( scalar(keys %{ $isa_spp{$asv} }) > 1){
            print $asv, "\t";
            foreach my $cat (sort {$a cmp $b} keys %{ $isa_spp{$asv} }){
                print $cat, "-";
                foreach my $level (@{ $isa_spp{$asv}{$cat} }){
                    print $level;
                }
                print "\t";
            }
            print "\n";
        }
    }
}

#Main
{
    #read ISA files
    my @inFiles = @ARGV;
    my $isa_sppRef = parse_files_to_hash(@inFiles);
    find_multi_indicators($isa_sppRef);
}

#!/usr/bin/perl
#Eric Morrison
#04082020
#Pass a dir of vsearch --usearch_global mapping files that are named by sample name. This script counts the numer of sequences that hit, and returns each hit with counts as a single sample ASV table with db IDs as ASV names. This can subsequently be joined in R using dplyr::full_join to get all hits across samples.

use strict;
use warnings;

sub process_dir{
    my $dir = $_[0];
    opendir(DIR, $dir) || die "Can't open directory\n";
    while(my $file = readdir(DIR)){
        print $file, "\n";
    }
}

#MAIN
{
    my $dir = $ARGV[0];
    process_dir($dir);
}

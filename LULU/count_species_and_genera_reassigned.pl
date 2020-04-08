#!/usr/bin/perl
#Eric Morrison
#04082020

use strict;
use warnings;

sub process_taxa{
    my $file = $_[0];
    open(IN, "$file") || die "Can not open input.\n";
    chomp(my @in = <IN>);
    shift(@in);
    my %taxa;
    foreach my $line (@in){
        my @line = split("\t", $line);
        my @child = split(";", $line[0]);
        my @parent = split(";", $line[1]);
        #$line[-2] is indexing hash by child name (one index per "discarded" ASV)
        #-1 index of child and index is species name
        #-2 is genus
        $taxa{$line[-2]}{"species"} = [$child[-1], $parent[-1], $line[-2], $line[-4], $line[-3]];
        $taxa{$line[-2]}{"genus"} = [$child[-2], $parent[-2], $line[-2], $line[-4], $line[-3]];
    }
    return(\%taxa);
}

sub num_reassigned{
    my $taxaRef = $_[0];
    my %taxa = %$taxaRef;
    my $sppReassigned = 0;
    my $genReassigned = 0;
    foreach my $child (sort{$a cmp $b} keys %taxa){
        if($taxa{$child}{"species"}[0] ne "NA" & $taxa{$child}{"species"}[0] ne $taxa{$child}{"species"}[1]){
            $sppReassigned++;
        }
        if($taxa{$child}{"genus"}[0] ne "NA" & $taxa{$child}{"genus"}[0] ne $taxa{$child}{"genus"}[1]){
            $genReassigned++;
        }
    }
    print "\nspp reassgined = $sppReassigned\ngenera reassigned = $genReassigned\n\n";
}

sub which_reassigned{
    my $taxaRef = $_[0];
    my $out = $_[1];
    open(OUT, ">$out") || die "Can't open output.\n";
    my %taxa = %$taxaRef;
    print OUT "tax_level\tchild_id\tsparent_id\teq_sim\tn_samples_child\tchild_tax\tparent_tax\n";
    foreach my $child (sort{$a cmp $b} keys %taxa){
        if($taxa{$child}{"species"}[0] ne "NA" & $taxa{$child}{"species"}[0] ne $taxa{$child}{"species"}[1]){
                print OUT "species", "\t", $taxa{$child}{"species"}[2], "\t", $taxa{$child}{"species"}[3], "\t", $taxa{$child}{"species"}[4], "\t", $taxa{$child}{"genus"}[0],";",$taxa{$child}{"species"}[0], "\t", $taxa{$child}{"genus"}[1],";",$taxa{$child}{"species"}[1], "\n";
        }elsif($taxa{$child}{"genus"}[0] ne "NA" & $taxa{$child}{"genus"}[0] ne $taxa{$child}{"genus"}[1]){
                print OUT "genus", "\t", $taxa{$child}{"genus"}[2], "\t", $taxa{$child}{"genus"}[3], "\t", $taxa{$child}{"genus"}[4], "\t", $taxa{$child}{"genus"}[0],";",$taxa{$child}{"species"}[0], "\t", $taxa{$child}{"genus"}[1],";",$taxa{$child}{"species"}[1],  "\n";
        }
    }
}

#MAIN
{
    my $file = $ARGV[0];
    my $out = $ARGV[1];
    my $taxaRef = process_taxa($file);
    num_reassigned($taxaRef);
    which_reassigned($taxaRef, $out);
}

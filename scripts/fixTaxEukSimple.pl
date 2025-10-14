#! /usr/bin/env perl

## Clean up taxonomy file, first removing all -_ annotations, then filling in missing taxonomic levels with "Unclassified" labels.
## Additionally, if some genera are found within multiple families, the genus level is removed to avoid ambiguity.

use strict;
use warnings;

open(my $fh, '<', $ARGV[0]) or die "Could not open input file: $!";

my $outfile = $ARGV[1] // "output.txt";
open(my $out, '>', $outfile) or die "Could not open output file: $!";

print $out join("\t", "contigs", "predictions"), "\n";

my %fam_genus;
my @output_lines;

while (my $line = <$fh>) {

    my %taxlevels = map { $_ => "" } qw(k p c o f g s);

    chomp $line;
    $line =~ s/;-_[^;]+//g;
    $line =~ s/-_[^;]+;//g;
    $line =~ s/-_[^;]+$//g;

    my @elements = split(/\t/, $line);
    if(@elements == 1){
        $elements[1] = "unknown";
    }
    next if $elements[1] eq 'predictions';

    my @taxonomy = split(/;/, $elements[1]);
    
    for my $taxon (@taxonomy) {
        if ($taxon !~ /^-/) {
            my $lev = $taxon;
            $lev =~ s/^([kpcofgs])_.*$/$1/;
            $taxlevels{$lev} = $taxon;
        }
    }

    my $fam_genus_conc = join(";", @taxlevels{qw(f g)});

    $fam_genus{$fam_genus_conc} = $taxlevels{g};
    
    push @output_lines, join("\t", @elements). "\n";

}

my %valcount;
$valcount{$_}++ for values %fam_genus;

# Find which values are duplicates
my @dupes = grep { $valcount{$_} > 1 } keys %valcount;

for my $i (@output_lines){
    for my $d (@dupes){
        next if $d eq "";
        if ($i =~ /$d/){
            $i =~ s/;g_.*$//;
        }
    }

    print $out $i;
}

close($fh);
close($out);
#!/usr/bin/env perl

use Bio::TreeIO;
use strict;
use warnings;
use Data::Dumper;

die "$0: <tree file> <name.list>" unless @ARGV == 2;

my $in = Bio::TreeIO->new(-file=>shift @ARGV);
my $out = Bio::TreeIO->new();
my $tree = $in->next_tree();

my $file = shift @ARGV;
my (%name);
open IN,  "<" . $file;
while (<IN>){
    if (/^\s+.+(S\d+).+=>\s+'(\S+)',*\s*$/){
	$name{$1} = $2;
    }
}
close IN;

my @nodes = $tree->get_nodes();
foreach my $nd (@nodes){
    $nd->id($name{$nd->id}) if $nd->is_Leaf;
}

#print Dumper(\%name);

$out->write_tree($tree);

exit;

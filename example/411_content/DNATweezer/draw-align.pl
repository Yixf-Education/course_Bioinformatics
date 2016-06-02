#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use Bio::AlignIO;
use Bio::Align::Graphics;
use Data::Dumper;

my $file = shift;
my $in   = new Bio::AlignIO( -file => $file, -format => 'clustalw' );
my $aln  = $in->next_aln();
my ( @domain_start, @domain_end, @domain_color );
for ( my $i = 1; $i < $aln->length(); $i += 6 ) {
    push @domain_start, $i;
    push @domain_end,   $i + 2;
    push @domain_color, 'gray';
}

print STDERR Dumper( \@domain_start );

my $print_align = Bio::Align::Graphics->new(
    align      => $aln,
    pad_bottom => 5,

    out_format         => "png",
    block_space        => 1,
    block_size         => 30,
    reference          => 1,
    show_nonsynonymous => 1,
    wrap               => 102
);

$print_align->draw();

exit;

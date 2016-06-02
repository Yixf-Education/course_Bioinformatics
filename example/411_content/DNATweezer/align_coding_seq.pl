#!/usr/bin/perl -w

# -- OLD DOCUMENTATION --
# Looks like this is no longer true...
# From CLUSTALW documentation: #define MAXCHAR is now set to 30
# ----------------
# Align cds according to aa alignment
# Use bioperl
# BUG warning: CLUSTALW has id length restrictions (clips off the end makes names non-consistent between *.aa.aln and *.cds.fas files)
# Solution: Restore names in *.aa.aln before executing this script
# ----------------
# // OLD DOCUMENTATION //

use strict;
use warnings;
use 5.010;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::AlignIO;
use Bio::SeqIO;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;

my %opts = ();
GetOptions( \%opts, "help|h", "man|m", "debug", "fasta-alignment|a" )
    or pod2usage(2);

pod2usage(1) if $opts{"help"};
pod2usage( -exitstatus => 0, -verbose => 2 ) if $opts{"man"};

pod2usage(2)
    if ( @ARGV < 2 )
    ;    # Show a short usage statement if missing required arguments
my ( $pro, $cds ) = @ARGV;

my $cds_in = Bio::SeqIO->new( -file => "<$cds", '-format' => 'Fasta' )
    || die "Error opening file $cds.\n";
my %seqs
    ; # Will be a hash of display_id=>Bio::Seq object key/value pairs (needed for aa_to_dna_aln)

while ( my $seq = $cds_in->next_seq() ) {

    #	my $label = $seq->display_id();
    #	$label =~ s/:/_/g;  # Not needed??
    #	$seq->display_id($label);
    my $seq_string = $seq->seq();
    $seq_string =~ tr/[A-Z]/[a-z]/;
    $seq->seq($seq_string);
    $seqs{ $seq->display_id } = $seq;
}

my $aa_in = Bio::AlignIO->new( -file => "<$pro", '-format' => 'clustalw' )
    || die "Error opening file $pro for reading.\n";
my $cds_aln_out = Bio::AlignIO->new( '-format' => 'clustalw' );

# For optional FASTA, same output name, add .fas extension
my $cds_fasaln_out;
if ( $opts{'fasta-alignment'} ) {
    $cds_fasaln_out = Bio::AlignIO->new(
        -file     => ">$pro" . "cds.fas",
        '-format' => 'Fasta'
    ) || die "Error opening file ${pro}.fas for writing.\n";
}

while ( my $aln = $aa_in->next_aln ) {
    my $dna_align = aa_to_dna_aln( $aln, \%seqs );
    $dna_align->set_displayname_flat()
        ;    # Used in example: removes start/end info
    $cds_aln_out->write_aln($dna_align);

    if ( $opts{'fasta-alignment'} ) {
        $cds_fasaln_out->write_aln($dna_align);
    }
}

__END__

=head1 NAME

align_coding_seq.pl - Back aligns amino acid sequence alignments to nucleotide
sequence alignments. 

=head1 SYNOPSIS

B<align_coding_seq.pl> [OPTIONS] <CLUSTALW aa alignment file> <FASTA nt seq file> > foo.nuc.aln

=head1 DESCRIPTION

align_coding_seq takes as arguments an amino acid sequence alignment in CLUSTALW
format and a multi-FASTA file of nucleotide sequences (in that order) and
produces the corresponding nucleotide sequence alignment in CLUSTALW format. The 
sequence display ids in both input files must agree (they must be corresponding 
sequences) in order for output to be correct.

As an option, align_coding_seq.pl can also provide a multi-FASTA file of the
nucleotide alignment, which is named <input alignment>.cds

The output nucleotide alignment is printed to STDOUT.

=head1 OPTIONS

=over 4

=item B<--help, -h>

Print a brief help message and exits.

=item B<--man, -m>

Prints the man page and exits.

=item B<--debug>

Debugging output is printed to STDERR. Same effect if the environment variable
$DEBUG is defined and nonzero.

=item B<--fasta-alignment, -a>

Also produce a FASTA file of the alignment. The filename is <input alignment>.cds

=back

=head1 SEE ALSO

 none

=head1 AUTHOR

YE<246>zen HernE<225>ndez (or Yozen Hernandez, if rendering is wrong) (yzhernand at gmail dot com)

=cut

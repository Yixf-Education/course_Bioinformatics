#!/usr/bin/perl 

use strict;
use warnings;

use Getopt::Std;
use Bio::SeqIO;

&main;
exit;

sub main {
    &usage if ( @ARGV < 1 );
    my $command  = shift(@ARGV);
    my %function = (
        format  => \&format,
        revcom  => \&revcom,
        length  => \&length,
        content => \&content,
        search  => \&search,
        subseq  => \&subseq,
        random  => \&random
    );
    die("Unknown command \"$command\".\n")
      if ( !defined( $function{$command} ) );
    &{ $function{$command} };
}

sub format {
    my %opts = ( w => 60 );
    getopts( 'i:o:w:', \%opts );
    die(
        qq/
Usage:    seqTools.pl format [options]

Options:  -i STR   Input file in FASTA format.
          -o STR   Output file in FASTA format.
          -w INT   The line width for FASTA output. [$opts{w}]
\n/
    ) unless ( exists $opts{i} && exists $opts{o} );

    my ( $fi, $fo, $width ) = ( $opts{i}, $opts{o}, $opts{w} );
    my $in  = Bio::SeqIO->new( -file => "$fi",  -format => 'fasta' );
    my $out = Bio::SeqIO->new( -file => ">$fo", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        $out->width($width);
        $out->write_seq($seq);
    }
}

sub revcom {
    my %opts = ( w => 60 );
    getopts( 'i:o:w:', \%opts );
    die(
        qq/
Usage:    seqTools.pl revcom [options]

Options:  -i STR   Input file in FASTA format.
          -o STR   Output file in FASTA format.
          -w INT   The line width for FASTA output. [$opts{w}]
\n/
    ) unless ( exists $opts{i} && exists $opts{o} );

    my ( $fi, $fo, $width ) = ( $opts{i}, $opts{o}, $opts{w} );
    my $in  = Bio::SeqIO->new( -file => "$fi",  -format => 'fasta' );
    my $out = Bio::SeqIO->new( -file => ">$fo", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        my $rc_seq = $seq->revcom;
        $out->width($width);
        $out->write_seq($rc_seq);
    }
}

sub length {
    my %opts;
    getopts( 'i:', \%opts );
    die(
        qq/
Usage:    seqTools.pl length [options]

Options:  -i STR   Input file in FASTA format.
\n/
    ) unless ( exists $opts{i} );

    my $fi = $opts{i};
    my $in = Bio::SeqIO->new( -file => "$fi", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        print $seq->id, "\t", $seq->length, "\n";
    }
}

sub content {
    my %opts;
    getopts( 'i:', \%opts );
    die(
        qq/
Usage:    seqTools.pl content [options]

Options:  -i STR   Input file in FASTA format.
\n/
    ) unless ( exists $opts{i} );

    my $fi = $opts{i};
    my $in = Bio::SeqIO->new( -file => "$fi", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        print $seq->id, "\t";
        my $length   = $seq->length;
        my $sequence = $seq->seq;
        my $num_a    = $sequence =~ tr/A//;
        my $num_c    = $sequence =~ tr/C//;
        my $num_g    = $sequence =~ tr/G//;
        my $num_t    = $sequence =~ tr/T//;
        my $num_o    = $length - $num_a - $num_c - $num_g - $num_t;
        my $gc       = sprintf( "%.2f", ( $num_c + $num_g ) / $length * 100 );
        my $per_a    = sprintf( "%.2f", $num_a / $length * 100 );
        my $per_c    = sprintf( "%.2f", $num_c / $length * 100 );
        my $per_g    = sprintf( "%.2f", $num_g / $length * 100 );
        my $per_t    = sprintf( "%.2f", $num_t / $length * 100 );
        my $per_o    = sprintf( "%.2f", $num_o / $length * 100 );
        print
"Length=$length\tGC=${gc}%\tA=$num_a(${per_a}%), C=$num_c(${per_c}%), G=$num_g(${per_g}%), T=$num_t(${per_t}%), Others=$num_o(${per_o}%)\n";
    }
}

sub search {
    my %opts;
    getopts( 'i:p:', \%opts );
    die(
        qq/
Usage:    seqTools.pl search [options]

Options:  -i STR   Input file in FASTA format.
          -p STR   Pattern you want to search.

Notes:    1. You can use RegExp supported by Perl.
          2. It uses case-sensitive mode.
\n/
    ) unless ( exists $opts{i} && $opts{p} );

    my ( $fi, $pattern ) = ( $opts{i}, $opts{p} );
    my $in = Bio::SeqIO->new( -file => "$fi", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        print $seq->id, "\t";
        my @starts;
        my $sequence = $seq->seq;
        while ( $sequence =~ m/($pattern)/g ) {
            my $start = pos($sequence) - CORE::length($1) + 1;
            pos($sequence) = $start;
            push @starts, $start;
        }
        my $num = @starts;
        my $str = join ",", @starts;
        print "$num\t$str\n";
    }
}

sub subseq {
    my %opts = ( n => ".+", s => 1, e => 0, );
    getopts( 'i:o:n:s:e:', \%opts );
    die(
        qq/
Usage:    seqTools.pl subseq [options]

Options:  -i STR   Input file in FASTA format.
          -o STR   Output file in FASTA format.
          -n STR   Name|ID of the FASTA record. [$opts{n}]
          -s INT   Start of the subsequence. [$opts{s}]
          -e INT   End of the subsequence. [$opts{e}]

Notes:    1. You can use RegExp supported by Perl in the Name|ID (-n Name|ID) to match more than one record. But the RegExp must be in the double quotation marks, for example, "chr\\d+".
          2. The RegExp use case-sensitive mode.
          3. "$opts{n}" for the Name|ID means all records in the input file.
          4. Start and End use the 1-based coordinate system, where the first base is 1 and the number is inclusive, ie 1-2 are the first two bases of the sequence.
          5. "$opts{e}" for the End means to the end of the record, that is to say, the End will be set to the Length.
          6. If Start OR End is larger than the length of the record, they will be set to the Length.
\n/
    ) unless ( exists $opts{i} && exists $opts{o} );

    my ( $fi, $fo, $name, $start, $end ) =
      ( $opts{i}, $opts{o}, $opts{n}, $opts{s}, $opts{e} );
    my $in  = Bio::SeqIO->new( -file => "$fi",  -format => 'fasta' );
    my $out = Bio::SeqIO->new( -file => ">$fo", -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {
        my $id   = $seq->id;
        my $desc = $seq->desc;
        my $subseq;
        if ( $start > $seq->length ) {
            $start = $seq->length;
        }
        if ( CORE::length $desc > 0 ) {
            $desc .= " ";
        }
        if ( $id =~ /^$name$/ ) {
            if ( $end == 0 || $end > $seq->length ) {
                $end = $seq->length;
                $subseq = $seq->subseq( $start, $end );
                if ( $start == 1 ) {
                    $desc .= "fullseq:$start-$end";
                }
                else {
                    $desc .= "subseq:$start-$end";
                }
                $end = 0;
            }
            else {
                $subseq = $seq->subseq( $start, $end );
                $desc = $desc . "subseq:$start-$end";
            }
            my $record = Bio::Seq->new(
                -display_id => "$id",
                -desc       => "$desc",
                -seq        => "$subseq",
            );
            $out->write_seq($record);
        }
    }
}

sub random {
    use Tie::File::AnyData::Bio::Fasta;
    use Fcntl qw/O_RDONLY O_RDWR O_CREAT/;
    use List::Util qw/shuffle/;

    my %opts = ( n => 1 );
    getopts( 'i:o:n:', \%opts );
    die(
        qq/
Usage:    seqTools.pl random [options]

Options:  -i STR   Input file in FASTA format.
          -o STR   Output file in FASTA format. [STDOUT]
          -n INT   The number of records. [$opts{n}]
\n/
    ) unless ( exists $opts{i} );

    my ( $fi, $n ) = ( $opts{i}, $opts{n} - 1 );
    my ( $fo, @out );
    tie my @in, 'Tie::File::AnyData::Bio::Fasta', $fi, mode => O_RDONLY
      or die $!;
    if ( $opts{o} ) {
        $fo = $opts{o};
        tie @out, 'Tie::File::AnyData::Bio::Fasta', $fo,
          mode => O_RDWR | O_CREAT
          or die $!;
    }
    @out = ( shuffle @in )[ 0 .. $n ];
    untie @in;
    if ( $opts{o} ) {
        untie @out;
    }
    else {
        foreach (@out) {
            print "$_\n";
        }
    }
}

sub usage {
    die(
        qq/
Usage:   seqTools.pl <command> [<arguments>]

Command: format   Format the FASTA record(s).
         revcom   Reverse and complement the FASTA record(s).
         length   Get the length of FASTA record(s).
         content  Calculate the GC content of FASTA record(s).
         search   Find subseq in the FASTA record(s).
         subseq   Extract sub|full-seq from the FASTA record(s).
         random   Take record(s) from a FASTA file randomly.

Prerequisites:    1. Getopt::Std. [ALL]
                  2. Bio::Seq. [except RANDOM]
                  3. Tie::File::AnyData::Bio::Fasta, Fcntl, List::Util. [RANDOM]

Author:  Yixf, yixf1986\@gmail.com
Version: v3.1

Notes:  1. A similar webtool developed by lh3: http:\/\/lh3lh3.users.sourceforge.net\/fasta.shtml
        2. When you search, you can use RegExp supported by Perl.
        3. When you extract the sub|full-seq, you can use RegExp in the Name|ID (-n Name|ID) to match more than one record. But the RegExp must be in the double quotation marks, for example, "chr\\d+".
        4. The RegExp use case-sensitive mode.

Versions: v1.0, 2011-06-01
          v2.0, 2011-11-03
          v3.0, 2012-02-27
		  v3.1, 2013-08-22
\n/
    );
}

#!/usr/bin/perl -w

# This script accepts fasta alignment files as arguments and generates a
# concatenated alignment; An id in one alignment file must match an id in
# another alignment file in order for the two sequences to be concatenated.
# This script is designed to accept alignment files with unique sequences,
# generating gaps in the concatenated output.

use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Data::Dumper;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;

my $output_format = 'clustalw';
my $output_filename = "concat_output";
my $delim = '|';
my %opts = ();
GetOptions (\%opts, "help|h",
					"man|m",
					"debug",
					"format|f=s" => \$output_format,
					"delim|d=s" => \$delim,
					"outfile|o=s" => \$output_filename,
                    "ref|r=s",
					"safe|safe-names|s",
					"lc:1",
					"uc") or pod2usage(2);

pod2usage(1) if $opts{"help"};
pod2usage(-exitstatus => 0, -verbose => 2) if $opts{"man"};

pod2usage(2) if (@ARGV == 0); # Show a short usage statement if missing required arguments

my %seq_ids; # used to be %taxon, used to keep track of *all* seq ids
my @aln_array; # This will be an array of Bio::SimpleAlign objects
my $concat_align;
my $lowercase = ($opts{"uc"}) ? 0 : $opts{"lc"};
my $unsafe_chars = "[-.\_]"; #This is a regex pattern, so some characters need to be escaped.

#At this point, @ARGV should only consist of filenames.

# <<Loop through command line args and populate @aln_array with alignments>>
my $aln_count = 0; # Assign a unique id for each alignment
while ( my $filename = shift @ARGV ) {
	my $in = Bio::AlignIO->new('-file' => $filename);
	my $lseq_array;
	print STDERR "$filename\n" if ($ENV{DEBUG} || $opts{'debug'});
	
	while ( my $aln_obj=$in->next_aln() ) {
		# Keep track of all sequence ids.
		foreach my $seq ($aln_obj->each_seq) {
			# Try to find a unique part of the id, assuming it consists of fields.
			# FIXME This this id parsing is not universal
			my ($unique_id) = split /\Q$delim\E/, $seq->id;
			$unique_id =~ s/$unsafe_chars//g if ($opts{'safe'});
			
			# We copy the sequence and create a new LocatableSeq with the unique id
			my $new_lseq = Bio::LocatableSeq->new('-id' => $unique_id, '-seq' => ($lowercase != 0) ? lc($seq->seq) : uc($seq->seq));
			
			# Add sequence id to the seq_ids hash if we've never seen it before
			$seq_ids{$unique_id} = 1 if ( !(exists $seq_ids{$unique_id}) );
			
			# Push it onto an array reference to place in a new SimpleAlign object
			push @$lseq_array, $new_lseq;
		}
		
		my $new_aln = Bio::SimpleAlign->new('-seqs' => $lseq_array);
		push @aln_array, $new_aln;
	}
	$aln_count++; # Increment counter for next alignment.
}

# By this point, we know all seq_ids that should ever appear in these files.

### DEBUGGING STATEMENTS ###
print STDERR "Sequences in last aln: ",Dumper( $aln_array[$aln_count-1] ) if ($ENV{DEBUG} || $opts{'debug'});
print STDERR "The aln_array: ",Dumper(\@aln_array) if ($ENV{DEBUG} || $opts{'debug'});
print STDERR "All seq_ids: ",Dumper(sort keys %seq_ids) if ($ENV{DEBUG} || $opts{'debug'});
############################

# <<Check all alignments for 'missing' sequences>>
$aln_count = 0; # Using this counter again, so reset to 0
foreach my $aln ( @aln_array ) {
	print STDERR "Aln id $aln_count\n" if ($ENV{DEBUG} || $opts{'debug'});
	
	# Get all seq_ids and save them in a hash with the id as key and seqobj as value. (makes lookup fast)
	my %seq_id_hash = map { $_->id => $_ } $aln->each_seq;
	
	print STDERR "seq_id_hash\n", Dumper(\%seq_id_hash) if ($ENV{DEBUG} || $opts{'debug'});
	
	# <<Add gap sequences for any missing sequences>>
	# Check each known seq_id
	foreach my $seqid ( keys %seq_ids ) {
		# If the seq_id is missing in this alignment, add a "missing data" sequence of the same length
		if ( !(exists $seq_id_hash{$seqid}) ) {
			my $missing = '?' x $aln->length;
			# New sequence with '?'s
			my $new_seq = Bio::LocatableSeq->new( '-seq' => $missing, '-id' => $seqid);
			$aln->add_seq($new_seq);
			$seq_id_hash{$seqid} = $new_seq;
		}
	}
	
#### At this point, this alignment will have had any missing sequences filled. ####
	
	# <<Concatenate all sequences together into one large alignment>>
	# If its the first alignment, it will be the one to which we concatenate the
	# remaining sequences, so we save it.
	if ( $aln_count == 0 ) {
		$concat_align = $aln;
	} else { # Else, we concatenate other sequences to those in concat_align and replace it
		my $lseq_array = [];
		foreach my $concat_seq ($concat_align->each_seq) {
			my $id = $concat_seq->id;
			my $old_seq = $concat_seq->seq;
			
			print STDERR "display id: ", $id, "\n" if ($ENV{DEBUG} || $opts{'debug'});
			
			my $seq_to_add = $seq_id_hash{$id}->seq; # We use the hash of the current alignment
			my $new_seq_obj = Bio::LocatableSeq->new( '-seq' => $old_seq . $seq_to_add, '-id' => $id);
			push @$lseq_array, $new_seq_obj;
		}
		
		$concat_align = Bio::SimpleAlign->new('-seqs' => $lseq_array);
	}
	
	$aln_count++;
}

print STDERR Dumper( $concat_align ) if ($ENV{DEBUG} || $opts{'debug'});

# <<Write out alignment to requested/default format>>
my $out_aln = Bio::AlignIO->new( -file => ">$output_filename", -format => $output_format);
$concat_align->set_displayname_flat(); # Suppress positions
#print STDERR "Setting new reference: ", $opts{"ref"}, "\n";
#print Dumper($concat_align);
$concat_align = $concat_align->set_new_reference($opts{"ref"}) if ($opts{"ref"});
#print Dumper($concat_align);
$out_aln->write_aln($concat_align);

exit;


__END__

=head1 NAME

align_concatenate.pl

=head1 SYNOPSIS

B<align_concatenate.pl> [OPTIONS] <input file> [input files...]

=head1 DESCRIPTION

align_concatenate.pl accepts alignment files as arguments and generates a
concatenated alignment in a single file. Sequence ids in the alignment files
must contain a matching portion across all input files in order for the two
sequences to be concatenated. For example, if two corresponding sequences
beloning to Homo sapiens were to be concatenated, the ids in each alignment
could be:
 >h._sapiens

or
 >h._sapiens|chr1
 >h._sapiens|chr2

where the matching portion is the beginning of the id and the '|' acts as a
delimiter. The delimiter can be changed on the command line with the -d switch.

If an alignment is missing some sequence, then it is added and the sequence data
is filled with '?'s (missing data, in mrbayes).

Currently, expected input is clustalw, although align_concatenate.pl will
happily accept any format known to Bio::AlignIO. It has only been tested with
clustalw input.  

align_concatenate.pl will take the output file format as an option on the
command line. The default is clustalw format.

=head1 OPTIONS

=over 4

=item B<--help, -h>

Print a brief help message and exits.

=item B<--man, -m>

Prints the man page and exits.

=item B<--debug>

Debugging output is printed to STDERR. Same effect if the environment variable
$DEBUG is defined and nonzero.

=item B<--delim, -d> '<character>'

Change the sequence id delimiter to '<character>'. By default, the delimiter is
'|'.

=item B<--format, -f>

Output file format. By default, this is "clustalw"

Specify output format.

=item B<--outfile, -o>

File name or full pathname of output. 

Specify the name/full path of the ouput file. Default is a file called
'concat_output' in the current directory.

=item B<--ref, -r> 'seq_id'

Change the reference sequence to be seq_id. seq_id must be one of the unique ids
or id parts required by this script.

=item B<--safe, --safe-names, -s>

If this option is given on the command line, output filenames and sequence
display id's will be "sanitized" so that no characters with possibly special
meaning in other programs will be included. So far, characters which this option
will remove include: ., -, and _.

=item B<--lc, --uc>

Lowercase and uppercase, respectively. If either of these options is given, then
sequence characters are converted to either lowercase or uppercase.

The default is lowercase (NCBI standard for nucleotide sequences). --uc
overrides --lc, so that if both options are given, then the sequences are all
uppercase.

=back

=head1 RESTRICTIONS (Bugs we won't/can't fix)

Even when specifying a unique portion of the input sequence name, this part won't
be displayed in the concatenated alignment file. Instead, whatever name the taxon
had in the first alignment supplied will be used.

=head1 TODO

[Priority = Low] Add more support to notations for different formats (ie NEXUS)

=head1 SEE ALSO

 none

=head1 AUTHOR

James Haven

YE<246>zen HernE<225>ndez (or Yozen Hernandez, if rendering is wrong) (yzhernand at gmail dot com)

=cut

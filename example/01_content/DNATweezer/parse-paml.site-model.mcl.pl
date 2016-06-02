#!/usr/bin/perl -w
#
# parse PAML4 site model output
#

use strict;
use warnings;
use Bio::Seq;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Std;

my %opts;
getopts('hrwb', \%opts);

if ($opts{h}) {
  print "\tUsage:$0 [-hrwb] <paml site-model mlc file>\n";
  print "\t-h:\t print help\n";
  print "\t-r:\t print to result table\n";
  print "\t-w:\t print to omega table\n";
  print "\t-b:\t print to bayes table\n";
  exit;
}

my $infile = shift @ARGV;
my $start_seq = 0;
my $model_0 = 0;
my $model_1 = 0;
my $model_2 = 0;
my $in_bayes = 0;
my @bayes;

my $result_0 = {
		'lnL'=>undef,
		'np'=>undef,
		'tree_length_dN'=>undef,
		'tree_length_dS'=>undef,
		'kappa'=>undef,
		'omega'=>undef,
	       };

my $result_1 = {
		'lnL'=>undef,
		'np'=>undef,
		'tree_length'=>undef,
		'kappa'=>undef,
		'p0' => undef,
		'w0' => undef,
		'p1' => undef,
		'w1' => undef,
	       };

my $result_2 = {
		'lnL'=>undef,
		'np'=>undef,
		'tree_length'=>undef,
		'tree_string'=>undef,
		'kappa'=>undef,
		'p0' => undef,
		'w0' => undef,
		'p1' => undef,
		'w1' => undef,
		'p2' => undef,
		'w2' => undef,
		'bayes' => undef,
	       };

my $result = {
	      'model_0' => $result_0,
	      'model_1' => $result_1,
	      'model_2' => $result_2,
	      'n_tax' => undef,
	      'n_codon' => undef,
	      'ref_seq' => undef,
	     };

open (my $in, "<", $infile) or die "can't open file: $!";

while (<$in>) {

  if (/^seed used/) {
    $start_seq = 1;
    next;
  }

  if (/^\s+(\d+)\s+(\d+)\s*$/ && $start_seq) {
    $result->{n_tax} = $1;
    $result->{n_codon} = $2;
   }

  if (/^(\S+)\s+(.+)\s*$/ && $start_seq) {
    my ($id, $str) = ($1, $2);
    $str =~ s/\s//g;
    $result->{ref_seq} = Bio::Seq->new(-id=>$id, -seq=>$str);
    $start_seq = 0;
  }

  $model_0 = 1 if /^Model 0: one-ratio/;
  if ($model_0) {
    if (/^lnL.+np:\s+(\d+)\):\s+(-\d+\.\d+)\s+\S+\s*$/) { 
      $result_0->{np} = $1;
      $result_0->{lnL} = $2;
    }

    if (/^kappa\s+.+=\s+([\d\.]+)\s*$/) {
      $result_0->{kappa} = $1;
    }

    if (/^omega\s+.+=\s+([\d\.]+)\s*$/) {
      $result_0->{omega} = $1;
    }

    if (/^tree length for dN:\s+([\d\.]+)\s*$/) {
      $result_0->{tree_length_dN} = $1;
    }

    if (/^tree length for dS:\s+([\d\.]+)\s*$/) {
      $result_0->{tree_length_dS} = $1;
      $model_0 = 0;
    }
  }

  $model_1 = 1 if /^Model 1: NearlyNeutral/;
  if ($model_1) {
    if (/^lnL.+np:\s+(\d+)\):\s+(-\d+\.\d+)\s+\S+\s*$/) { 
      $result_1->{np} = $1;
      $result_1->{lnL} = $2;
    }

    if (/^kappa\s+.+=\s+([\d\.]+)\s*$/) {
      $result_1->{kappa} = $1;
    }

    if (/^tree length\s+=\s+([\d\.]+)\s*$/) {
      $result_1->{tree_length} = $1;
    }

    if (/^p:\s+([\d\.]+)\s+([\d\.]+)\s*$/) {
      $result_1->{p0} = $1;
      $result_1->{p1} = $2;
    }

    if (/^w:\s+([\d\.]+)\s+([\d\.]+)\s*$/) {
      $result_1->{w0} = $1;
      $result_1->{w1} = $2;
      $model_1 = 0;
    }
  }

  $model_2 = 1 if /^Model 2: PositiveSelection/;
  if ($model_2) {
    if (/^lnL.+np:\s+(\d+)\):\s+(-\d+\.\d+)\s+\S+\s*$/) { 
      $result_2->{np} = $1;
      $result_2->{lnL} = $2;
    }

    if (/^kappa\s+.+=\s+([\d\.]+)\s*$/) {
      $result_2->{kappa} = $1;
    }

    if (/^tree length\s+=\s+([\d\.]+)\s*$/) {
      $result_2->{tree_length} = $1;
    }

    if (/^\(\D.+\);\s*$/) {
      chomp;
      $result_2->{tree_string} = $_;
    }

    if (/^p:\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s*$/) {
      $result_2->{p0} = $1;
      $result_2->{p1} = $2;
      $result_2->{p2} = $3;
    }

    if (/^w:\s+([\d\.]+)\s+(1\.00000)\s*([\d\.]+)\s*$/) {
      $result_2->{w0} = $1;
      $result_2->{w1} = $2;
      $result_2->{w2} = $3;
    }

    if (/^Bayes Empirical Bayes/) {
      $in_bayes = 1; # print "getting into the bayes block ....\n";
    }

    if (/^The grid/) {
      $result_2->{bayes} = \@bayes; # print ".........getting out of bayes.\n";
      $in_bayes = 0;
      $model_2 = 0;
    }

    if (/^\s+(\d+)\s+(\D)\s+([\d\.]+)\**\s+([\d\.]+)\s+\S\S\s+([\d\.]+)\s*$/ && $in_bayes) {
      my ($pos, $res, $post, $mean, $se) = ($1, $2, $3, $4, $5);
      my $site = {
		  'position' => $pos,
		  'residue' => $res,
		  'post_prob' => $post,
		  'w' => $mean,
		  'w_se' => $se,
	   };
      push @bayes, $site; 
    }
  }
}

close $in;
#print Dumper($result); exit;
&print_result_table() if $opts{r};
&print_omega_table() if $opts{w};
&print_bayes_table() if $opts{b};
exit;

sub print_result_table {
  print join "\t", (
		    $infile,
		    1, # model 0
		    $result->{model_0}->{lnL}, 
		    $result->{ref_seq}->id(),
		    $result->{n_tax}, 
		    $result->{n_codon}, 
		    $result->{model_0}->{np},
		    $result->{model_0}->{kappa},
		    $result->{model_0}->{tree_length_dN}+$result->{model_0}->{tree_length_dS},
		  );
  print "\n";
  
  print join "\t", (
		    $infile,
		    2, # model 0
		    $result->{model_1}->{lnL}, 
		    $result->{ref_seq}->id(),
		    $result->{n_tax}, 
		    $result->{n_codon}, 
		    $result->{model_1}->{np},
		    $result->{model_1}->{kappa},
		    $result->{model_1}->{tree_length},
		  );
  print "\n";

  print join "\t", (
		    $infile,
		    3, # model 0
		    $result->{model_2}->{lnL}, 
		    $result->{ref_seq}->id(),
		    $result->{n_tax}, 
		    $result->{n_codon}, 
		    $result->{model_2}->{np},
		    $result->{model_2}->{kappa},
		    $result->{model_2}->{tree_length},
		  );
  print "\n";
}

sub print_omega_table {
  print join "\t", ( # one w for model 0
		    $infile,
		    1, # model class
		    0, # negative, one-ratio
		    't', # all foreground
		    1.0, # 100%
		    $result->{model_0}->{omega},
		  );
  print "\n";
  
  print join "\t", ( # 2 w's for model 1
		    $infile,
		    2, # model 1
		    0, # negative
		    't',
		    $result->{model_1}->{p0},
		    $result->{model_1}->{w0},
		  );
  print "\n";

  
  print join "\t", (
		    $infile,
		    2, # model 1
		    1, # neutral
		    't',
		    $result->{model_1}->{p1},
		    $result->{model_1}->{w1},
		  );
  print "\n";

  print join "\t", ( # 3 w's for model 2
		    $infile,
		    3, # model 2
		    0, # negative
		    't',
		    $result->{model_2}->{p0},
		    $result->{model_2}->{w0},		    
		  );
  print "\n";

  print join "\t", ( # 3 w's for model 2
		    $infile,
		    3, # model 2
		    1, # neutral
		    't',
		    $result->{model_2}->{p1},
		    $result->{model_2}->{w1},		    
		  );
  print "\n";

  print join "\t", ( # 3 w's for model 2
		    $infile,
		    3, # model 2
		    2, # positive
		    't',
		    $result->{model_2}->{p2},
		    $result->{model_2}->{w2},		    
		  );
  print "\n";

}

sub print_bayes_table {
  exit unless defined @{$result->{model_2}->{bayes}};
  for my $site (@{$result->{model_2}->{bayes}}) { 
    print join "\t", (
		      $infile,
		      3, # model 2
		      $site->{position}, 
		      $site->{residue}, 
		      $site->{post_prob},
		      $site->{w},
		      $site->{w_se},   
		     );
    print "\n";
  }
}


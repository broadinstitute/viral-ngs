#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Easel::MSA;
use Bio::Easel::SqFile;
require "sqp_seq.pm";
require "sqp_utils.pm";

my $usage;
$usage  = "fasta-trim-terminal-ambigs.pl\n\n";
$usage .= "Usage:\n\n";
$usage .= "perl fasta-trim-terminal-ambigs.pl [OPTIONS] <fasta file>\n\n";
$usage .= "OPTIONS:\n";
$usage .= "  --minlen <n> : min allowed sequence length after trimming [default (df): 1]\n";
$usage .= "  --maxlen <n> : max allowed sequence length after trimming [df: 1Gb]\n";
$usage .= "  --sfx <s>    : suffix to add to each sequence name is <s> [df: do not change names]\n";
$usage .= "  --strict     : die (instead of skipping a seq) if any seq is not within minlen..maxlen range after trimming\n";
$usage .= "\nAlternatively, use the --3rules option to trim using a more sophisticated method that removes\n";
$usage .= "terminal nucleotides (nts) using the following three rules:\n\n";

$usage .= "Rule 1. remove the 10 terminal nts if > <n1> ambiguous nt exist in 10 terminal nt (<n1> from --ten <n1>)\n";
$usage .= "Rule 2. remove the 50 terminal nts if > <n2> ambiguous nt exist in 50 terminal nt (<n2> from --fifty <n2>)\n";
$usage .= "Continue to enforce rule 1 until it no longer applies, then enforce rule 2 and if rule 2 applies\n";
$usage .= "we go back and try to apply rule 1 again, and continue until neither rule 1 or 2 applies. Then\n";
$usage .= "move on to rule 3.\n";

$usage .= "Rule 3. remove all terminal ambiguous nt (such that sequence starts/ends with non-ambiguous nt)\n\n";
$usage .= "Then output sequence only if it has length between the minimum and maximum (from --minlen and\n";
$usage .= "--maxlen) and is less than <f> fraction Ns (<f> from --maxfrac <f>).\n\n";

$usage .= "Options related to the alternative '3 rules' strategy:\n";
$usage .= "  --3rules     : use the 3 rules strategy described above [df: do not]\n";
$usage .= "  --ten <n>    : max number of ambiguous nucleotides allowed in first/final 10 [df: 5]\n";
$usage .= "  --fifty <n>  : max number of ambiguous nucleotides allowed in first/final 50 [df: 15]\n";
$usage .= "  --maxfrac <f>: max allowed fraction of sequence that can be Ns after trimming [df: 0.5]\n";

# set defaults, for SC2 processing, GenBank used min of 50 and max of 30000 (as of May 2021)
my $minlen          = 1;
my $maxlen          = 1000000000; # 1Gb

my $sfx             = undef;
my $do_strict       = 0;
my $do_3rules       = 0; # set to 1 if --3rules used

# variables only relevant if --3rules is used
my $df_ten_max_ambig   = 5;
my $df_fifty_max_ambig = 15;
my $df_maxfrac_Ns      = 0.5;
my $ten_max_ambig      = undef;
my $fifty_max_ambig    = undef;
my $maxfrac_Ns         = undef;

&GetOptions( "minlen=s" => \$minlen,
             "maxlen=s" => \$maxlen,
             "sfx=s"    => \$sfx, 
             "strict"   => \$do_strict,
             "3rules"   => \$do_3rules,
             "ten=s"    => \$ten_max_ambig,
             "fifty=s"  => \$fifty_max_ambig,
             "maxfrac=s"=> \$maxfrac_Ns);

if(scalar(@ARGV) != 1) { die $usage; }
my ($fasta_file) = @ARGV;

# enforce --3rules related options only used if --3rules also used
if(! $do_3rules) {
  if(defined $ten_max_ambig) {
    die "ERROR, using --ten only makes sense if --3rules is also used";
  }
  if(defined $fifty_max_ambig) {
    die "ERROR, using --fifty only makes sense if --3rules is also used";
  }
  if(defined $maxfrac_Ns) { 
    die "ERROR, using --maxfrac only makes sense if --3rules is also used";
  }
}
if(! defined $ten_max_ambig) {
  $ten_max_ambig = $df_ten_max_ambig;
}
if(! defined $fifty_max_ambig) {
  $fifty_max_ambig = $df_fifty_max_ambig;
}
if(! defined $maxfrac_Ns) { 
  $maxfrac_Ns = $df_maxfrac_Ns;
}

# enforce ranges that make sense for options
if($minlen < 0) { 
  die "ERROR with --minlen <n>, <n> must be >= 0";
}
if($maxlen < 0) { 
  die "ERROR with --maxlen <n>, <n> must be >= 0";
}
if(($ten_max_ambig < 0) || ($ten_max_ambig >= 10)) { 
  die "ERROR with --ten <n>, <n> must in range [0..9]";
}
if(($fifty_max_ambig < 0) || ($fifty_max_ambig >= 50)) { 
  die "ERROR with --ten <n>, <n> must in range [0..49]";
}
if($minlen > $maxlen) { 
  die "ERROR with --minlen <n1> and --maxlen <n2>, <n1> must not be greater than <n2>";
}
if($ten_max_ambig > $fifty_max_ambig) { 
  die "ERROR with --ten <n1> and --fifty <n2>, <n1> must not be greater than <n2>";
}
if(($maxfrac_Ns < -0.00001) || ($maxfrac_Ns > 1.00001)) { 
  die "ERROR with --maxfrac <f>, <f> must be between 0 and 1";
}

if(! -s $fasta_file) { 
  die "ERROR fasta file $fasta_file does not exist or is empty";
}

# determine number of sequences in file
my $nseq = `grep ^\\> $fasta_file | wc -l`;
chomp $nseq;
$nseq =~ s/^\s+//; # remove leading whitespace
$nseq =~ s/\s+$//; # remove trailing whitespace
if($nseq !~ /^\s*\d+\s*$/) { 
  die "ERROR could not determine number of sequences in file using grep and wc";
}

# open Bio::Easel SqFile object
my $sqfile  = Bio::Easel::SqFile->new({ fileLocation => $fasta_file }); # the sequence file object

my $fasta_seq = undef;
my @fasta_seq_A = ();
my ($header, $sqstring, $out_header) = (undef, undef, undef);

# for each sequence, remove Ns and output to STDOUT
for(my $i = 0; $i < $nseq; $i++) { 
  $fasta_seq = $sqfile->fetch_consecutive_seqs(1, "", -1); # -1: unlimited line length
  @fasta_seq_A = split(/\n/, $fasta_seq);
  if(scalar(@fasta_seq_A) != 2) { 
    die "ERROR fetching sequence, could not split into header and sequence lines:\n$fasta_seq\n"; 
  }
  ($header, $sqstring) = (@fasta_seq_A);

  # add sfx to sequence name, if --sfx used
  if(defined $sfx) { 
    if($header =~ /^\>(\S+)(\s*.*)$/) { 
      $out_header = ">" . $1 . $sfx . $2;
    }
    else { 
      die "ERROR unable to parse header line:\n$header\n";
    }      
  }        
  else { # --sfx not used
    $out_header = $header;
  }

  # two modes:
  # ! do_3rules  (--3rules not used):  remove leading and trailing ambigs only
  # do_3rules    (--3rules):           enforce the three rules

  if(! $do_3rules) { 
    $sqstring =~ s/^[^ACGTUacgtu]+//; # remove any 5'-terminal ambiguous nts
    $sqstring =~ s/[^ACGTUacgtu]+$//; # remove any 3'-terminal ambiguous nts
  }
  else { 
    # --3rules used, enforce the 3 rules 
    # We trim the sequence at 5' end, then reverse it and trim at 3' end
    # (this makes it so we can reuse same code for both ends)
    # then reverse it back and output it (if there's any sequence left after trimming)
    $sqstring = trim_5p_end_using_three_rules($sqstring, $ten_max_ambig, $fifty_max_ambig);
    ##printf("5' trimmed length: %d\n", length($sqstring));
    if($sqstring ne "") { 
      $sqstring = reverse($sqstring);
      $sqstring = trim_5p_end_using_three_rules($sqstring, $ten_max_ambig, $fifty_max_ambig);
      ##printf("3' trimmed length: %d\n", length($sqstring));
      if($sqstring ne "") { # reverse it back to original forward direction
        $sqstring = reverse($sqstring);
      }
    }
  }
  my $seqlen = length($sqstring);
  my $n_N    = () = $sqstring =~ /[Nn]/g; # count Ns
  my $frac_N = ($seqlen > 0) ? ($n_N / $seqlen) : 0.;

  if(($seqlen < $minlen) ||     # too short (after trimming)
     ($seqlen > $maxlen) ||    # too long  (after trimming)
     (($do_3rules) && ($frac_N > $maxfrac_Ns))) { # --3rules used and too many Ns (after trimming)
    # do not output sequence
    if(! $do_strict) { 
      ; # --strict not used, it's ok to skip seqs, do nothing
    }
    else { 
      # --strict used, it's not ok to skip seqs, die
      if($seqlen < $minlen) { 
        die "ERROR sequence with the header line below is too short (after trimming) length is $seqlen minlen is $minlen.\n$header\n";
      }
      elsif($seqlen > $maxlen) { 
        die "ERROR sequence with the header line below is too short (after trimming) length is $seqlen minlen is $minlen.\n$header\n";
      }
      else { 
        die "ERROR sequence with the header line below has too many Ns (after trimming) fraction Ns is $frac_N, max allowed is $maxfrac_Ns\n$header\n";
      }
    }
  }
  else { 
    # $sqstring is not too short, not too long, and doesn't have too many Ns (if --3rules)
    print $out_header . "\n";
    print seq_SqstringAddNewlines($sqstring, 60);
  }
}

exit 0;

#################################################################
# Subroutine:  trim_5p_end_using_three_rules()
# Incept:      EPN, Wed Apr  7 07:22:43 2021
#
# Purpose:   Given a sequence string <$sqstring>, remove ambiguous
#            nts from the 5' end the same way that GenBank processing 
#            does, using 3 rules:
#            while any rule results in trimming:
#            - rule 1: remove the 10 5'-most nt if > $ten_max_ambig 
#                      ambiguous nt exist in 10 5'-most nt [5 by default]
#            - rule 2: remove the 50 5'-most nt if > $fifty_max_ambig 
#                      ambiguous nt exist in 50 5'-most nt [15 by default]
#            - rule 3: after rule 1 and rule 2 trimming is finished,
#                      trim terminal ambiguous nucleotides before returning.
#
# Arguments:
#   $sqstring:        the nucleotide sequence string
#   $ten_max_ambig:   max number of ambiguous nucleotides allowed in first 10
#   $fifty_max_ambig: max number of ambiguous nucleotides allowed in first 50
#
# Returns:    the input <$sqstring> trimmed as above, or "" if none is left after trimming
#
#################################################################
sub trim_5p_end_using_three_rules { 
  my $sub_name = "trim_5p_end_using_three_rules";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($sqstring, $ten_max_ambig, $fifty_max_ambig) = @_;

  my $keep_going = 1;
  my ($next_10, $next_50, $nambig, $trim_offset); 
  my $ntrimmed = 0;
  while(($keep_going) && ($sqstring ne "")) { 
    $keep_going = 0; # set to 1 if we use any rules below
    # rule 1: remove the 10 5'-most nt if > $ten_max_ambig ambiguous nt exist in 10 5'-most nt
    $next_10 = substr($sqstring, 0, 10);
    $nambig = () = $next_10 =~ /[^ACGTUacgtu]/g;
    if($nambig > $ten_max_ambig) { 
      # determine how many to trim, stop trimming at final ambig char in region
      # remove all non-ambiguous nt from the end to determine length to trim
      $next_10 =~ s/[ACGTUacgtu]+$//;
      $trim_offset = length($next_10);

      $ntrimmed += $trim_offset;
      ##printf("nambig: $nambig, trimming $trim_offset: %s (ntrimmed: $ntrimmed)\n", substr($sqstring, 0, $trim_offset));

      $sqstring = substr($sqstring, $trim_offset);
      $keep_going = 1;
    }
    # rule 2: remove the 50 5'-most nt if > $fifty_max_ambig ambiguous nt exist in 50 5'-most nt
    if(! $keep_going) { # only enforce this rule if we didn't use rule 1
      $next_50 = substr($sqstring, 0, 50);
      $nambig = () = $next_50 =~ /[^ACGTUacgtu]/g;
      if($nambig > $fifty_max_ambig) { 
        # determine how many to trim, stop trimming at final ambig char in region
        # remove all non-ambiguous nt from the end to determine length to trim
        $next_50 =~ s/[ACGTUacgtu]+$//;
        $trim_offset = length($next_50);

        $ntrimmed += $trim_offset;
        ##printf("nambig: $nambig, trimming $trim_offset: %s (ntrimmed: $ntrimmed)\n", substr($sqstring, 0, $trim_offset));
        $sqstring = substr($sqstring, $trim_offset);
        $keep_going = 1;
      }
    }
  }

  # finally, trim any leading (terminal) ambiguities
  $sqstring =~ s/^[^ACGTUacgtu]+//;

  return $sqstring;
}


#!/usr/bin/perl

# Compare a sequence to a reference, measure similarity etc. If variant calling 
# result is specified for the assembly, translate the coordinates onto the ref
# prerequisite: specify the path of muscle

use strict;
use Getopt::Long;
use Data::Dumper;
 
my %option = (
	h         		=> '',
	ref						=> '',
	seq						=> '',
	oPrefix				=> '',
	mlen					=> 20,
	msim					=> 0.9,
	vc 			 			=> '',
	silent				=> '',
);

GetOptions(
	"h"								=> \$option{h},
  "ref=s"          	=> \$option{ref},
  "seq=s"   		    => \$option{seq},
  "oPrefix=s"				=> \$option{oPrefix},
  "mlen=i"					=> \$option{mlen},
  "msim=f"					=> \$option{msim},
  "vc=s"            => \$option{vc},
  "silent"					=> \$option{silent},
) || printHelp (); 

if ($option{h}) { printHelp();}

my $oPrefix = $option{oPrefix};

unless ($option{ref} && $option{seq} && $oPrefix) { printHelp();}

# specify the path of muscle 
my $musclepath = "/seq/annotation/bio_tools/muscle/3.8/muscle";

my $cat_file = $oPrefix.".cat.fa";
my $aln_file = $oPrefix.".muscle.afa";

unless ($option{silent}) {
	print "\ncat $option{ref} $option{seq} > $cat_file\n\n";
}	
system ("cat $option{ref} $option{seq} > $cat_file");
unless ($option{silent}) {
	print "\n$musclepath -in $oPrefix.cat.fa -out $aln_file -quiet\n\n";
}	
system("$musclepath -in $oPrefix.cat.fa -out $aln_file -quiet");

system ("rm $cat_file");

# read in alignment file
open (ALN, "<$aln_file") or die "unable to open file $aln_file to read\n";
my $ref = "";
my $seq = "";
my $cnt = 0;
while (<ALN>){
	if (/>/) {
		$cnt ++;
		if ($cnt > 2) {
			print "ERR: more than 2 symbol of > found\n";
			exit;
		}
	} else {
		s/\s+//g; # remove whitespace
		if ($cnt == 1) {
	    $ref .= $_; # add sequence
	  } else {
	  	$seq .= $_; 
	  }
	}
}
close (ALN);

uc($ref);
uc($seq); # convert to upper case


# ---- identify alignment region ----
my $aln_start = -1;
my $aln_end = -1;
identify_aln_region($ref, $seq, $option{mlen}, $option{msim}, \$aln_start, \$aln_end);

unless ($option{silent}) {
	print "aln_start , end = $aln_start, $aln_end \n\n";
}	
#-------------------------------------------------------------------------------------------------------------------------------------------------
# translates the coordinates of $seq to $ref, calculate alignment similarity etc.
my $prefix = substr($ref, 0, $aln_start);
my $suffix = substr($ref, $aln_end + 1, length($ref) - $aln_end - 1); 
$prefix =~ s/-//g;
$suffix =~ s/-//g;
my $idx_ref = length ($prefix);
my $ref_len = length($prefix) + length($suffix);

$prefix = substr($seq, 0, $aln_start);
my $suffix = substr($seq, $aln_end + 1, length($seq) - $aln_end - 1); 
$prefix =~ s/-//g;
$suffix =~ s/-//g;
my $idx_seq = length ($prefix);
my $seq_len = length($prefix) + length($suffix);

my $num_indel = 0;
my $num_subst = 0;


# check aligned region

my %coord_map; # map coordinate of seq to ref

unless ($option{silent}) {
	print "\nAliged Region: aln_start = $aln_start, aln_end = $aln_end \n\n";
}

for (my $i = $aln_start; $i <= $aln_end; ++ $i) {
	my $ref_base = substr($ref, $i, 1);
	my $seq_base = substr($seq, $i, 1);
	
	if ($ref_base ne '-' && $seq_base ne '-') {
		
		$coord_map{$idx_seq} = $idx_ref;
		
		if ($ref_base ne $seq_base) {
			unless ($ref_base eq 'N' || $seq_base eq 'N') { # disregard "N"s
				++ $num_subst;
			}
		} 
		
		++ $idx_seq;
		++ $idx_ref;
		if ($ref_base ne 'N') { ++ $ref_len; }
		if ($seq_base ne 'N') {	++ $seq_len; }
		
	} else {
 		if ($ref_base ne '-') {
 			if ($ref_base ne 'N') {	++ $ref_len; }
			++ $idx_ref;
		} 	
		if ($seq_base ne '-') {
			if ($seq_base ne 'N') {	++ $seq_len; }
			++ $idx_seq;
		} 
		
		++ $num_indel; 	
	}
}

my $aln_len = $aln_end - $aln_start + 1;
my $similarity = 100-100*($num_indel + $num_subst)/$aln_len;
my $span = 100*$aln_len/$ref_len;
$span = sprintf "%.2f", $span;
$similarity = sprintf "%.2f", $similarity;
#print "\nRef_len\tSeq_len\tAln_len\tSpan(%)\tSubsts\tIndels\tSimilarity(%)\n";
#print "$ref_len\t$seq_len\t$aln_len\t$span\t$num_subst\t$num_indel\t$similarity\n\n";

# read in variant calling file and convert to the ref coordinates; the first column of vc file has to be the coordinates

my $num_var = 0;
if ($option{vc}) {

	unless ($option{silent}) {
		print "\nConvert Variant calling output to the ref coordinates, \n Note: the first column registers the coordinates\n\n";
	}
	my $vc_output = $oPrefix.".ref.vc.txt";
	open (OUTPUT, ">$vc_output") or die "unable to open $vc_output to write\n";
	
	open (VPH, "<$option{vc}") or die "unable to open file $option{vc} to read\n";

	while (<VPH>) {
    next if /^(\s)*$/;  # skip blank lines
		my $line = $_;
		if ($line !~ /#/) {
			++ $num_var;
			my @entries = split ("\t",$line);
			my $num = $#entries + 1;

			my $pos = $entries[0];
			-- $pos;
			if (exists $coord_map{$pos}) {
				my $map_pos = $coord_map{$pos} + 1;
				print OUTPUT $map_pos."\t";
				for (my $i = 1; $i < $num; ++ $i) {
					chomp($entries[$i]);
					print OUTPUT $entries[$i]."\t";
				}
				print OUTPUT "\n";
			} else {
				
				#print "\nErr: Coordinate $pos doesn't exists !\n";
				#exit;
			}	
		}	
	} # while (<VPH>)
	close (VPH);
	close (OUTPUT);
}
# ----------------------------------
if (!$option{silent}) {
	print "\nRef_len\tSeq_len\tAln_len\tSpan(%)\tSubsts\tIndels\tSim(%)\tnum_var\n";
	print "$ref_len\t$seq_len\t$aln_len\t$span\t$num_subst\t$num_indel\t$similarity\t$num_var\n";
} else {
	print "$span\t$similarity\t";
}

unless ($option{silent}) {
	print "\nDONE !\n\n";
}
# ---------------------------------------------------------------------------------------------------------
# sub functions
# ---------------------------------------------------------------------------------------------------------
sub printHelp {
		print "\n--------------------------------------------------------------------\n";
		print "usage: ./align2ref.pl -ref [ref.fa] -seq [seq.fa] -oPrefix [output_prefix]\n\n";
		
		print "\t-ref: reference in fasta format onto which assembly will be mapped\n";				
		print "\t-seq: sequence to map to reference\n";
		print "\t-oPrefix: output Prefix\n";
		print "\t-mlen: default 20; required to specify the valid alignment region b/t ref and seq\n";
		print "\t-msim: default 0.9; required to specify the valid alignment region similarity\n";		
		print "\t\t where no indels occur in the mlen bp of prefix and suffix region\n";
		print "\t-vc: OPT; variant calling result file by vphaser 2 with respect to the sequence\n";
		print "\n--------------------------------------------------------------------\n";

		exit;
}

# @brief -- find the start and end coordinates for a given pairwise alignment satisfying
#						the alignment shall not contain any indel for	at least min_len bases on both ends
sub identify_aln_region {
	my $seq = shift;
	my $seq2 = shift;
	my $min_len = shift;
	my $min_sim = shift;
	my $start = shift;
	my $end = shift;
	my $len = length ($seq);
	if ($len != length($seq2)) {
		print "In identify_aln_region: two input sequences should have identical length";
		exit;
	}
	
	uc ($seq);
	uc ($seq2);
	my $cnt = 0;
	my $diff = 0;
	# prefix
	for (my $i = 0; $i < $len; ++ $i) {
		my $base = substr($seq, $i, 1);
		my $base2 = substr($seq2, $i, 1); 
		
		if ($base ne '-' && $base2 ne '-') { 			
			++ $cnt; 
			if ($base ne $base2) { ++ $diff; }
			if ($cnt >= $min_len && ($diff/$cnt <= 1 - $min_sim)) {
				$$start = $i;
				last;
			}
		}	else { 
			$cnt = 0; 
			$diff = 0;
		}	
	}
	$$start -= ($min_len - 1) ;
	#suffix 
	$cnt = 0; 
	$diff = 0;
	for (my $i = $len - 1; $i >= 0; -- $i) {
		my $base = substr($seq, $i, 1);
		my $base2 = substr($seq2, $i, 1); 
		if ($base ne '-' && $base2 ne '-') { 			
			++ $cnt; 
			if ($base ne $base2) { ++ $diff; }
			if ($cnt >= $min_len && ($diff/$cnt <= 1 - $min_sim)) {
				$$end = $i;
				last;
			}
		}	else { 
			$cnt = 0; 
			$diff = 0;
		}	
	}	
	$$end += ($min_len - 1);
} # identify_aln_region

 

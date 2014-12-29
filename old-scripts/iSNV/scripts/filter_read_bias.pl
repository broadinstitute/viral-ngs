open FILE, $ARGV[0];

$POS_IDX = 0;
$VAR_IDX = 1;
$REF_IDX = 2;
$CT_IDX = 6;
$MIN_READS_FWD = 5;
$MIN_READS_REV = 5;
$STRAND_BIAS_PERC = 10;
$STRAND_BIAS_FRAC = 1.0/$STRAND_BIAS_PERC;
$RATIO_STRAND_BIAS_PERC = 5;
$RATIO_STRAND_BIAS_FRAC = 1.0/$RATIO_STRAND_BIAS_PERC;
while ($line = <FILE>) {
    chomp $line;
    @parts = split /\t/, $line;
      for ($i = $CT_IDX; $i <= $#parts; $i++) {
	#print "$parts[$i]\n";
	@read_cts = split /:/, $parts[$i];
	#print "$parts2[0]\n";
	#if we are looking at the number of reads for the variant
	if ($read_cts[0] eq $parts[$REF_IDX]) {
	    #print "REF:$read_cts[0]\t$parts[$REF_IDX]\t$read_cts[1]\t$read_cts[2]\n";
	    $ref_fwd = $read_cts[1];
	    $ref_rev = $read_cts[2];
	}

	if ($read_cts[0] eq $parts[$VAR_IDX]) {
	    $var_fwd = $read_cts[1];
	    $var_rev = $read_cts[2];
	}
    }
    
    #print "$parts[$VAR_IDX]\t$parts[$REF_IDX]\t$ref_fwd\t$ref_rev\t$var_fwd\t$var_rev\n";

    #pseudocount to get rid of zeroes
    $ref_strand_bias = (1.0 * $ref_fwd + 1) / (1.0 * $ref_rev + 1);
    $var_strand_bias = (1.0 * $var_fwd + 1) / (1.0 * $var_rev + 1);
    $ratio = $ref_strand_bias / $var_strand_bias;
   
    #print "$ratio\t$STRAND_BIAS_PERC\t$STRAND_BIAS_PERC_2\t$var_fwd\t$var_rev\t$MIN_READS\n";
    if ($ratio < $RATIO_STRAND_BIAS_PERC && $ratio > $RATIO_STRAND_BIAS_FRAC && $var_strand_bias < $STRAND_BIAS_PERC && $var_strand_bias > $STRAND_BIAS_FRAC  && $var_fwd > $MIN_READS_FWD && $var_rev > $MIN_READS_REV) {
	print "$line\n";
    }
}


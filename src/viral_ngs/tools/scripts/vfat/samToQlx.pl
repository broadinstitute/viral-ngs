#!/usr/bin/env perl

# Copyright © 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.
#
# This code was developed by Patrick Charlebois <patrickc@broadinstitute.org>

use strict;
use Getopt::Long;

my %option = (
	nqsmq     	=> 20,
	nqsaq     	=> 15,
	nqssize   	=> 5,
	nqvalue		=> 1,
	nosuffix	=> '',
	fakequals	=> 0,
	h				=> '',
);


GetOptions(
  "nqsmq=i"		=> \$option{nqsmq},
  "nqsaq=i"		=> \$option{nqsaq},
  "nqssize=i"		=> \$option{nqssize},
  "nqvalue=i"		=> \$option{nqvalue},
  "nosuffix"		=> \$option{nosuffix},
  "fakequals=i"	=> \$option{fakequals},
  "h"					=> \$option{h},
) || die("Problem processing command-line options: $!\n");

my $saminput = shift  || printHelp();
my $ref = shift || printHelp();
my $output = shift || printHelp();

if($option{h}){printHelp();}

open(SAMFILE, $saminput) || die("Unable to open $saminput\n");
open(REF, $ref) || die("Unable to open $ref\n");
open(OUTPUT, ">$output.qlx") || die;
open(UNMAPPED, ">$output"."_unmappedIDs.txt") || die;

my %refseqs;
my $currefname = '';
my $currefseq = '';
while(my $line = <REF>)
{
	chomp $line;
	if($line =~ />(.+)/)
	{
		if($currefseq)
		{
			$refseqs{$currefname} = $currefseq;
		}
		$currefname = $1;
		next;
	}else{
		$currefseq .= $line;
	}
}
$refseqs{$currefname} = $currefseq;

my $curreadid = '';
my %readseqs;
my %readquals;
my %nqsquals;

my @flagkeys; ### Used for conversion of the FLAG integer in column 2 of SAM format
$flagkeys[0] = 'multsegments';
$flagkeys[1] = 'propaligned';
$flagkeys[2] = 'unmapped';
$flagkeys[3] = 'nextunmapped';
$flagkeys[4] = 'revcomp';
$flagkeys[5] = 'nextrevcomp';
$flagkeys[6] = 'firstseg';
$flagkeys[7] = 'lastseg';
$flagkeys[8] = 'secondaryalign';
$flagkeys[9] = 'badqual';
$flagkeys[10] = 'duplicate';

my $readid;
my $readstart;
my $readend;
my $readlen;
my $strand;
my $refstart;
my $refend;
my $cigarstr = '';
my $refstr = '';
my $readstr = '';
my $qualstr = '';
my $mode;
my $rawreadstr;
my $rawqualstr;
my %unaligned;
my %flags;

while (my $line = <SAMFILE>)
{
	chomp $line;
	if($line =~ /^\@.*/)
	{
		next;
	}
		
	if(length($line) > 1)
	{
		my @linedata = split(/\t/,$line);
		$readid = $linedata[0];
		if($readid =~ /(.+) (.*)/)
		{
			$readid = $1;
		}
		
		undef %flags;
		assignFlags($linedata[1]);
		if($flags{multsegments} == 1 && $flags{firstseg} == 1 && !($option{nosuffix}))
		{
			$readid .= "/1";
		}elsif($flags{multsegments} == 1 && $flags{lastseg} == 1 && !($option{nosuffix}))
		{
			$readid .= "/2";
		}

		if($flags{unmapped})
		{
			$unaligned{$readid} = 1;
			print UNMAPPED "$readid\n";
			next;
		}
		
		$ref = $linedata[2];
		$refstart = $linedata[3];
		if($flags{revcomp} == 1){$strand = '-'}else{$strand = '+'};
		$cigarstr = $linedata[5];
		$rawreadstr = $linedata[9];
		$rawqualstr = $linedata[10];
		
		$readstr = '';
		$refstr = '';
		$qualstr = '';
		$refend = $refstart - 1;
		$readlen = 0;
		$readstart = 1;

		my $refpos = $refstart - 1;
		my $readpos = 0;
		my $readextra = 0;
		while(length($cigarstr) > 0)
		{
			$cigarstr =~ /([0-9]+?)([A-Z]+?)(.*)/;
			my $len = $1;
			my $letter = $2;
			$cigarstr = $3;
		
			if($letter eq 'M' || $letter eq '=' || $letter eq 'X')
			{
				$readlen += $len;
				$refend += $len;
				$readstr .= substr($rawreadstr, $readpos, $len);
				$qualstr .= substr($rawqualstr, $readpos, $len);
				$refstr .= substr($refseqs{$ref}, $refpos, $len);
				$readpos += $len;
				$refpos += $len;
			}elsif($letter eq 'D')
			{
				$refend += $len;
				$readstr .= "-" x $len;
				$qualstr .= " " x $len;
				$refstr .= substr($refseqs{$ref}, $refpos, $len);
				$refpos += $len;
			}elsif($letter eq 'I')
			{
				$readlen += $len;
				$refstr .= "-" x $len;
				$readstr .= substr($rawreadstr, $readpos, $len);
				$qualstr .= substr($rawqualstr, $readpos, $len);
				$readpos += $len;
			}elsif($letter eq 'S')
			{
				$readlen += $len;
				if($readlen == $len)
				{
					$readstart += $len;
					$readpos += $len;
				}else{
					$readextra = $len;
				}
			}elsif($letter eq 'H')
			{
			        $readlen += $len;
				if($readlen == $len)
				{
				        $readstart += $len;
					# we don't advance $readpos here
					# because the rawreadstr is clipped
				}else{
				        $readextra = $len;
				}
			}
		}
		
		my $nqsstr;
		if($option{fakequals})
		{
			my $newqualstr = '';
			for(my $pos = 0; $pos < length($qualstr); $pos++)
			{
				if(substr($qualstr, $pos, 1) ne ' ')
				{
					$newqualstr .= chr($option{fakequals} + 33);
					$nqsstr .= '1';
				}else{
					$newqualstr .= ' ';
					$nqsstr .= '8';
				}
			}
			$qualstr = $newqualstr;
		}else{
			$nqsstr = generateNQS($qualstr);
		}
		$readend = $readlen - $readextra;

		# handle reads that have clipping and are reverse oriented

		if ($strand eq '-') {
		  $readend = $readlen - $readstart + 1;
		  $readstart = $readextra + 1;
		}
		
		print OUTPUT ">Read $readid $readstart $readend $readlen $strand $ref $refstart $refend\n";
		print OUTPUT $refstr."\n".$readstr."\n".$nqsstr."\n".$qualstr."\n\n";
	}
}

sub generateNQS
{
	my $qualstr = shift;
	my $nqsstr = '';
	my $next5 = '';
	my $prev5 = '';
	my $prevqual = 0;
	
	for(my $qualpos = 0; $qualpos < length($qualstr); $qualpos++)
	{	
		if(substr($qualstr, $qualpos, 1) eq ' '){
			if($prevqual == 0){$nqsstr .= '9';}else{$nqsstr .= '8';}
			next;
		}
		my $curqual = ord(substr($qualstr, $qualpos, 1)) - 33;
		
		$next5 = '';
		
		my $qualpos2 = $qualpos;
		while(length($next5) < $option{nqssize})
		{
			$qualpos2++;
			if($qualpos2 >= length($qualstr)){last;}
			if(substr($qualstr, $qualpos2, 1) eq ' '){next;}
			if((ord(substr($qualstr, $qualpos2, 1)) - 33) == $option{nqvalue}){next;}
			
			if((ord(substr($qualstr, $qualpos2, 1)) - 33) < $option{nqsaq})
			{
				$next5 .= "0";
			}else{
				$next5 .= "1";
			}
		}
		if($curqual < $option{nqsmq})
		{
			$nqsstr .= "0";
			$prevqual = 0;
		}elsif($next5 =~ /0/ || $prev5 =~ /0/)
		{
			$nqsstr .= "0";
			$prevqual = 0;
		}else{
			$nqsstr .= "1";
			$prevqual = 1;
		}
		if($curqual == $option{nqvalue}){next;}
		if($curqual >= $option{nqsaq}){$prev5 .= "1";}else{$prev5 .= "0";}
		if(length($prev5) > $option{nqssize}){$prev5 = substr($prev5, 1);}
	}
	return $nqsstr;
}

sub assignFlags
{
	my $intflag = shift;
	my $binno = sprintf("%b",$intflag);
		
	for(my $i = 0; $i <= 10; $i++)
	{
		if(length($binno) > $i && substr($binno, length($binno) - ($i+1), 1) == 1)
		{
			$flags{$flagkeys[$i]} = 1;
		}else{
			$flags{$flagkeys[$i]} = 0;
		}
	}
}

sub printHelp
{
	print "Usage : perl samToQlx.pl input.sam refseq.fa outputBaseName\n";
	print "Options:\n";
	print "-nqsmq : NQS (Neighborhood Quality Score) main base minimum quality (20 default)\n";
	print "-nqsaq : NQS adjacent bases minimum quality (10 default)\n";
	print "-nqssize : NQS window size for adjacent bases (5 default)\n";
	print "-nqvalue : Quality value for inserted Ns that should not influence adjacent bases NQS score (1 default)\n";
	print "-fakequals : fakes qualities in the qlx output to a given amount (off by default)\n";
	print "-nosuffix : By default, in paired reads data the read name will have a suffix added of /1 or /2 depending on the .sam file. This option removes this addition\n";
	exit();
}

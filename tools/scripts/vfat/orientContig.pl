#!/usr/bin/env perl

# Copyright © 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

use strict;
use Getopt::Long;
my %option = (
	minlongcont		=> 100, # Minimum length of a raw contig. Any contig below this length will be filtered out from the start
	mincontlen		=> 350,
	maxorigaplen	=> 10,
	h				=> '',
);

GetOptions(
	"minlongcont=i"		=> \$option{minlongcont},
	"mincontlen=i"		=> \$option{mincontlen},
	"maxorigaplen=i"	=> \$option{maxorigaplen},
	"h"					=> \$option{h},
);

if($option{h})
{
	print "Usage : perl orientContig.pl inputContigs.fa reference.fa outputBaseName\n";
	print "orientContig.pl requires MUSCLE to run\n";
	print "Options:\n";
	print "-mincontlen(350)\tMinimum length of contigs that will be analyzed\n";
	print "-minlongcont(100)\tMinimum length for the longest continuous stretch of the contig aligning to reference\n";
	exit();
}

my $inputcontIn = shift || die("Usage : perl orientContig.pl inputContigs.fa reference.fa outputBaseName\nFor more details about the options use -h\n");
my $refIn = shift || die("Usage : perl orientContig.pl inputContigs.fa reference.fa outputBaseName\nFor more details about the options use -h\n");
my $outputIn = shift || die("Usage : perl orientContig.pl inputContigs.fa reference.fa outputBaseName\nFor more details about the options use -h\n");

my $musclepath = "/seq/annotation/bio_tools/muscle/3.8/";

my @contigIDs;

orientContig($inputcontIn, $refIn, $outputIn."_orientedContigs");

my $reffastaname = "";

sub orientContig
{
	my $inputcont = shift;
	my $reffile = shift;
	my $output = shift;
	
	my %contList;

	open(REFFILE, $reffile) || die("Unable to open $reffile\n");
	open(OUTPUT, ">$output.fa") || die("invalid output $output");
	
	my $refseq = '';
	while(my $line = <REFFILE>)
	{
		chomp $line;
		if($line =~ />(.+)/)
		{
			$reffastaname = $1;
		}else{
			$refseq .= $line;
		}
	}

	open(INPUTCONT, $inputcont) || die("invalid contig file\n");
	my $contId = '';
	my $contSeq = '';
	my $contflag = 0;
	while(my $line = <INPUTCONT>)
	{
		chomp $line;
		if($line =~ />(.+)/)
		{
			if($contSeq)
			{
				if(length($contSeq) < $option{mincontlen}){
					print "Rejected $contId (length ".length($contSeq)." < ".$option{mincontlen}.")\n";
				}else{
					$contList{$contId} = $contSeq;
				}
				$contSeq = '';
			}
			
			$contId = 'contig_'.$contflag;
			$contigIDs[$contflag] = $1;
			$contflag++;
		}else{
			$contSeq .= $line;
		}
	}
	if(length($contSeq) < $option{mincontlen}){
		print "Rejected $contId (length ".length($contSeq)." < ".$option{mincontlen}.")\n";
	}else{
		$contList{$contId} = $contSeq;
	}


	foreach my $curCont (sort keys %contList)
	{
		open(CURCONTOUT, ">$output"."_$curCont.orimfa");
		print CURCONTOUT ">$curCont\n".$contList{$curCont}."\n>$reffastaname\n$refseq\n";
		close CURCONTOUT;

		system($musclepath."muscle -in $output"."_$curCont.orimfa -out $output"."_$curCont.oriafa -quiet");
		my $inputlen = length($contList{$curCont});
	
		my $longestcont = 0;
	
		(my $initpctgap, my $initasspctgap, $longestcont) = calcPctGap($output."_$curCont.oriafa",$longestcont);
		my $longestfowcont = $longestcont;
		$longestcont = 0;
		
		my $curRevSeq = reverseSeq($contList{$curCont});
		
		open(CURREVOUT, ">$output"."_$curCont.orirevmfa");
		print CURREVOUT ">$curCont\n$curRevSeq\n>$reffastaname\n$refseq\n";
		close CURREVOUT;
		
		system($musclepath."muscle -in $output"."_$curCont.orirevmfa -out $output"."_$curCont.orirevafa -quiet");
	
		(my $revpctgap, my $revasspctgap, $longestcont) = calcPctGap($output."_$curCont.orirevafa",$longestcont);
		my $longestrevcont = $longestcont;
		if($longestfowcont > $longestcont){$longestcont = $longestfowcont;}
	
		my $pctlongcont = ($longestcont/$inputlen);
	
		my $goodrevori = 0;
		my $ambiguous = 0;
		my $veryambiguous = 0;
	
		if($longestfowcont > ($longestrevcont * 2) && $pctlongcont > 0.5){
			$goodrevori = 0;
		}elsif($longestrevcont > ($longestfowcont * 2) && $pctlongcont > 0.5)
		{
			$goodrevori = 1;
		}elsif($initasspctgap < ($revasspctgap / 2))
		{
			$goodrevori = 0;
		}elsif($revasspctgap < ($initasspctgap / 2))
		{
			$goodrevori = 1;
		}elsif($initasspctgap < $revasspctgap && $longestfowcont > $longestrevcont)
		{
			$goodrevori = 0;
			$ambiguous = 1;
		}elsif($initasspctgap > $revasspctgap && $longestfowcont < $longestrevcont)
		{
			$goodrevori = 1;
			$ambiguous = 1;
		}elsif($longestfowcont < $longestrevcont)
		{
			$goodrevori = 1;
			$veryambiguous = 1;
		}else{
			$goodrevori = 0;
			$veryambiguous = 1;
		}
		
		if($longestfowcont < $option{minlongcont} && $longestrevcont < $option{minlongcont})
		{
			print "$curCont : Bad contig, both orientation can't give a continuous stretch of greater length than ".$option{minlongcont}."\n";
			system("rm $output"."_$curCont.orimfa");
			system("rm $output"."_$curCont.orirevmfa");
			system("rm $output"."_$curCont.orirevafa");
			system("rm $output"."_$curCont.oriafa");
			next;
		}
		
		
		if($goodrevori){
			print "$curCont : Reverse ($revasspctgap vs $initasspctgap % contig gaps, longest contig $longestrevcont vs $longestfowcont)\n";
			print OUTPUT ">$curCont\n$curRevSeq\n";
			
			system("rm $output"."_$curCont.orimfa");
			system("rm $output"."_$curCont.oriafa");
			system("rm $output"."_$curCont.orirevmfa");
			system("mv $output"."_$curCont.orirevafa $output"."_$curCont"."_orientalign.afa");
		}else{
			print "$curCont : Forward ($initasspctgap vs $revasspctgap % contig gaps, longest contig $longestfowcont vs $longestrevcont)\n";
			print OUTPUT ">$curCont\n".$contList{$curCont}."\n";
			system("rm $output"."_$curCont.orimfa");
			system("rm $output"."_$curCont.orirevmfa");
			system("rm $output"."_$curCont.orirevafa");
			system("mv $output"."_$curCont.oriafa $output"."_$curCont"."_orientalign.afa");
		}
	}
}


sub calcPctGap
{
	my $alignfilename = shift;
	my $longestcont = shift;

	my $maxgaplen = $option{maxorigaplen};

	open(AFA, "$alignfilename");

	my $refseq = "";
	my $otherseq = "";
	my $curseq = "";
	
	while (my $line = <AFA>)
	{
		chomp $line;
		if($line =~ />(.+)/)
		{
			if($curseq)
			{
				if($1 eq $reffastaname)
				{
					$otherseq = $curseq;
				}else{
					$refseq = $curseq;
				}
				$curseq = "";
			}
		}else
		{
			$curseq .= $line;
		}
	}
	
	if($refseq)
	{
		$otherseq = $curseq;
	}else{
		$refseq = $curseq;
	}
	
	my $refpos = 0;
	my $otherpos = 0;
	my $nbrefgap = 0;
	my $nbassgap = 0;
	my $startedref = 0;
	my $startedass = 0;
	my $currefgaplen = 0;
	my $curassgaplen = 0;
	
	my $curcontlen = 0;
	
	for(my $seqpos = 0; $seqpos < length($refseq); $seqpos++)
	{
		my $curref = substr($refseq, $seqpos, 1);
		my $curother = substr($otherseq, $seqpos, 1);
		
		if($curref eq "-")
		{
			$startedass = 1;
			$currefgaplen++;
			if($startedref){$nbrefgap++;}
			$otherpos++;
			$curassgaplen = 0;
			$curcontlen++;
		}else{
			$startedref = 1;
			if($curother ne "-")
			{
				$curcontlen++;
				$startedass = 1;
				$otherpos++;
				$refpos++;
				$curassgaplen = 0;
			}else{
				if($startedass){$nbassgap++;}
				$curassgaplen++;
				if($curassgaplen > $maxgaplen){
					#print $curassgaplen."\n";
					if($curcontlen > $longestcont)
					{	
						$longestcont = $curcontlen;
#						print "Long at $seqpos ".$longestcont."\n";
					}
					$curcontlen = 0;
				}
				$refpos++;
			}
			$currefgaplen = 0;
		}
	}
	if($curcontlen > $longestcont)
	{
		$longestcont = $curcontlen;
#		print "Long at ".length($refseq)." $longestcont\n";
	}
	$nbrefgap -= $currefgaplen;
	$nbassgap -= $curassgaplen;
	my $pctrefgap = $nbrefgap/length($refseq);
	my $pctassgap = $nbassgap/length($otherseq);
	return ($pctrefgap, $pctassgap, $longestcont);	
}

sub reverseSeq
{
	my $seq = shift || die;

	$seq =~ tr/ACGTacgt/TGCAtgca/;
	my $rev = reverse $seq;
	return $rev;
}
			

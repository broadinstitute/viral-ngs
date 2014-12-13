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
use warnings;

my %option = (
	mincontlen 	=> 350,	# Minimum length of contigs that will be analyzed
	readfa				=> '',		# Input file containing the reads in fasta format
	readq				=> '',		# Input file containing the read quals
	readfa2				=> '',		# Input file containing the reads in fasta format
	readq2				=> '',		# Input file containing the read quals
	readfq				=> '',		# Input file containing the reads in fasta format
	readfq2				=> '',		# Input file containing the read quals
	maxseggap			=> 30,
	maxsegins			=> 60,
	minseglen			=> 50,
	maxsegdel			=> 0.25,
	longcont			=> 1500,
	contigcov			=> 2,
	allcont				=> '',
	sequencer			=> 'illumina',
	fakequals			=> 30,
	h					=> '',
);

GetOptions(
	"mincontlen=i"	=> \$option{mincontlen},
	"readfa=s"			=> \$option{readfa},
	"readq=s"			=> \$option{readq},
	"readfa2=s"			=> \$option{readfa2},
	"readq2=s"			=> \$option{readq2},
	"readfq=s"			=> \$option{readfq},
	"readfq2=s"			=> \$option{readfq2},
	"maxseggap=i"		=> \$option{maxseggap},
	"maxsegins=i"		=> \$option{maxsegins},
	"minseglen=i"		=> \$option{minseglen},
	"maxsegdel=f"		=> \$option{maxsegdel},
	"longcont=i"		=> \$option{longcont},
	"contigcov=i"		=> \$option{contigcov},
	"allcont"			=> \$option{allcont},
	"sequencer=s"		=> \$option{sequencer},
  	"fakequals=i"		=> \$option{fakequals},
	"h"					=> \$option{h},
) || die("Problem processing command-line options: $!\n");

if($option{h})
{
	print "Usage : perl contigMerger.pl orientContigsOutputName reference.fa outputBaseName\n";
	print "contigMerger.pl requires	MUSCLE to run, and Mosaik if reads are used\n";
	print "Options:\n";
	print "-readfa\t\tReads file in fasta format\n";
	print "-readq\t\tReads quality file in fasta format (spaced integers)\n";
	print "-readfa2\t\tReads file in fasta format (2nd mate if paired)\n";
	print "-readq2\t\tReads quality file in fasta format (spaced integers) (2nd mate if paired)\n";
	print "-readfq\t\tReads file in fastq format\n";
	print "-readfq2\t\tReads file in fastq format (2nd mate if paired)\n";
	print "-sequencer(illumina)\tSets which sequencer was used to determine which version of Mosaik to run. Currently supports illumina, 454\n";
	print "-fakequals()\tFake all quality scores to a given value in qlx files. This won't affect the results of any script in this package and will speed it up, but qlx files will not be valid for other softwares using them like V-Phaser\n";
	print "-maxseggap(30)\tMaximum gap length before splitting segments\n";
	print "-maxsegins(60)\tMaximum insertion length before splitting segments\n";
	print "-minseglen(50)\tMinimum length of a segment\n";
	print "-maxsegdel(0.25)\tMax % of deletion allowed in a segment\n";
	print "-longcont(1500)\tMin length to be a 'long contig', which will always be considered\n";
	print "-contigcov(2)\tContigs will be considered as long as you have a coverage below contigcov\n";
	print "-allcont\tAnalyze all contigs, regardless of other filters\n";
	exit();
}

my $orientContOut = shift || die ("Usage : perl contigMerger.pl orientContigsOutputName reference.fa outputBaseName\nFor option details and how to incorporate read data use the option -h\n");
my $reffile = shift || die ("Usage : perl contigMerger.pl orientContigsOutputName reference.fa outputBaseName\nFor option details and how to incorporate read data use the option -h\n");
my $output = shift || die ("Usage : perl contigMerger.pl orientContigsOutputName reference.fa outputBaseName\nFor option details and how to incorporate read data use the option -h\n");

# We definitely need the MUSCLE aligner
my $musclepath = "/seq/annotation/bio_tools/muscle/3.8/";
# This is to call runMosaik2.pl, which is necessary
my $scriptpath = "/idi/sabeti-scratch/kandersen/bin/VfatSoftwarePackage/";
# We can ignore R (and drawSegments()) if we don't care about the graphical output
my $Rpath = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_2.15.1/bin/";

my $hasreads = '';

if($option{readfa} || $option{readfq}){$hasreads = 1;}

my $mosaikParam = " -hgop 20 -gop 40 -gep 10 -bw 29 -st illumina";
if($option{sequencer} eq '454')
{
	$mosaikParam = " -hgop 4 -gop 15 -gep 6.66 -bw 0 -st 454";
}
if($option{fakequals})
{
	$mosaikParam .= " -fakequals ".$option{fakequals};
}

my $reffastaname = '';
(my $refid, my $refseq) = readFasta($reffile);

my $maxgaplen = $option{maxseggap};
my $mincontlen = $option{minseglen};
my $nbCont = 0;

my  %contigSeqs;
my %contigAlign;
readOrientContig($orientContOut.".fa");

my %segments;
my %segOverlap;
my %segInserts;
my @segMap;

# Build segments, draw map, filter contigs
open(INDELS, ">$output"."_largeDeletions.txt");
print INDELS "ContigNo\tContigStart\tContigStop\tReferenceStart\tReferenceStop\n";
buildSegments();
drawSegments();
#die;
buildOverlapMap();

open(OUTPUT, ">$output"."_assembly.fa");
open(OUTALN, ">$output"."_vsRef.mfa");
open(CONTIGLIST, ">$output"."_mergingList.txt");

my $consensusSeq;
my %ntfreq;
if($hasreads)
{
	alignReads();
	calcNtFreq();
	$consensusSeq = mergeContigsReads();
}else{
	$consensusSeq = mergeContigsNoRead();
}

print OUTPUT ">Consensus\n$consensusSeq\n";
print OUTALN ">Consensus\n$consensusSeq\n>$refid\n$refseq\n";
close OUTALN;
close OUTPUT;
system($musclepath."muscle -in $output"."_vsRef.mfa -out $output"."_vsRef.afa -quiet");

if($option{readfa} && $option{readq})
{
	system("perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -ref $output"."_assembly.fa -o $output"."_readAlignment -qlx -qlxonly".$mosaikParam);
}

###### SUBS ######

sub calcNtFreq
{
	foreach my $contid (keys %segments)
	{
		fillNtFreq("$output"."_$contid"."_aligned.qlx", $contid);
	}
}

sub alignReads
{
	foreach my $contid (keys %segments)
	{
		my $file = "$output"."_$contid.fa";
		open(TMPFA, ">$file");
		print TMPFA ">$contid\n".$contigSeqs{$contid}{seq}."\n";
		close TMPFA;
		if($option{readfa} && $option{readq})
		{
			if($option{readfa2} && $option{readq2})
			{
				system("perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -fa2 ".$option{readfa2}." -qual2 ".$option{readq2}." -ref $file -o $output"."_$contid"."_aligned -qlx -qlxonly".$mosaikParam);
			}else{
				system("perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -ref $file -o $output"."_$contid"."_aligned -qlx -qlxonly".$mosaikParam);
			}
		}elsif($option{readfq})
		{
			if($option{readfq2})
			{
				system("perl $scriptpath"."runMosaik2.pl -fq ".$option{readfq}." -fq2 ".$option{readfq2}." -ref $file -o $output"."_$contid"."_aligned -qlx -qlxonly".$mosaikParam);
			}else{
				system("perl $scriptpath"."runMosaik2.pl -fq ".$option{readfq}." -ref $file -o $output"."_$contid"."_aligned -qlx -qlxonly".$mosaikParam);
			}
		}
	}
}

sub mergeContigsReads
{
	my $consseq = '';
	my @contseq;
	my $prevstop = 0;
	my $prevstart = 0;
	my %allUsedCont;
	foreach my $curseg(@segMap)
	{
		#print $$curseg{refstart}."\t".$$curseg{refstop}."\n";
		if($$curseg{refstart} > $prevstop + 1 && $prevstop > 0)
		{
			$consseq .= 'N' x ($$curseg{refstart} - $prevstop - 1);
			push(@contseq, ($prevstop + 1)."\t".($$curseg{refstart} - 1)."\t".($$curseg{refstart} - $prevstop - 1)."N");
		} 
		$prevstop = $$curseg{refstop};
		if($$curseg{nbcontigs} == 1)
		{
			foreach my $curcont (sort {$$curseg{contigs}{$b}{length} <=> $$curseg{contigs}{$a}{length}} keys %{$$curseg{contigs}})
			{
				my $contstart = $$curseg{contigs}{$curcont}{start};
				my $contstop = $$curseg{contigs}{$curcont}{stop};

				if($$curseg{refstart} == 1){
					$contstart = 1;
				}
				if($$curseg{refstop} == length($refseq))
				{
					$contstop = length($contigSeqs{$curcont}{seq});
				}
	
	#			print ">$curcont\n".$contigSeqs{$curcont}{seq}."\n";
				my $curcontseq = substr($contigSeqs{$curcont}{seq}, ($contstart - 1), ($contstop - $contstart + 1));
	
				$consseq .= $curcontseq;
				push(@contseq, $$curseg{refstart}."\t".$$curseg{refstop}."\t".$curcont."(".$contstart."-".$contstop.")");
				$allUsedCont{$curcont} = 1;
				last;
			}
		}else{
			my $curfilename = "$output"."_segmap_".$$curseg{refstart}."-".$$curseg{refstop};
			my %usedContigs;
			undef(%usedContigs);
			open(TMPOUT, ">$curfilename.mfa");
			foreach my $curcont (sort {$$curseg{contigs}{$b}{length} <=> $$curseg{contigs}{$a}{length}} keys %{$$curseg{contigs}})
			{
				print TMPOUT ">$curcont\n".$$curseg{contigs}{$curcont}{seq}."\n";
			}
			close TMPOUT;
			system($musclepath."muscle -in $curfilename.mfa -out $curfilename.afa -quiet");
			
			my %curpos;
			my %alnseq;
			foreach my $curcont (keys %{$$curseg{contigs}})
			{
				$curpos{$curcont} = $$curseg{contigs}{$curcont}{start};
				#print $curcont."\t".$$curseg{contigs}{$curcont}{start}."\t".$$curseg{contigs}{$curcont}{stop}."\n";
			}
			
			open(TMPALN, "$curfilename.afa");
			my $curid = '';
			my $curseq = '';
			while(my $line = <TMPALN>)
			{
				chomp $line;
				if($line =~ />(.+)/)
				{
					$curid = $1;
				}else{
					$alnseq{$curid} .= $line;
				}
			}
			close TMPALN;
			my $lenaln = length($alnseq{$curid});
			my $cursegseq = '';
			my $lastcont;
			my %prevbase;
			for(my $alnpos = 0; $alnpos < $lenaln; $alnpos++)
			{
				my %curres;
				my $prevres = '';
				my $diff = 0;
				my $hasgap = 0;
				my %gapConts;
				my %nogapConts;
				my %pregap;
				foreach my $curcont (keys %{$$curseg{contigs}})
				{
					$curres{$curcont} = substr($alnseq{$curcont}, $alnpos, 1);
					if($curres{$curcont} ne '-'){
						$curpos{$curcont}++;
						$nogapConts{$curcont} = 1;
					}else{
						if($prevbase{$curcont} ne '-')
						{
							$hasgap = 1;
							$gapConts{$curcont} = 1;
							$pregap{$curcont} = $prevbase{$curcont};
						}
					}
					if($prevres){
						if($curres{$curcont} ne $prevres){
							$diff = 1;
						}
					}else{$prevres = $curres{$curcont};}
					$prevbase{$curcont} = $curres{$curcont};
					$lastcont = $curcont;
				}
				if($diff)
				{
					if($hasgap)
					{
						my $dominantcontig = '';
						my $dominantcount = 0;
						my $dominantseq = '';
						my %sortcontgap;
						my $gaplen = 0;
						foreach my $contgap (keys %gapConts)
						{
							my $tmppos = $alnpos;
							$gaplen = 0;
							while(substr($alnseq{$contgap}, $tmppos, 1) eq '-')
							{
								$gaplen++;
								$tmppos++;
							}
							#print $contgap."\t".$gaplen."\t".$pregap{$contgap}."\t".$ntfreq{$contgap}{$curpos{$contgap} - 1}{$pregap{$contgap}}{count}."\n";
							
							my $origcontgap = $ntfreq{$contgap}{$curpos{$contgap} - 1}{$pregap{$contgap}}{count};
							
							foreach my $nogapcont (keys %nogapConts)
							{
								$sortcontgap{$contgap}{count} = $origcontgap;
								$sortcontgap{$contgap}{seq} = '';
								my $insert = substr($alnseq{$nogapcont}, $alnpos, $gaplen);
								$insert =~ s/\-//g;
								#print "INSERT\t$insert\n";
								my $countinsert = 0;
								if($ntfreq{$contgap}{$curpos{$contgap} - 2}{'I'}{$insert})
								{
									$countinsert = $ntfreq{$contgap}{$curpos{$contgap} - 2}{'I'}{$insert};
								}
								my $avgcov = 0;
								for(my $tmppos = $curpos{$nogapcont}; $tmppos < $curpos{$nogapcont} + length($insert); $tmppos++)
								{
									$avgcov += $ntfreq{$nogapcont}{$tmppos - 1}{substr($insert, $tmppos - $curpos{$nogapcont}, 1)}{count};
									#print "".($tmppos - 1)."\t".substr($insert, $tmppos - $curpos{$nogapcont}, 1)."\t".$ntfreq{$nogapcont}{$tmppos - 1}{substr($insert, $tmppos - $curpos{$nogapcont}, 1)}{count}."\t".$ntfreq{$nogapcont}{$tmppos - 1}{coverage}{count}."\n";
								}
								$avgcov = $avgcov/length($insert);
								
								#print $insert."\t".$avgcov."\n";
								my $fulldelcov = 0;
								# Adjust count of reads with insertion to remove deleted reads
							#	print "$nogapcont\t".($curpos{$nogapcont} - 1)."\t".'D'.length($insert)."\n";
								if($ntfreq{$nogapcont}{$curpos{$nogapcont} - 1}{'D'.length($insert)}{count})
								{
									 $fulldelcov = $ntfreq{$nogapcont}{$curpos{$nogapcont} - 1}{'D'.length($insert)}{count};
								}
							#	print "Del cov : $fulldelcov\n";
								$sortcontgap{$nogapcont}{count} = $avgcov - $fulldelcov;
							#	print "Count : ".$sortcontgap{$nogapcont}{count}."\n";
								$sortcontgap{$nogapcont}{seq} = $insert;
								
								# Adjust count of reads with deletion to remove inserted reads
								if($ntfreq{$contgap}{$curpos{$contgap} - 1}{I}{$insert})
								{
									$sortcontgap{$contgap}{count} -= $ntfreq{$contgap}{$curpos{$contgap} - 1}{I}{$insert};
								}
							}
							
							foreach my $sortedcont (sort {$sortcontgap{$b}{count} <=> $sortcontgap{$a}{count}} keys %sortcontgap)
							{
								if($sortcontgap{$sortedcont}{count} > $dominantcount)
								{
									$dominantcontig = $sortedcont;
									$dominantcount = $sortcontgap{$sortedcont}{count};
									$dominantseq = $sortcontgap{$sortedcont}{seq};
								}
								last;
							}
							
							#print $dominantcontig."\t".$dominantcount."\t".$dominantseq."\n";
						}
						$cursegseq .= $dominantseq;
						$usedContigs{$dominantcontig} = 1;
						#print "LOOP\n";
						foreach my $curcont (keys %{$$curseg{contigs}})
						{
							if($gaplen > 0 && $curres{$curcont} ne '-'){
								$curpos{$curcont}--;
							}
						}
						for(my $tmppos = $alnpos; $tmppos < $alnpos + $gaplen; $tmppos++)
						{
							foreach my $curcont (keys %{$$curseg{contigs}})
							{
								$curres{$curcont} = substr($alnseq{$curcont}, $tmppos, 1);
								#print $curcont."\t".$curpos{$curcont}."\t".$curres{$curcont}."\n";
								if($curres{$curcont} ne '-'){
									$curpos{$curcont}++;
								}
							}
						}
						$alnpos += $gaplen - 1;
					}else{
						my %contcov;
						foreach my $curcont (keys %{$$curseg{contigs}})
						{
							$contcov{$curcont} = $ntfreq{$curcont}{$curpos{$curcont} - 1}{$curres{$curcont}}{count};
							unless($contcov{$curcont}){$contcov{$curcont} = 0;}
							#print $curcont."\t".($curpos{$curcont}-1)."\t".$curres{$curcont}."\t".$contcov{$curcont}."\n";
						}
						foreach my $curcont (sort {$contcov{$b} <=> $contcov{$a}} keys %contcov)
						{
							unless($curres{$curcont} eq '-'){
								$cursegseq .= $curres{$curcont};
								$usedContigs{$curcont} = 1;
							}
							last;
						}
					}
				}else{
					$cursegseq .= $curres{$lastcont};
				}
			}
			#print "CONS\n".$cursegseq."\n\n";
			$consseq .= $cursegseq;
			my $contstr = '';
			if(scalar(keys %usedContigs) > 0)
			{
				foreach my $usedcont (sort keys %usedContigs)
				{
					if($contstr){
						$contstr .= '-';
					}elsif(scalar(keys %usedContigs) > 1)
					{
						$contstr .= "MIX OF ";
					}
					$allUsedCont{$usedcont} = 1;
					$contstr .= $usedcont."(".$$curseg{contigs}{$usedcont}{start}."-".$$curseg{contigs}{$usedcont}{stop}.")";
				}
			}else{
				foreach my $curcont (keys %{$$curseg{contigs}})
				{
					$allUsedCont{$curcont} = 1;
					if($contstr){
						$contstr .= '-';
					}
					$contstr .= $curcont."(".$$curseg{contigs}{$curcont}{start}."-".$$curseg{contigs}{$curcont}{stop}.")";
				}
				$contstr .= " (all identical)";
			}
			push(@contseq, $$curseg{refstart}."\t".$$curseg{refstop}."\t$contstr");
	#		die;
		}
	}
	foreach my $cseq (@contseq)
	{
		print CONTIGLIST $cseq."\n";
	}
	print CONTIGLIST "\n*****************\nNB CONTIG USED:\t".scalar(keys %allUsedCont)."\n\n";
	foreach my $allusedcont (sort keys %allUsedCont)
	{
		print CONTIGLIST $allusedcont."\n";
	}
	return $consseq;
}

sub mergeContigsNoRead
{
	my $consseq = '';
	my @contseq;
	my $prevstop = 0;
	my $prevstart = 0;
	my $prevcont = '';

	my %allUsedCont;
	
	foreach my $curseg(@segMap)
	{
#		print "Merger : ".$$curseg{refstart}."\t".$$curseg{refstop}."\n";
		if($$curseg{refstart} > $prevstop + 1 && $prevstop > 0)
		{
			$consseq .= 'N' x ($$curseg{refstart} - $prevstop - 1);
			push(@contseq, ($prevstop + 1)."\t".($$curseg{refstart} - 1)."\t".($$curseg{refstart} - $prevstop - 1)."N");
		}
		$prevstop = $$curseg{refstop};
		
		my $segdone = 0;
		if($$curseg{nbcontigs} > 1 && $prevcont ne '')
		{
#			print "Contigs : ".$$curseg{allcont}."\n";
			my $count = 1;
			my $midpoint = int($$curseg{refstart} + ($$curseg{refstop} - $$curseg{refstart})/2);
#			print "Midpoint : ".$midpoint."\n";
			my @ovContigs;

			foreach my $curcont (sort {$$curseg{contigs}{$b}{length} <=> $$curseg{contigs}{$a}{length}} keys %{$$curseg{contigs}})
			{
#				print "MULTI : ".$curcont."\t".$$curseg{contigs}{$curcont}{start}."\t".$$curseg{contigs}{$curcont}{stop}."\n";
				if($curcont eq $prevcont)
				{
					$ovContigs[0] = $curcont;
				}
			}
			if($ovContigs[0])
			{
				foreach my $curcont (sort {$$curseg{contigs}{$b}{length} <=> $$curseg{contigs}{$a}{length}} keys %{$$curseg{contigs}})
				{
					if($curcont ne $ovContigs[0])
					{
						$ovContigs[1] = $curcont;
					}
				}
				
				(my $contstart, my $contstop, my $contseq) = getSegmentSubstr($ovContigs[0], $$curseg{refstart}, $midpoint);
#				print $ovContigs[0]."\t".$$curseg{refstart}."\t".$midpoint."\t$contstart\t$contstop\n";
				$consseq .= $contseq;
				push(@contseq, $$curseg{refstart}."\t".$midpoint."\t".$ovContigs[0]." ".$contstart."-".$contstop);
				$allUsedCont{$ovContigs[0]} = 1;

				($contstart, $contstop, $contseq) = getSegmentSubstr($ovContigs[1], ($midpoint+1), $$curseg{refstop});
#				print $ovContigs[1]."\t".($midpoint+1)."\t".$$curseg{refstop}."\t$contstart\t$contstop\n";
				$consseq .= $contseq;
				push(@contseq, ($midpoint + 1)."\t".$$curseg{refstop}."\t".$ovContigs[1]." ".$contstart."-".$contstop);
				$prevcont = $ovContigs[1];
				$allUsedCont{$ovContigs[1]} = 1;
				$segdone = 1;
			}
		}
		if($segdone){next;}

		foreach my $curcont (sort {$$curseg{contigs}{$b}{length} <=> $$curseg{contigs}{$a}{length}} keys %{$$curseg{contigs}})
		{
			my $contstart = $$curseg{contigs}{$curcont}{start};
			my $contstop = $$curseg{contigs}{$curcont}{stop};
#			print "SINGLE : ".$curcont."\t".$$curseg{contigs}{$curcont}{start}."\t".$$curseg{contigs}{$curcont}{stop}."\n";

			if($$curseg{refstart} == 1){
				$contstart = 1;
			}
			if($$curseg{refstop} == length($refseq))
			{
				$contstop = length($contigSeqs{$curcont}{seq});
			}

#			print ">$curcont\n".$contigSeqs{$curcont}{seq}."\n";
			my $curcontseq = substr($contigSeqs{$curcont}{seq}, ($contstart - 1), ($contstop - $contstart + 1));

			$consseq .= $curcontseq;
			push(@contseq, $$curseg{refstart}."\t".$$curseg{refstop}."\t".$curcont." ".$contstart."-".$contstop);
			$allUsedCont{$curcont} = 1;
			$prevcont = $curcont;
			last;
		}
	}

	foreach my $cseq (@contseq)
	{
		print CONTIGLIST $cseq."\n";
	}
	print CONTIGLIST "\n*****************\nNB CONTIG USED:\t".scalar(keys %allUsedCont)."\n\n";
	foreach my $allusedcont (sort keys %allUsedCont)
	{
		print CONTIGLIST $allusedcont."\n";
	}
	return $consseq;
}


sub buildOverlapMap
{
	my %validSegments;
	my $count = 0;
	my $previousstart = 10000000;
	my $previousstop = 0;
	foreach my $cursegcont (sort keys %segments)
	{
		foreach my $segno (sort {$segments{$cursegcont}{$a}{start} <=> $segments{$cursegcont}{$b}{start}} keys %{$segments{$cursegcont}})
		{
			if($segments{$cursegcont}{$segno}{end} - $segments{$cursegcont}{$segno}{start} + 1 >= $option{minseglen})
			{	
				unless($option{readfa} && $option{readq})
				{
					if($segments{$cursegcont}{$segno}{start} >= $previousstart && $segments{$cursegcont}{$segno}{end} <= $previousstop){next;}
				}
				$validSegments{$count}{start} = $segments{$cursegcont}{$segno}{start};
				$validSegments{$count}{end} = $segments{$cursegcont}{$segno}{end};
				$validSegments{$count}{contig} = $cursegcont;
				$validSegments{$count}{segno} = $segno;
				$previousstart = $validSegments{$count}{start};
				$previousstop = $validSegments{$count}{end};
				$count++;
			}
		}
	}
	undef(%segOverlap);
	my $overcount = 1;
	foreach my $valseg (sort {$validSegments{$a}{start} <=> $validSegments{$b}{start}} keys %validSegments)
	{
		my $hasover = 0;
#		print "VALSEG : ".$validSegments{$valseg}{contig}."\t".$validSegments{$valseg}{start}."\t".$validSegments{$valseg}{end}."\t$hasover\n";
		my $curmin = $validSegments{$valseg}{start};
		my $curmax = $validSegments{$valseg}{end};
		foreach my $valseg2 (sort {$validSegments{$a}{start} <=> $validSegments{$b}{start}} keys %validSegments)
		{
			if($valseg == $valseg2){next;}
			if($validSegments{$valseg2}{start} > $validSegments{$valseg}{end}){last;}
			
			my $arrref = getOverlap($validSegments{$valseg}{start}, $validSegments{$valseg}{end}, $validSegments{$valseg2}{start}, $validSegments{$valseg2}{end});
			unless(scalar(@$arrref) == 0){$hasover = 1;}
			
			my @curOverlap = @$arrref;
			
			foreach my $frag (@curOverlap)
			{
				# Check if overlap fragment exists
				my $newover = 1;
				foreach my $segover (keys %segOverlap)
				{
					if($segOverlap{$segover}{start} == $$frag{start} && $segOverlap{$segover}{end} == $$frag{end})
					{
						$newover = 0;
					}
				}
				if($newover)
				{
					$segOverlap{$overcount}{start} = $$frag{start};
					$segOverlap{$overcount}{end} = $$frag{end};
					$overcount++;
				}
			}
		}
		if($hasover == 0)
		{
			$segOverlap{$overcount}{start} = $validSegments{$valseg}{start};
			$segOverlap{$overcount}{end} = $validSegments{$valseg}{end};
			$overcount++;
#			print $validSegments{$valseg}{contig}." HAS NO OVER ($overcount)!\n";
		}
	}
	
	# Get the list of all overlap start/stops
	my @starts;
	my @stops;
	foreach my $segover (sort {$segOverlap{$a}{start} <=> $segOverlap{$b}{start}} keys %segOverlap)
	{
		my %usednames;
		push(@starts, $segOverlap{$segover}{start});
		push(@stops, $segOverlap{$segover}{end});
#		print "OVERLAP : ".$segOverlap{$segover}{start}."\t".$segOverlap{$segover}{end}."\n";
	}

	# Determine shortest segments
	@starts = sort {$a <=> $b} (@starts);
	@stops = sort {$a <=> $b} (@stops);
	my @shortestseg;
	my $prevstart = 0;
	my $prevstop = 0;
	my $segflag = 0;
	while(scalar(@starts) > 0)
	{
		my $curstart = shift(@starts);
		if($curstart == $prevstart){next;}
		$prevstart = $curstart;
		$shortestseg[$segflag]{start} = $curstart;
		my $curstop = 0;
		while($curstop < $curstart && scalar(@stops) > 0)
		{
			$curstop = shift(@stops);
		}
		
		$shortestseg[$segflag]{stop} = $curstop;
		$segflag++;
	}

	# Get contigs spanning each segment
	foreach my $shortseg (@shortestseg)
	{
		my $curcontlist = '';
		foreach my $curseg (keys %validSegments)
		{
			if($validSegments{$curseg}{start} <= $$shortseg{start} && $validSegments{$curseg}{end} >= $$shortseg{stop})
			{
				if($curcontlist){$curcontlist .= "/";}
				$curcontlist .= $validSegments{$curseg}{contig};#.".".$validSegments{$curseg}{segno};
			}
		}
		$$shortseg{contigs} = $curcontlist;
	}
	
	my $arrflag = 0;
	foreach my $shortseg (@shortestseg)
	{
#		print "SHORTSEG: ".$$shortseg{start}."\t".$$shortseg{stop}."\t".$$shortseg{contigs}."\n";
		my @contiglist = split(/\//, $$shortseg{contigs});
		$segMap[$arrflag]{refstart} = $$shortseg{start};
		$segMap[$arrflag]{refstop} = $$shortseg{stop};
		$segMap[$arrflag]{allcont} = $$shortseg{contigs};
		foreach my $curcontlist (@contiglist)
		{
			(my $contstart, my $contstop, my $contseq) = getSegmentSubstr($curcontlist, $$shortseg{start}, $$shortseg{stop});
#			print "CONTLIST : ".$contstart."\t".$contstop."\t$curcontlist\n";
			$segMap[$arrflag]{contigs}{$curcontlist}{start} = $contstart;
			$segMap[$arrflag]{contigs}{$curcontlist}{stop} = $contstop;
			$segMap[$arrflag]{contigs}{$curcontlist}{seq} = $contseq;
			$segMap[$arrflag]{contigs}{$curcontlist}{length} = $contstop - $contstart + 1;
			$segMap[$arrflag]{nbcontigs}++;
		}
		$arrflag++;
	}
}

sub getSegmentSubstr
{
	my $contigname = shift;
	my $start = shift;
	my $stop = shift;
	
	open(CONTALN, $orientContOut."_".$contigname."_orientalign.afa");
	
	my $contalnseq = '';
	my $refalnseq = '';
	my $mode = '';
	
	while(my $line = <CONTALN>)
	{
		chomp $line;
		if($line =~ />$contigname/)
		{
			$mode = 'cont';
		}elsif($line =~ />.+/)
		{
			$mode = 'ref';
		}elsif($mode eq 'cont')
		{
			$contalnseq .= $line;
		}elsif($mode eq 'ref')
		{
			$refalnseq .= $line;
		}
	}
	
	my $contpos = 0;
	my $refpos = 0;
	
	my $contstartpos = 0;
	my $contstoppos = 0;
	
	for(my $pos = 0; $pos < length($contalnseq); $pos++)
	{
		my $curcont = substr($contalnseq, $pos, 1);
		my $curref = substr($refalnseq, $pos, 1);
		
		if($curcont ne '-'){$contpos++;}
		if($curref ne '-'){$refpos++;}
		
		if($refpos == $start)
		{
			$contstartpos = $contpos;
		}elsif($refpos == $stop)
		{
			$contstoppos = $contpos;
			last;
		}
	}
	
	return($contstartpos, $contstoppos, substr($contigSeqs{$contigname}{seq}, $contstartpos - 1, $contstoppos-$contstartpos+1));
}

sub getOverlap
{
	my $start1 = shift;
	my $end1 = shift;
	my $start2 = shift;
	my $end2 = shift;
	
	my @overlapFrags;
	# Check if no overlap
	if($end1 < $start2 || $end2 < $start1)
	{
		return \@overlapFrags;
	}
	
	my $overnb = 0;
	if($start1 == $start2)
	{
		$overlapFrags[0]{start} = $start1;
		$overlapFrags[0]{contig} = "both";
	}elsif($start1 < $start2)
	{
		$overlapFrags[$overnb]{start} = $start1;
		$overlapFrags[$overnb]{end} = $start2 - 1;
		$overlapFrags[$overnb]{contig} = "c1";
		$overnb++;
		$overlapFrags[$overnb]{start} = $start2;
		$overlapFrags[$overnb]{contig} = "both";
	}elsif($start1 > $start2)
	{
		$overlapFrags[$overnb]{start} = $start2;
		$overlapFrags[$overnb]{end} = $start1 - 1;
		$overlapFrags[$overnb]{contig} = "c2";
		$overnb++;
		$overlapFrags[$overnb]{start} = $start1;
		$overlapFrags[$overnb]{contig} = "both";
	}

	if($end1 == $end2)
	{
		$overlapFrags[$overnb]{end} = $end2;
	}elsif($end1 < $end2)
	{
		$overlapFrags[$overnb]{end} = $end1;
		$overnb++;
		$overlapFrags[$overnb]{start} = $end1 + 1;
		$overlapFrags[$overnb]{end} = $end2;
		$overlapFrags[$overnb]{contig} = "c2";
	}else{
		$overlapFrags[$overnb]{end} = $end2;
		$overnb++;
		$overlapFrags[$overnb]{start} = $end2 + 1;
		$overlapFrags[$overnb]{end} = $end1;
		$overlapFrags[$overnb]{contig} = "c1";
	}
	
	return \@overlapFrags;
}


sub buildSegments
{
	my %tmpSegments;
	foreach my $contig (keys %contigSeqs)
	{
		my $maxInsertLen = $option{maxsegins};
		open(AFA, $orientContOut."_".$contig."_orientalign.afa");
		my $contalignseq = '';
		my $refalignseq = '';
		my $curseq = '';
		while(my $line = <AFA>)
		{
			chomp $line;
			if($line =~ />(.+)/)
			{
				if($1 eq $contig)
				{
					$curseq = 'cont';
				}else{$curseq = 'ref';}
			}else{
				if($curseq eq 'cont'){$contalignseq .= $line;}
				else{$refalignseq .= $line;}
			}
		}
		
		my $refpos = 0;
		my $ingap = 0;
		my $gaplen = 0;
		my $segno = 0;
		my $inseg = 0;
		my $refstart = 0;
		my $segstart = 0;
		my $segend = 0;
		my $segpos = 0;
		my $insertlen = 0;
		my $curseggap = 0;
		
		for(my $pos = 0; $pos < length($refalignseq); $pos++)
		{
			my $curref = substr($refalignseq, $pos, 1);
			my $curcont = substr($contalignseq, $pos, 1);
			if($curcont ne '-'){$segpos++;}
			if($curref ne '-'){
				if($insertlen >= $maxInsertLen){
					$segInserts{$contig}{$segno}{refpos} = $refpos;
					$segInserts{$contig}{$segno}{seq} = substr($contalignseq, $pos - $insertlen, $insertlen);
					$inseg = 0;
				}
				$insertlen = 0;
				$refstart=1;
				$refpos++;
				$contigAlign{$contig}{$refpos} = $segpos;
			}else{
				if($refstart)
				{
					$insertlen++;
				}
			}
			unless($refstart){next;}
			if($curcont ne '-'){
				$gaplen = 0;
				if($inseg == 0 && $insertlen == 0)
				{
					$segno++;
					$inseg = 1;
					unless($segments{$contig}){$nbCont++;}
					$tmpSegments{$contig}{$segno}{start} = $refpos;
					$tmpSegments{$contig}{$segno}{segstart} = $segpos;
					$tmpSegments{$contig}{$segno}{end} = $refpos;
					$tmpSegments{$contig}{$segno}{segend} = $segpos;
					$tmpSegments{$contig}{$segno}{alnpos} = 1;
					if($curref eq '-')
					{
						$tmpSegments{$contig}{$segno}{nbgap} = 1;
					}else{
						$tmpSegments{$contig}{$segno}{nbgap} = 0;
					}
				}elsif($insertlen == 0){
					$tmpSegments{$contig}{$segno}{segend} = $segpos;
					$tmpSegments{$contig}{$segno}{end} = $refpos;
					if($curref eq '-')
					{
						$tmpSegments{$contig}{$segno}{nbgap}++;
					}
					$tmpSegments{$contig}{$segno}{alnpos}++;
				}
			}else{
				if($inseg){
					$gaplen++;
					$tmpSegments{$contig}{$segno}{nbgap}++;
					$tmpSegments{$contig}{$segno}{alnpos}++;
					if($gaplen > $maxgaplen){
						$inseg = 0;
					}
				}
			}
		}	
	}
	
	#Filter out segments too short or with too high internal gap
	foreach my $contig(keys %tmpSegments)
	{
		foreach my $seg (sort keys %{$tmpSegments{$contig}})
		{
			my $segdelpct = $tmpSegments{$contig}{$seg}{nbgap}/$tmpSegments{$contig}{$seg}{alnpos};
			my $seglen = $tmpSegments{$contig}{$seg}{end} - $tmpSegments{$contig}{$seg}{start} + 1;
			if($segdelpct <= $option{maxsegdel} && $seglen >= $option{minseglen})
			{
				$segments{$contig}{$seg}{start} = $tmpSegments{$contig}{$seg}{start};
				$segments{$contig}{$seg}{segstart} = $tmpSegments{$contig}{$seg}{segstart};
				$segments{$contig}{$seg}{end} = $tmpSegments{$contig}{$seg}{end};
				$segments{$contig}{$seg}{segend} = $tmpSegments{$contig}{$seg}{segend};
#				print $contig."\t".$seg."\t".$tmpSegments{$contig}{$seg}{start}."-".$tmpSegments{$contig}{$seg}{end}."\t".$tmpSegments{$contig}{$seg}{segstart}."-".$tmpSegments{$contig}{$seg}{segend}."\t$segdelpct"."(".$tmpSegments{$contig}{$seg}{nbgap}."/".$tmpSegments{$contig}{$seg}{alnpos}.")\n";
			}
		}
	}
	
	# Select only segments required to add coverage
	my %covered;
#	print "\n-----\n";
	foreach my $contig(sort {$contigSeqs{$b}{length} <=> $contigSeqs{$a}{length}} keys %contigSeqs)
	{
		if($segments{$contig} && !($option{allcont}))
		{
			# Determine if Contig should be added
			my $addCont = 0;
			foreach my $seg (keys %{$segments{$contig}})
			{
				for(my $i = $segments{$contig}{$seg}{start}; $i <= $segments{$contig}{$seg}{end}; $i++)
				{
					if(!($covered{$i}) || ($contigSeqs{$contig}{length} >= $option{longcont} && $covered{$i} < $option{contigcov}))
					{
#						print $contig."\t$seg\t$i\n";
						$addCont = 1;
						last;
					}
				}
			}
			# If yes, add coverage for all segments of Contig
#			print "addCont $addCont\n";
			if($addCont)
			{
				foreach my $seg (keys %{$segments{$contig}})
				{
#					print $contig."\t".$seg."\t".$contigSeqs{$contig}{length}."\n";
					for(my $i = $segments{$contig}{$seg}{start}; $i <= $segments{$contig}{$seg}{end}; $i++)
					{
						$covered{$i}++;
					}
				}
			}else{
				delete $segments{$contig};
			}
		}
	}
#	die;
}

sub drawSegments
{
	my @segcolors = ("blue","green","orange","brown","pink","purple");

	open(ROUTPUT, ">$output"."_contigsMap.R");
	print ROUTPUT "x=c(";

	#print X axis values
	for(my $i = 1; $i < length($refseq); $i++)
	{
		print ROUTPUT $i.",";
	}
	print ROUTPUT length($refseq).")\n";
	print ROUTPUT "y=c(0";
	
	# Fill coverage matrix
	my $maxycov = 0;
	my %posCovered;
	for(my $curpos = 2; $curpos <= length($refseq); $curpos++)
	{
#		if($coverages{$curpos})
#		{
#			if($coverages{$curpos}{allcov})
#			{
#				print ROUTPUT " ,".$coverages{$curpos}{allcov};
#				if($coverages{$curpos}{allcov} > $maxycov){$maxycov = $coverages{$curpos}{allcov};}
#			}else{
		print ROUTPUT " ,0";
#			}
#		}else{
#			print ROUTPUT " ,0";
#		}
	}
	print ROUTPUT ")\n";
	my $contigSpace = 1;
#	my $ampSpace = int(0.01 * $maxycov);
#	if($ampSpace < 3){$ampSpace = 3};
	my $miny = (scalar(keys %segments)*(-$contigSpace));
	if($miny > -5){$miny = -5;}
	print ROUTPUT "plot(x,y,col=\"red\", xlab=\"Position\",ylab=\"Contigs\",type=\"l\",lwd=2,ann=T,cex.lab=0.8, ylim=c(".$miny.",max(y)), main=\"$output\")\n";

	# Draw contigs
	my $contflag = 1;
	my $used = 0;
	foreach my $contno (sort {$segments{$a} <=> $segments{$b}} keys %segments)
	{
		my $nbseg = 1;
		my $lastsegend = 0;
		foreach my $seg (sort {$segments{$contno}{$a}{start} <=> $segments{$contno}{$b}{start}} keys %{$segments{$contno}})
		{
			my $segdelpct = (($segments{$contno}{$seg}{end} - $segments{$contno}{$seg}{start} + 1) - ($segments{$contno}{$seg}{segend} - $segments{$contno}{$seg}{segstart} + 1))/($segments{$contno}{$seg}{end} - $segments{$contno}{$seg}{start} + 1);
			
			if(($segments{$contno}{$seg}{end} - $segments{$contno}{$seg}{start} + 1) >= $mincontlen && $segdelpct <= $option{maxsegdel})
			{
				if($nbseg > 1)
				{
					print INDELS "$contno\t".$contigAlign{$contno}{($lastsegend+1)}."\t".$contigAlign{$contno}{$segments{$contno}{$seg}{start}}."\t".($lastsegend+1)."\t".($segments{$contno}{$seg}{start}-1)."\n";
				}
				$used = 1;
				print ROUTPUT "segments(".$segments{$contno}{$seg}{start}.",".($contflag*(-$contigSpace)).",".$segments{$contno}{$seg}{end}.",".($contflag*(-$contigSpace)).",col=\"".$segcolors[($contflag+5)%6]."\",lwd=2)\n";
#				print "ContNo : $contno\tsegments(".$segments{$contno}{$seg}{start}.",".($contflag*(-5)).",".$segments{$contno}{$seg}{end}.",".($contflag*(-5)).",col=\"".$segcolors[($contflag+5)%6]."\",lwd=2)\n";
				for(my $i = $segments{$contno}{$seg}{start}; $i <= $segments{$contno}{$seg}{end}; $i++)
				{
					$posCovered{total}{$i}++;
					$posCovered{$contno}{$i}++;
					$posCovered{$contno}{length}++;
				}
				$lastsegend = $segments{$contno}{$seg}{end};
				$nbseg++;
			}
		}
		if($used){
			$contflag++;
			$used = 0;
		}
	}

	# Draw amplicons
#	$ampno = 0;
#	foreach my $curamp(@refAmps)
#	{
#		print ROUTPUT "segments(".$$curamp{start}.",max(y)+".((scalar(@refAmps) - $ampno) * $ampSpace).",".$$curamp{end}.",max(y)+".((scalar(@refAmps) - $ampno) * $ampSpace).",col=\"black\",lwd=2)\n";
#		$ampno++;
#	}


	my $pdfname = "$output"."_contigsMap.pdf";
	print ROUTPUT "pdf(\"$pdfname\")\n";
	system($Rpath."R < $output"."_contigsMap.R --no-save --slave -q");
	`mv Rplots.pdf $pdfname`;
	
#	foreach my $contig(keys %contigSeqs)
#	{
#		system("rm $output"."_$contig"."*");
#	}
}


sub readOrientContig
{
	my $inputcont = shift;
	
	open(INPUTCONT, $inputcont) || die("invalid contig file\n");
	my $contId = '';
	my $contSeq = '';
	while(my $line = <INPUTCONT>)
	{
		chomp $line;
		if($line =~ />(.+)/)
		{
			if($contSeq)
			{
				$contigSeqs{$contId}{seq}  = $contSeq;
				$contigSeqs{$contId}{length} = length($contSeq);
				$contSeq = '';
			}
			$contId = $1;
		}else{
			$contSeq .= $line;
		}
	}
	$contigSeqs{$contId}{seq} = $contSeq;
	$contigSeqs{$contId}{length} = length($contSeq);
}


sub readFasta
{
	my $fastafile = shift;
	open(FASTA, $fastafile) || die("invalid $fastafile file\n");

	my $id = '';
	my $seq = '';
	
	while(my $line = <FASTA>)
	{
		chomp $line;
		if($line =~ />(.+)/)
		{
			$id = $1;
		}else{
			$seq .= $line;
		}
	}
	return $id, $seq;
}

sub fillNtFreq
{
	my $qlxfile = shift;
	my $curid = shift;
	
	open(QLXFILE, $qlxfile);

	my $refstart;
	my $refend;
	my $readstring;
	my $refstring;
	my $qualstring;
	my $state = 0;

	while (my $line = <QLXFILE>)
	{
		chomp $line;
		if($line =~ />/)
		{
			my @query = split(/\s/, $line);
			$refstart = $query[7];
			$refend = $query[8];
			$readstring = "";
			$refstring = "";
			$qualstring = "";
			$state = 1;
		}elsif($line =~ /[ATGC]/ && $state == 1)
		{
			$refstring .= $line;
			$state = 2;
		}elsif($line =~ /[ATGC]/ && $state == 2)
		{
			$readstring .= $line;
			$state = 3;
		}elsif($line =~ /[0-9]/ && $state == 3)
		{
			$qualstring = $line;
			evalReadNtFreq($readstring, $refstring, $refstart, $curid);
		}
	}
}

sub evalReadNtFreq
{
	my $readstr = shift;
	my $refstr = shift;
	my $refstart = shift;
	my $curid = shift;

	#READ NTS
	my $ntpos = $refstart;
	my $curinsert = '';
	my $curdel = '';
	my $delstart = 0;
	my $insertstart = 0;
	
	for(my $currefpos = 0; $currefpos < length($refstr); $currefpos++)
	{
		my $curnt = substr($readstr, $currefpos, 1);
	#	my $curntqual = substr($qualstr, $currefpos, 1);
		my $curref = substr($refstr, $currefpos, 1);
		
		if($curref eq "-")
		{
			$curinsert .= $curnt;
		}elsif($curnt eq "-")
		{
			unless($curdel){$delstart = $ntpos;}
			$curdel .= '-';
		}else{
			if($curinsert)
			{
				$ntfreq{$curid}{$ntpos - 1}{I}{$curinsert}++;
				$ntfreq{$curid}{$ntpos - 1}{I}{count}++;
			}
			$curinsert = '';

			if($curdel)
			{
				my $del = "D".length($curdel);
				for(my $delpos = $delstart; $delpos < $ntpos; $delpos++)
				{
					if($delpos == $delstart)
					{
						$ntfreq{$curid}{$delpos}{'-start'}{count}++;
						$ntfreq{$curid}{$delpos}{'-'}{count}++;
						$ntfreq{$curid}{$delpos}{$del}{count}++;
					}else{
						$ntfreq{$curid}{$delpos}{'-'}{count}++;
					}
					$ntfreq{$curid}{$delpos}{coverage}{count}++;
				}
				$curdel = '';
				$delstart = 0;
			}
			$ntfreq{$curid}{$ntpos}{$curnt}{count}++;
			$ntfreq{$curid}{$ntpos}{coverage}{count}++;
		}
		if($curref ne "-")
		{
			$ntpos++;
		}
	}
}


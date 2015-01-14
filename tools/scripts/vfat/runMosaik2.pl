#!/usr/bin/perl perl -w

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
# This wrapper is designed to run the software Mosaik which can be found at 
# http://bioinformatics.bc.edu/marthlab/Mosaik

use strict;
use Getopt::Long;
use warnings;
use File::Basename;

my %option = (  
	hs			=> 10,
	act			=> 15,
	mmp			=> 0.25,
	minp		=> 0.25,
	mms			=> -9,
	ms			=> 10,
	hgop 		=> 20,
	gop			=> 40,
	gep			=> 10,
	nqsmq		=> 20,
	nqsaq		=> 15,
	nqssize   	=> 5,
	nqvalue		=> 1,
	bw			=> 29,
	m			=> "unique",
	st			=> "illumina",	#sequencing technology. Accepted : sanger, 454, illumina, illumina_long, helicos, solid
	fakequals	=> 0,
	mfl 		=> 600, 	#median fragment length
	ref			=> '',
	fa			=> '',
	qual		=> '',
	fa2			=> '',
	qual2		=> '',
	fq			=> '',
	fq2			=> '',
	o			=> '',
	paramillu	=> '',
	annpe 		=> '',
	annse 		=> '',
	mosaikpath  => "/gsap/garage-viral/viral/analysis/xyang/external_programs/MOSAIK-2.1.33-source/bin",
	mosaiknetworkpath => "/gsap/garage-viral/viral/analysis/xyang/external_programs/MOSAIK-2.1.33-source/networkFile",
);

GetOptions(
  "hs=i"		=> \$option{hs},
  "act=i"		=> \$option{act},
  "mmp=f"		=> \$option{mmp},
  "minp=f"		=> \$option{minp},
  "mms=f"		=> \$option{mms},
  "ms=f"		=> \$option{ms},
  "hgop=f"		=> \$option{hgop},
  "gop=f"		=> \$option{gop},
  "gep=f"		=> \$option{gep},
  "nqsmq=i"		=> \$option{nqsmq},
  "nqsaq=i"		=> \$option{nqsaq},
  "nqssize=i"	=> \$option{nqssize},
  "nqvalue=i"	=> \$option{nqvalue},
  "bw=i"		=> \$option{bw},
  "m=s"			=> \$option{m},
  "st=s"		=> \$option{st},
  "fakequals=i" => \$option{fakequals},
  "mfl=s"		=> \$option{mfl},
  "ref=s"		=> \$option{ref},
  "fa=s"		=> \$option{fa},
  "qual=s"		=> \$option{qual},
  "fa2=s"		=> \$option{fa2},
  "qual2=s"		=> \$option{qual2},
  "fq=s"		=> \$option{fq},
  "fq2=s"		=> \$option{fq2},
  "o=s"			=> \$option{o},
  "paramillu" 	=> \$option{paramillu},
  "annpe=s"		=> \$option{annpe},
  "annse=s"		=> \$option{annse},
  "mosaikpath=s" => \$option{mosaikpath},
  "mosaiknetworkpath=s" => \$option{mosaiknetworkpath}
) || die("Problem processing command-line options: $!\n");

my $output = $option{o};

unless($option{ref} && $option{o} && ($option{fq} || ($option{fa} && $option{qual})))
{
	printHelp();
}
if($option{h})
{
	printHelp();
}

my $scriptpath = dirname(__FILE__);
my $mosaikpath = $option{mosaikpath};
my $mosaiknetworkpath = $option{mosaiknetworkpath};

if($mosaiknetworkpath && substr($mosaiknetworkpath, length($mosaiknetworkpath) - 1) ne '/'){$mosaiknetworkpath .= '/';}
if($mosaikpath && substr($mosaikpath, length($mosaikpath) - 1) ne '/'){$mosaikpath .= '/';}
if($scriptpath && substr($scriptpath, length($scriptpath) - 1) ne '/'){$scriptpath .= '/';}

unless($option{annpe})
{
	$option{annpe} = $mosaiknetworkpath."2.1.26.pe.100.0065.ann"
}
unless($option{annse})
{
	$option{annse} = $mosaiknetworkpath."2.1.26.se.100.005.ann"
}


if($option{paramillu})
{
	$option{gop} = 40;
	$option{hgop} = 20;
	$option{gep} = 10;
}

my $refdat = $output.".refdat";
my $readdat = $output.".readdat";

my $uniquesort = '';
if($option{m} eq 'unique')
{
	$uniquesort = '-u';
}
 
my $bw = '';
if($option{bw})
{
	$bw = " -bw ".$option{bw};
}

#Build Ref and Reads .dat with Mosaik

system($mosaikpath."MosaikBuild -fr ".$option{ref}." -oa $refdat");

if($option{fq} && $option{fq2})
{
	system($mosaikpath."MosaikBuild -q ".$option{fq}." -q2 ".$option{fq2}." -out $readdat -st ".$option{st}." -mfl ".$option{mfl});
}elsif($option{fq})
{
	system($mosaikpath."MosaikBuild -q ".$option{fq}." -out $readdat -st ".$option{st});
}elsif($option{fa} && $option{qual} && $option{fa2} && $option{qual2})
{
	system($mosaikpath."MosaikBuild -fr ".$option{fa}." -fr2 ".$option{fa2}." -fq ".$option{qual}." -fq2 ".$option{qual2}." -out $readdat -st ".$option{st}." -mfl ".$option{mfl});
}elsif($option{fa} && $option{qual})
{
	system($mosaikpath."MosaikBuild -fr ".$option{fa}." -fq ".$option{qual}." -out $readdat -st ".$option{st});
}else{
	die("Must enter either a fastq file (-fq), paired fq files (-fq -fq2), reads and quals (-fa -qual) or paired reads and quals (-fa -fa2 -qual -qual2) for read inputs\n");
}

system($mosaikpath."MosaikAligner -in $readdat -out $output -ia $refdat -hs ".$option{hs}." -act ".$option{act}." -mm 500 -mmp ".$option{mmp}." -minp ".$option{minp}." -ms ".$option{ms}." -mms ".$option{mms}." -gop ".$option{gop}." -hgop ".$option{hgop}." -gep ".$option{gep}."$bw -m ".$option{m}." -annpe ".$option{annpe}." -annse ".$option{annse});

my $sam2qlxopt = '';
if($option{fakequals})
{
	$sam2qlxopt .= " -fakequals ".$option{fakequals};
}
system("perl ".$scriptpath."samToQlx.pl $output.sam ".$option{ref}." $output -nqsmq ".$option{nqsmq}." -nqsaq ".$option{nqsaq}." -nqssize ".$option{nqssize}." -nqvalue ".$option{nqvalue}.$sam2qlxopt);
system("rm $output.sam");

system("rm $output.readdat");
system("rm $output.refdat");
system("rm $output.stat");

sub printHelp
{
	print "Usage : perl runMosaik2.pl -ref reference.fa -o outputBaseName [REQUIRES ONE OF : -fq reads.fq OR -fa reads.fa -qual reads.qual]\n";
	print "Requires Mosaik (currently /seq/viral/analysis/xyang/external_programs/MOSAIK-2.1.33-source/bin/, modify ".'$mosaikpath'." if changes)\n";
	print "Requires Samtools (use Samtools)\n";
	print "Options (see Mosaik documentation for more details):\n";
	print "-hs : set hashsize\n";
	print "-act : specifies the alignment candidate threshold\n";
	print "-mmp : maximum mismatch percentage (in fraction, i.e. 0.25)\n";
	print "-minp : minimum percentage of read length that is aligned (in fraction, i.e. 0.25)\n";
	print "-ms : match score (10 default)\n";
	print "-mms : mismatch penalty (negative number, -9 default)\n";
	print "-gop : gap opening penalty (40 default)\n";
	print "-hgop : homopolymer gap opening penalty (20 default)\n";
	print "-gep : gap extension penalty (10 default)\n";
	print "-nqsmq : NQS (Neighborhood Quality Score) main base minimum quality (20 default)\n";
	print "-nqsaq : NQS adjacent bases minimum quality (10 default)\n";
	print "-nqssize : NQS window size for adjacent bases (5 default)\n";
	print "-nqvalue : Quality value for inserted Ns that should not influence adjacent bases NQS score (1 default)\n";
	print "-bw : banded Smith-Waterman bandwith (29 default)\n";
	print "-m : alignment mode (default unique)\n";
	print "-st : sequencing technology (default illumina)\n";
	print "-fakequals : fakes qualities in the qlx output to a given amount (off by default)\n";
	print "-mfl : medium fragment length for paired reads\n";
	print "-ref : reference sequence fasta\n";
	print "-fa : reads fasta\n";
	print "-qual : reads qual\n";
	print "-fa2 : reads second mate fasta if paired\n";
	print "-qual2 : reads second mate qual if paired\n";
	print "-fq : reads fastq format\n";
	print "-fq2 : reads fastq second mate if paired\n";
	print "-o : output base name\n";
	print "-paramillu : sets Smith-Waterman parameters for illumina sequencing technology\n";
	print "-annpe : Network file. This file is actually included in the Mosaik distribution. It should be located in the networkFile folder. For Mosaik 2.1.26, this file is named: 2.1.26.pe.100.0065.ann\n";
	print "-annse : Network file. This file is actually included in the Mosaik distribution. It should be located in the networkFile folder. For Mosaik 2.1.26, this file is named: 2.1.26.se.100.005.ann\n";
	print "-h print help\n";
	exit();
}
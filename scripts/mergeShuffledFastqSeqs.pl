#!/usr/bin/perl -w 

#
#    Copyright 2012 Adhemar Zerlotini Neto 
#    Author: Adhemar Zerlotini Neto (azneto@gmail.com)
#
#    mergeShuffledFastqSeqs.pl is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    mergeShuffledFastqSeqs.pl is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
 

use strict;
use Getopt::Long;

my $usage = "

$0 -f1 <file 1> -f2 <file 2> -r <regular expression> -o <output file name> -t
$0 -h

-f1 <file 1>                    : Fastq input file containing sequences identified as 1
-f2 <file 2>                    : Fastq input file containing sequences identified as 2
-r <regular expression>         : Regular expression to locate the ID without the 1 and 2 in each sequence's description line
                                  Eg. '^@(\\S+)/[1|2]\$' to indicate the ID is after the @ symbol and before a / symbol
-o <output file name>           : Prefix name for the result files
                                  Eg. 'mergedSequences'
-t                              : Save two files rather than a interleaved one
-h                              : Help message

";


$| = 1;

my $file1;
my $file2;
my $regex;
my $outname;
my $two;
my $help;

GetOptions ("f1=s" => \$file1,
            "f2=s" => \$file2,
            "r=s" => \$regex,
            "o=s" => \$outname,
            "t!"  => \$two,
            "h!"  => \$help);

if ($help) {
        die $usage;
}

die "\nAll parameters are required",$usage unless( $file1 && $file2 && $regex && $outname);

open (FILEA, "< $file1") or die "\nCouldn't open file1 ($file1).",$usage;
open (FILEB, "< $file2") or die "\nCouldn't open file2 ($file2).",$usage;

if(defined $two) {
 
  open (FASTQ1OUT, "> $outname.1.fastq") or die "\nCouldn't open fastq 1 outfile ($outname.1.fastq)",$usage;
  open (FASTQ2OUT, "> $outname.2.fastq") or die "\nCouldn't open fastq 2 outfile ($outname.2.fastq)",$usage;
  open (NOMATCHOUT, "> $outname.nomatch.fastq") or die "\nCouldn't open no match outfile ($outname.nomatch.fastq)",$usage;
 
} else {
 
  open (FASTQOUT, "> $outname.merged.fastq") or die "\nCouldn't open fastq outfile ($outname.merged.fastq)",$usage;
  open (NOMATCHOUT, "> $outname.nomatch.fastq") or die "\nCouldn't open no match outfile ($outname.nomatch.fastq)",$usage;
 
}


my %idsFileA;
my %idsFileB;

print "\nLoading the first file...";

while(my $line = <FILEA>) {

  chomp $line;

  if($line =~ /$regex/) {
     
     $idsFileA{$1} .=  <FILEA>;
     my $discard = <FILEA>;
     $idsFileA{$1} .=  "+" . $1 . "/1\n";
     $idsFileA{$1} .=  <FILEA>;

  } else {

     die "The regular expression provided doesn't match sequence description.",$usage;

  }

}

print "\nMatching the entries...";

while(my $line = <FILEB>) {

    chomp $line;

    my $id;

    if($line =~ /$regex/) {

       $id = $1;

       $idsFileB{$1} .=  <FILEB>;
       my $discard = <FILEB>;
       $idsFileB{$1} .=  "+" . $id . "/2\n";
       $idsFileB{$1} .=  <FILEB>;

    } else {

       die "The regular expression provided doesn't match sequence description.",$usage;

    }

    if($idsFileA{$id} && $idsFileB{$id}) {

        if(defined $two) {

            print FASTQ1OUT "@" . $id . "/1\n";
            print FASTQ1OUT $idsFileA{$id};

            print FASTQ2OUT "@" . $id . "/2\n";
            print FASTQ2OUT $idsFileB{$id};


        } else {

            print FASTQOUT "@" . $id . "/1\n";
            print FASTQOUT $idsFileA{$id};

            print FASTQOUT "@" . $id . "/2\n";
            print FASTQOUT $idsFileB{$id};

        }

        delete $idsFileA{$id};
        delete $idsFileB{$id};

    } else {

        print NOMATCHOUT "@" . $id . "/2\n";
        print NOMATCHOUT $idsFileB{$id};

        delete $idsFileB{$id};

    }

}

print "\nSaving the remaining sequences...";

for my $id (keys %idsFileA) {

        print NOMATCHOUT "@" . $id . "/1\n";
        print NOMATCHOUT $idsFileA{$id};

        delete $idsFileA{$id};
}

if(defined $two) {

   close(FASTQ1OUT);
   close(FASTQ2OUT);

} else {

   close(FASTQOUT);

}

close(NOMATCHOUT);

print "\nDONE\n\n";



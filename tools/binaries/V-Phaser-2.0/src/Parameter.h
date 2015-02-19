//========================================================================
// Project     : DiversityCaller
// Name        : Parameter.h
// Author      : Xiao Yang
// Created on  : Apr 9, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================


#ifndef PARAMETER_H_
#define PARAMETER_H_

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#include "xutil.h"

class Parameter{
public:
	std::string iDirNm, oDirNm;
	int errModel; // 1: pileup + phase; 2: pileup
	int pSample, windowSz, delta;
	int ignoreBases; // # bases on both end of each read to be neglected
	// indicator value for the corresponding variant to be considered
	int var_qt, var_dt, var_cycle, var_matepair;
	double alpha; // stat sig value
	Parameter(int argc, char** argv): argnum(argc), arg(argv){
        initialize();
		for (int i = 1; i < argnum; i += 2) {
			std::string option = argv[i];
			if (option.compare("-h") == 0 ||	option.compare("--h") == 0) {
				printUsage (argv[0]);
			} else if (option.compare("-i") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				iDirNm = argv[i + 1];
			} else if (option.compare("-o") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				oDirNm = argv[i + 1];

			} else if (option.compare("-e") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				errModel = std::atoi(argv[i + 1]);
				if (errModel != 1 && errModel != 2) {
					std::cout << "-e illegal value\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-w") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				windowSz = std::atoi(argv[i + 1]);
				if (windowSz < 0) {
					std::cout << "-w illegal value\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-delta") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				delta = std::atoi(argv[i + 1]);
				if (delta < 0) {
					std::cout << "-delta illegal value\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-ps") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				pSample = std::atoi(argv[i + 1]);
				if (pSample <= 0 || pSample > 100) {
					std::cout << "-ps illegal value\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-dt") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				var_dt = std::atoi(argv[i + 1]);
				if (var_dt != 0 && var_dt != 1) {
					std::cout << "-dt illegal value\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-cy") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				var_cycle = std::atoi(argv[i + 1]);
				if (var_cycle != 0 && var_cycle != 1) {
					std::cout << "-cy illegal value\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-mp") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				var_matepair = std::atoi(argv[i + 1]);
				if (var_matepair != 0 && var_matepair != 1) {
					std::cout << "-mp illegal value\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-qual") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				var_qt = std::atoi(argv[i + 1]);
				if (var_qt < 0) {
					std::cout << "-qual illegal value\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-a") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				alpha = std::atof(argv[i + 1]);
				if (alpha < 0 || alpha > 100) {
					std::cout << "-a illegal value\n";
					printUsage(argv[0]);
				}
			} else if (option.compare("-ig") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				ignoreBases = std::atof(argv[i + 1]);
				if (ignoreBases < 0) {
					std::cout << "-ig illegal value\n";
					printUsage(argv[0]);
				}
			}
		} // for
		if (iDirNm.empty()) {
			std::cout << "input BAM file is unspecified\n";
			printUsage(argv[0]);
		}
		if (iDirNm.empty()) {
			std::cout << "output DIR is unspecified\n";
			printUsage(argv[0]);
		}

		printSpec();
	}
private:
	int argnum;
	char** arg;

	void initialize () {	// set default values
		alpha = 0.05;
		errModel = 1;
		var_qt = 20;
		var_dt = var_cycle = var_matepair = 1;
		pSample = 30;
		windowSz = 500;
		delta = 2;
		ignoreBases = 0;
	}
	void printUsage(char* exe) {
		std::cout << "\n--------------------------------------------------------\n";
		std::cout << "Usage: " << exe << "\n"
		  << "\t-i  [input.bam] -- input sorted bam file\n"
		  << "\t-o	[output DIR] -- output directory\n"
		  << "\t-e	[1 or 2] -- default 1; 1: pileup + phasing; 2: pileup\n"
		  << "\t-w	-- default 500; alignment window size\n"
		  <<	 "\t-ig	-- default 0; # of bases to ignore on both end of a read\n"
		  << "\t-delta	-- default 2; constrain PE distance by delta x fragsize_variation (auto measured by program)\n"
		  << "\t-ps	(0, 100] -- default 30; percentage of reads to sample to get stats.\n"
		  << "\t-dt	[0 or 1] -- default 1; 1: dinucleotide for err prob measure; 0: not\n"
		  << "\t-cy	[0 or 1] -- default 1; 1: read cycle for err calibr; 0: not\n"
		  << "\t-mp	[0 or 1] -- default 1; 1: mate-pair for err calibr; 0: not\n"
		  << "\t-qual [0, 40] -- default 20; quantile of qual for err calibr\n"
		  << "\t-a	-- default 0.05; significance value for stat test\n";
		std::cout << "----------------------------------------------------------\n\n";
		exit(1);
	}

	void printSpec() {
		std::cout << "\n--------------------------------------------------------\n";
		std::cout << "Program runs with the following Parameter setting:\n\n";
		std::cout << "\tinput BAM file\t=\t" << iDirNm << "\n";
		std::cout << "\toutput Directory\t=\t" << oDirNm << "\n";
		if (errModel == 1) std::cout << "\terrModel\t\t=\tpileup + phase\n";
		else std::cout << "\terrModel\t=\tpileup\n";
		std::cout << "\talpha\t\t=\t" << alpha << "\n";
		std::cout << "\tignoreBases \t=\t" << ignoreBases << "\n";
		std::cout << "\t(var_matepair, var_cycle, var_dt, var_qt)\t=\t"
		  << var_matepair << "," << var_cycle << "," << var_dt << "," << var_qt << "\n";
		std::cout << "\tpSample\t\t=\t" << pSample << "%\n";
		std::cout << "\twindowSz\t=\t" << windowSz << "\n";
		std::cout << "\tdelta\t=\t" << delta << "\n";
		std::cout << "\n--------------------------------------------------------\n\n";
	} // printSpec
};

/* input parameters */
/*
class Parameter{ // the version to read config file

public:
	std::string iDirNm, oDirNm;
	int errModel; // 1: pileup + phase; 2: pileup
	int pSample, windowSz, delta;
	int ignoreBases; // # bases on both end of each read to be neglected
	// indicator value for the corresponding variant to be considered
	int var_qt, var_dt, var_cycle, var_matepair;
	double alpha; // stat sig value
	Parameter(int argc, char** argv): argnum(argc), arg(argv){

		if (argnum != 2) printUsage(argv[0]);
         std::ifstream input(arg[1]);
        if (!input)
        { std::cout << "cannot open " << arg[1] << "\n"; exit(1);}
        initialize();
        getPara (input);
		printSpec();
		input.close();
	}

private:
	int argnum;
	char** arg;
	void printUsage(char* exe) {
		std::cout << "\n--------------------------------------------------------\n";
		std::cout << "Usage: " << exe << " [Config File]\n";
		std::cout << "\tPlease refer to config file for setting params.\n";
		std::cout << "----------------------------------------------------------\n\n";
		exit(1);
	}

	void initialize () {	// set default values
		alpha = 0.05;
		errModel = 1;
		var_qt = 20;
		var_dt = var_cycle = var_matepair = 1;
		pSample = 30;
		windowSz = 500;
		delta = 2;
		ignoreBases = 0;
	}

	void getPara (std::ifstream& input) {
	   std::string line, s1;
	   std::istringstream buf;
		while (std::getline(input, line)) {
			buf.clear();
			buf.str(line);

			if (buf >> s1) {
				if (s1 == "inputDirectory") {
					buf >> iDirNm;
					if (iDirNm.empty()) {
						std::cout << "err: pls specify input directory\n";
						exit(1);
					}
				} else if (s1 == "outputDirectory"){
					buf >> oDirNm;
					if (oDirNm.empty())
						abording("pls specify output directory");
				} else if (s1 == "ignoreBases") {
					buf >> ignoreBases;
				} else if (s1 == "errModel"){
					buf >> errModel;
					if (errModel != 1 && errModel != 2)
						abording ("Specify errModel -- 1: pileup + phase\t"
								"2: pileup");
				} else if (s1 == "alpha") {
					buf >> alpha;
					if (alpha >= 1 || alpha <= 0)
						abording ("alpha should be in the range (0,1)");
				} //else if (s1 == "platform") {
				  //	buf >> platform;
					//if (platform > 3 || platform < 1)
					//	abording ("Specify platform -- 1: Illumina \t"
					//			"2: 454 \t 3: Ion Torrent");
				//}
				else if (s1 == "var_qt") {
					buf >> var_qt;
					if (var_qt < 1) abording ("var_qt has to be > 0");
				} else if (s1 == "var_dt") {
					buf >> var_dt;
					if (var_dt != 0 && var_dt != 1) {
						abording ("var_dt should be set as 0 or 1");
					}
				} else if (s1 == "var_cycle") {
					buf >> var_cycle;
					if (var_cycle != 0 && var_cycle != 1) {
						abording ("var_cycle should be set as 0 or 1");
					}
				} else if (s1 == "var_matepair") {
					buf >> var_matepair;
					if (var_matepair != 0 && var_matepair != 1) {
						abording ("var_matepair should be set as 0 or 1");
					}
				} else if (s1 == "pSample") {
					buf >> pSample;
					if (pSample > 100 || pSample < 1) {
						abording ("pSample should be in the range (0,100]");
					}
				} else if (s1 == "windowSz") {
					buf >> windowSz;
				} else if (s1 == "delta") {
					buf >> delta;
				}
			}// if
		} // while
	} // getPara

	void printSpec() {
		std::cout << "\n--------------------------------------------------------\n";
		std::cout << "Program runs with the following Parameter setting:\n\n";
		std::cout << "\tinputDirectory\t=\t" << iDirNm << "\n";
		std::cout << "\toutputDirectory\t=\t" << oDirNm << "\n";
		if (errModel == 1) std::cout << "\terrModel\t\t=\tpileup + phase\n";
		else std::cout << "\terrModel\t=\tpileup\n";
		std::cout << "\talpha\t\t=\t" << alpha << "\n";
		std::cout << "\tignoreBases \t=\t" << ignoreBases << "\n";
		std::cout << "\t(var_matepair, var_cycle, var_dt, var_qt)\t=\t"
		  << var_matepair << "," << var_cycle << "," << var_dt << "," << var_qt << "\n";
		std::cout << "\tpSample\t\t=\t" << pSample << "%\n";
		std::cout << "\twindowSz\t=\t" << windowSz << "\n";
		std::cout << "\tdelta\t=\t" << delta << "\n";
		//std::cout << "\tFDR_corr\t=\t" << FDR_corr << "\n";
		std::cout << "\n--------------------------------------------------------\n\n";
	} // printSpec

}; */

/* record the mapping info for each read:
 * (which_bam_file, which_ref, start_on_ref, end)
 * groupId specifies which variant group this read belongs to, will be
 * specified later. groupID = 0 meaning this read is consistent with
 * the reference sequence
 */
struct MapRecord {
	int start;
	int end;
	uint8_t fileID;
	uint8_t refID;
	uint8_t groupID;
	bool isfirstMate;
	MapRecord (){ // indicating empty entry
		start = -1;
		isfirstMate = false;
	}
	MapRecord (int s, int e, uint8_t fId, uint8_t refId, bool fm) {
		start = s; end = e; refID = refId;
		fileID = fId; groupID = 0; isfirstMate = fm;
	}
};

/* the mapping location of the read in the second entry, if not empty,
 * should have coordinate greater equal than the first entry */
typedef std::pair<MapRecord, MapRecord> MapEntry;



/* for each refName: there exists a corresponding refID in
 * each bam file [fID]
 */
struct RefContent{
	std::string seq;
	int length;
	//int gID; 	// global id assigned to this reference
	std::vector<ipair_t> bfID_refID;
	RefContent (int l, int bamfileID, int refID) {
		//seq = "";
		length = l;
		bfID_refID = std::vector<ipair_t> (1, ipair_t (bamfileID, refID));
	};
	RefContent () {};
};

typedef std::map<std::string, RefContent> Ref;

/* information about jeb bucket */
class Bucket {
private:
	int err_S, sum_S, num_bk_S; // summation of total errors, summation of total cnts, num buckets
	int err_L, sum_L, num_bk_L;
	double bayes_prior_L, bayes_prior_S;
public:
	Bucket () {
		err_S = sum_S = num_bk_S = 0;
		err_L = sum_L = num_bk_L = 0;
	}
	void add_err_S (int num) { err_S += num; }
	void add_err_L (int num) { err_L += num; }
	void add_sum_S (int num) { sum_S += num; }
	void add_sum_L (int num) { sum_L += num; }
	void decre_err_S (int num) { err_S -= num; }
	void decre_err_L (int num) { err_L -= num; }
	void decre_sum_S (int num) { sum_S -= num; }
	void decre_sum_L (int num) { sum_L -= num; }

	void incre_bknum_L () { ++ num_bk_L; }
	void incre_bknum_S () { ++ num_bk_S; }
 	void calculate_prior () {
 		bayes_prior_L = (1.0 * err_L) / num_bk_L;
 		bayes_prior_S = (1.0 * err_S) / num_bk_S;
 	};

 	double get_prior_L () const { return bayes_prior_L; }
 	double get_prior_S () const { return bayes_prior_S; }

 	int average_L () const {	return sum_L/num_bk_L; 	}
 	int average_S () const {	return sum_S/num_bk_S;	}
 	double pe_L () { return 1.0 * err_L/sum_L; }
 	double pe_S () { return 1.0 * err_S/sum_S; }
 	void print () {
 		std::cout << "\t LV pe = " << err_L << "/" << sum_L
 				<< " = " << pe_L()  << "\n";
 		std::cout << "\t avg bucket size: " << average_L() << "\n";
 		std::cout << "\t SNPV pe = " << err_S << "/" << sum_S
 				<< " = " << pe_S() << "\n";
 		std::cout << "\t avg bucket size: " << average_S() << "\n";
 	}
};

struct GlobalParam{
	int avgFragSz;	// average fragment size
	int stdFragSz;	// standard deviation of fragment sz
	int maxRL;             // max read length
	uint32_t bonferroni;
	//int constant_pe;
	//double pe;
	//double lenpoly_pe;
	Bucket bkinfo;
	strvec_t bam_filelist;
	std::vector<strvec_t> alnfiles; // alnfiles[i] --> the ith ref genome
	std::string ebfile;
	/* names of ref sequences ordered in bam file for read alignment */
	Ref reference;

	/* quality score --> quantile map */
	imap_t qq;
	/* readName -> mapEntry */
	std::map<std::string, MapEntry> alnMap;

	/* determined by (pair, cycle, dint, quantile)
	 * pair -- 0: first pair, 1: second pair
	 * cycle -- sequencing cycle
	 * dint -- di-nucleotide content (see xny::EnumParser)
	 * quantile -- quality score quantile
	 * errBuckets[i].first -- # cnts where read_base != ref_base
	 * errBuckets[i].second -- total cnts
	 */
	//std::vector<ipair_t> errBucket;

	GlobalParam () { // initialize
		bonferroni = 0;
		//constant_pe = 10000; // to increase pe this fold, increase precision
		                     // when comparing small probabilities
		avgFragSz = stdFragSz = 0;
		maxRL = 0;
	}

};

#endif /* PARAMETER_H_ */

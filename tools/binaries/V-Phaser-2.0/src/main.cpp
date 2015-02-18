//========================================================================
// Project     : VariantCaller
// Name        : main.cpp
// Author      : Xiao Yang
// Created on  : Mar 20, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
// Assumptions :
//		  1) the number of reference sequences, input bam files,
//		      and variants at each loci is no more than 256
//		  2) each read is uniquely mapped
//  		  3) mapping information is not duplicated, i.e. when a reference
//			 sequence ref_i is present in more than 1 bam file, each read
//			 is aligned to ref_i in only one of these files.
//		  4) PCR duplicated read pairs have been marked
//========================================================================

#include <iostream>
#include "Parameter.h"
#include "xutil.h"
#include "bam_manip.h"
#include "xny/file_manip.hpp"
#include "format.h"
#include "em.h"

int main (int argc, char** argv){

	double timing = get_time();
	double start_time = timing;

	Parameter myPara (argc, argv);
	GlobalParam gParam;

	// -------------------------------------------------------------------

	//std::cout << "Get input bam file list...\n\n";
	//xny::getInputFilelist (gParam.bam_filelist, "bam", myPara.iDirNm);

	gParam.bam_filelist.push_back(myPara.iDirNm);

	if (gParam.bam_filelist.size() == 0) abording ("no .bam file found");
	else {
		int sz = gParam.bam_filelist.size();
		std::cout << "\n\t" << sz  << " bam file(s) found: \n";
		print_strvec ("\t\t", gParam.bam_filelist);
	}

	// -------------------------------------------------------------------

	std::cout << "\nParse bam header: get refSeq info & sanity check\n\n";
	parse_bam_header (gParam, myPara.pSample);

	// -------------------------------------------------------------------

	std::cout << "\nGet maxQ, minQ, maxReadLen, avgFragSz, "
			"stdFragSz from bam files ...\n\n";

	int min_qual, max_qual;
	sampling (min_qual, max_qual, gParam.maxRL, gParam.avgFragSz,
			gParam.stdFragSz, gParam.bam_filelist, myPara.pSample);
	if (max_qual - min_qual + 1 < myPara.var_qt) {
		myPara.var_qt = max_qual - min_qual + 1;
		std::cout << "\tqQuantile is reset to " << myPara.var_qt << "\n";
	}

	// -------------------------------------------------------------------
	std::cout << "\nGenerate qual -> quantile map ... \n\n";

	qqMap (gParam.qq, min_qual, max_qual, myPara.var_qt);

	// -------------------------------------------------------------------
	std::cout << "\nSet up paired read map arrays ... \n\n";
	set_rmap_array (gParam, myPara);
	//set_arrays (gParam, qq, myPara);

	// remeasure fragment size based on paired read alignment as the
	// BAM file did not properly record such information
	if (gParam.avgFragSz == 0) {
		ivec_t fragSz;
		std::map<std::string, MapEntry>::iterator it_am = gParam.alnMap.begin();
		for (; it_am != gParam.alnMap.end(); ++ it_am) {
			if (it_am->second.second.start >= it_am->second.first.start){
				fragSz.push_back(it_am->second.second.start -
						it_am->second.first.start);
			}
		}
		/* calculate mean frag sz and std */
		int sample_sz = fragSz.size();
		uint64_t sum = 0;
		for (int i = 0; i < sample_sz; ++ i) sum += fragSz[i];
		if (sum > 0 && sample_sz > 0) gParam.avgFragSz = sum/sample_sz;
		sum = 0;
		for (int i = 0; i < (int) fragSz.size(); ++ i) {
			sum += (gParam.avgFragSz - fragSz[i])*(gParam.avgFragSz - fragSz[i]);
		}
		if (sample_sz > 1) gParam.stdFragSz = std::sqrt(sum/sample_sz - 1);
		std::cout << "\trecalculated frag size and std" << gParam.avgFragSz
				<< ", " << gParam.stdFragSz << "\n";
	}

	// -------------------------------------------------------------------
	std::cout << "\nPrepare aln columns file...\n\n";
	prep_aln_file (gParam, myPara);

	// -------------------------------------------------------------------
	std::cout << "\nInference ...\n\n";

	EM (gParam, myPara);

	std::cout << "-----------------------------------------------------------\n\n";

    print_time("Whole program takes \t", start_time);
	std::cout << "DONE!\n";
	return (EXIT_SUCCESS);
}

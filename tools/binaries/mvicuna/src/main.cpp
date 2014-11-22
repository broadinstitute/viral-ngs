// Project     : M-Vicuna (Modularized-Vicuna)
// Name        : main.cpp
// Author      : Xiao Yang
// Created on  : May 13, 2013
// Version     : 1.0
// Copyright   : The Broad Institute
//  				 SOFTWARE COPYRIGHT NOTICE AGREEMENT
// 				 This software and its documentation are copyright (2013)
//				 by the Broad Institute. All rights are reserved.
//
// 				 This software is supplied without any warranty or
//				 guaranteed support whatsoever. The Broad Institute cannot
//				 be responsible for its use,	misuse, or functionality.
// Description :

#include <iostream>
#include "Parameter.h"
#include "xutil.h"
#include "DuplRm.h"
#include "MergeReadPair.h"
#include "Trim.h"
#include "SeqFrqEst.h"

#include "xny/seq_cmp.hpp"
#include "jaz/fastx_iterator.hpp"
#include "xny/sketch.hpp"
#include "jaz/hash.hpp"

int main (int argc, char** argv){
	/*strvec_t input {"output/a.fq", "output/b.fq"};
	std::string output = "output/test.fq";
	xny::deletefile(output);
	xny::append2file (output, input);
	return 0; */
	/*
	std::vector<std::string> myS = {"GCTGAAAGTTGCGAGCTGAC",
									"GCTGAAAGTTGCGAGCTGTC"};
	jaz::murmur264 myhash (546);

	for (auto s0 : myS) {

		xny::sketch_list mylist (7, true);
		std::vector<sketch_t> list = mylist (s0, myhash);
		std::cout << s0 << "\n";
		for (auto x: list) {
			std::cout << x.first << ", " << x.second << "\n";
		}
		std::cout << std::endl;
		xny::super_sketch ss (3);
		std::sort(list.begin(), list.end(), xny::cmp_sketch());

		for (auto x: list) {
			std::cout << x.first << ", " << x.second << "\n";
		}
		std::cout << std::endl;

		sketch_t my_ss;
		if (ss (my_ss, list, myhash)) {
			std::cout << "\n" << my_ss.first << ", " << my_ss.second << "\n";
		}
	}
	exit(1);
	*/
	/*{// to test alignment
		//scores (1, -2, -6)
		// >sequence one
		// AGTGCTGAAAGTTGCGCCAGC-TGAC-
		//>sequence two
		// AGTGCTGAA-GTT-CGCCAGCTTGACG
		std::string s0 = "GCTGAAAGTTGCGAGCTGAC",
					s1 = "GCTGAAGTTCGAGCTTGACG";
		//s0 = "cctaaa";
		//s1 = "ccctaaagg";

		s0 = "agttgg";
		s1 = "gtt";
		bio::global_alignment galn (1, -3, -5, -2);
		galn.set_alignment_type(1);
		galn (s0, s1);
		std::cout << galn.path() << "\n";
	}*/
	//GCTGAACAGTTGCGAGCTTGAC-
	//GCTGAAC-GTT-CGAGCTTGACG";

	/*{ // to test kmer_anchor_aln
		std::string s0 = "ATTTGACTGAACAGTTGCGAGCTTGACTT",
					s1 =    "GCACTGAACGTTCGAGCTTGACGTCCA";

		xny::kmer_anchor_aln aln (21, 90, 2);

		std::tuple<coord_t, std::string >
			rslt = aln (s0, s1);
		std::cout << "final result:\n" << std::get<1> (rslt) << "\n";
		exit(1);
	}*/
	//

	double timing = get_time();
	double start_time = timing;
	Parameter myPara (argc, argv);

	strset_t inputfiles (myPara.ipfq.begin(), myPara.ipfq.end());
	inputfiles.insert(myPara.isfq.begin(), myPara.isfq.end());
	inputfiles.insert(myPara.ifa.begin(), myPara.ifa.end());

	strset_t intermediate_files;

 	omp_set_num_threads(myPara.pthreads);

 	xny::low_complexity lc (myPara.lc_n, myPara.lc_mono, myPara.lc_di);

	for (auto& task: myPara.tasks) {

		switch (task) {
			/** Task: de novo duplicate removal for each pair of input fq pairs
			 * 	-- remove duplicated fragments
			 * 	-- remove low complexity fragments
			 *
			 *	Input: a list of paired fastq files
			 *	Output: a list of paired fastq files or 2 fastq files (as a pair)
			 */
		case 0: // duprm & low complexity frag removal
		{
			if(!myPara.silent) {
				std::cout << "Duplicate removal...\n";
				std::cout << "\tinput: ";
				for (auto& x : myPara.ipfq) std::cout << "\n\t\t" << x;
				std::cout << "\n\n";
			}

			duplicate_removal (myPara.ipfq, myPara.drm, myPara.w,
				myPara.w2, lc, myPara.batch, myPara.silent);

			// update the input paired fastq for next stage !
			myPara.ipfq = myPara.drm.op;

			// update intermediate_files that may be removed
			intermediate_files.insert(myPara.drm.op.begin(),	myPara.drm.op.end());

			if(!myPara.silent) {
				std::cout << "\n\toutput:";
				for (auto& x: myPara.drm.op) std::cout << "\n\t\t" << x;
				std::cout << "\n\n";

				print_time("duplicate removal complete !\t", timing);
			}

			break;
		}
		case 1: // trim
		{
			/** Task:
			 * 	 -- trim known primers
			 * 	 -- low quality bases
			 * 	 -- low complexity reads
			 * 	Input:
			 * 	 -- a list of paired fastq files to be trimmed;
			 * 	 -- a list of unpaired fastq files to be trimmed;
			 *	 -- a fasta file contains primers to apply trimming
			 *	Output:
			 *	 -- a list of output paired fastq files (>= 2)
			 *	 -- a list of output singleton fastq files (>= 1)
			 */
			if(!myPara.silent) {
				std::cout << "Trimming ...\n";
				std::cout << "\tinput: ";
				for (auto& x : myPara.ipfq) std::cout << "\n\t\t" << x;
				for (auto& x : myPara.isfq) std::cout << "\n\t\t" << x;
				std::cout << "\n\n";
			}

			trimming (myPara.ipfq, myPara.isfq, myPara.trm, lc,
					myPara.batch, myPara.silent);

			// update the input paired fastq files for next stage
			myPara.ipfq = myPara.trm.op;
			myPara.isfq.insert (myPara.isfq.end(), myPara.trm.os.begin(),
					myPara.trm.os.end());

			// update intermediate_files that may be removed
			intermediate_files.insert(myPara.trm.os.begin(),	myPara.trm.os.end());
			intermediate_files.insert(myPara.trm.op.begin(),	myPara.trm.op.end());

			if(!myPara.silent) {
				std::cout << "\n\toutput:";
				for (auto& x: myPara.trm.op) std::cout << "\n\t\t" << x;
				for (auto& x: myPara.trm.os) std::cout << "\n\t\t" << x;
				std::cout << "\n\n";

				print_time("trimming complete !\t", timing);
			}
			break;
		}
		case  2: // PairedReadMerge
			/** Task: Merge paired-end read.
			 * Input: a list of paired fastq files
			 * Output:
			 * -- a fasta file recording merged pairs
			 * -- a pair of fastq files recording unmerged pairs
			 * Note: the second pair is reverse complemented wrt the input
			 */
			if(!myPara.silent) {
				std::cout << "Merge paired-read ...\n";
				std::cout << "\tinput: ";
				for (auto& x : myPara.ipfq) std::cout << "\n\t\t" << x;
				std::cout << "\n\n";
			}

			merge_paired_read (myPara.ipfq, myPara.prm.op,
							   myPara.prm.os, myPara.batch);

			// update the input paired fastq and fasta for the next stage
			myPara.ipfq = strvec_t { myPara.prm.op.first, myPara.prm.op.second };
			myPara.isfq.push_back(myPara.prm.os);

			// update intermediate_files that may be removed
			intermediate_files.insert(myPara.prm.op.first);
			intermediate_files.insert(myPara.prm.op.second);
			intermediate_files.insert(myPara.prm.os);

			if(!myPara.silent) {
				std::cout << "\n\toutput:";
				std::cout << "\n\t\t" << myPara.prm.op.first;
				std::cout << "\n\t\t" << myPara.prm.op.second;
				std::cout << "\n\t\t" << myPara.prm.os;
				std::cout << "\n\n";

				print_time("merge paired-read complete!\t", timing);
			}
			break;
		case 3: // sequence frequency estimate -- statistical repeat masking !!
			if (!myPara.silent) {
				std::cout << "Sequence frequency estimate ...\n";
				std::cout << "\tinput: ";
				for (auto& x : myPara.ipfq) std::cout << "\n\t\t" << x;
				for (auto& x : myPara.isfq) std::cout << "\n\t\t" << x;
				std::cout << "\n\n";
			}

			estSeqFrq (myPara.ipfq, myPara.isfq, myPara.fe_k, myPara.batch, myPara.silent);
			if(!myPara.silent) print_time("seq frq estimate complete !\t", timing);
			break;
		default:
			break;

		} // switch (task)
	} // for (auto& x: myPara.tasks) {

	//---------- consolidate final output files -------------------

	if (!myPara.silent) std::cout << "\nconsolidate output files...\n";

	// identify paired fq and singleton fq files to consolidate
	strvec_t pfq1, pfq2, sfq;
	for (unsigned int i = 0; i < myPara.ipfq.size(); i += 2) {
		if (!inputfiles.count(myPara.ipfq[i])) {
			pfq1.push_back(myPara.ipfq[i]);
			pfq2.push_back(myPara.ipfq[i + 1]);
			intermediate_files.erase(myPara.ipfq[i]);
			intermediate_files.erase(myPara.ipfq[i + 1]);
		}
	}
	for (auto& x: myPara.isfq) { // singleton fq files
		if (!inputfiles.count(x)) {
			sfq.push_back(x);
			intermediate_files.erase(x);
		}
	}

    // consolidate
	if (pfq1.size() == 1 && !myPara.noclean) {
		std::cout << "\trenaming " << pfq1[0] << " to " << myPara.opfq[0] << "\n";
		std::rename(pfq1.front().c_str(), myPara.opfq.front().c_str());
		std::cout << "\trenaming " << pfq2[0] << " to " << myPara.opfq[1] << "\n";
		std::rename(pfq2.front().c_str(), myPara.opfq.back().c_str());
		pfq1.clear();
		pfq2.clear();
	} else { // merge to 2 output paired fq files
		xny::deletefile(myPara.opfq.front());
		xny::deletefile(myPara.opfq.back());
		xny::append2file(myPara.opfq.front(), pfq1);
		xny::append2file(myPara.opfq.back(), pfq2);
	}

	if (sfq.size() == 1 && !myPara.noclean) { // change name
		std::cout << "\trenaming " << sfq[0] << " to " << myPara.osfq << "\n";
		std::rename(sfq[0].c_str(), myPara.osfq.c_str());
		sfq.clear();
	} else { // merge to an output file
		xny::deletefile(myPara.osfq);
		xny::append2file(myPara.osfq, sfq);
	}

	// --------------- remove intermediate files -----------------------
	if (!myPara.noclean) {
		for (auto& x: pfq1) xny::deletefile (x);
		for (auto& x: pfq2) xny::deletefile (x);
		for (auto& x: sfq) xny::deletefile (x);
		for (auto& x: intermediate_files) xny::deletefile (x);
	}

	if (!myPara.silent) {
		print_time("Whole program takes \t", start_time);
		std::cout << "DONE!\n";
	}
	return (EXIT_SUCCESS);
}







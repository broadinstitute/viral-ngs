  /// Project     : M-Vicuna (Modularized-Vicuna)
/// Name        : Parameter.h
/// Author      : Xiao Yang
/// Created on  : May 13, 2013
/// Version     : 1.0
/// Copyright   : The Broad Institute
///  				 SOFTWARE COPYRIGHT NOTICE AGREEMENT
/// 				 This software and its documentation are copyright (2013)
///				 by the Broad Institute. All rights are reserved.
///
/// 				 This software is supplied without any warranty or
///				 guaranteed support whatsoever. The Broad Institute cannot
///				 be responsible for its use,	misuse, or functionality.
/// Description :


#ifndef PARAMETER_H_
#define PARAMETER_H_

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "xutil.h"
#include "jaz/string_add.hpp"

/* task: de novo duplicate removal */
struct drm_t {
 	strvec_t op; // duplicate removed output paired fq files
	int perc_sim;
//	int max_mismatch;
};

/* task: Trim */
struct trm_t {
	std::string vecfa;
	strvec_t op;  // output unmerged pairs in fastq format
	strvec_t os;  // output merged seq in single-end fastq format
	int min_match; // minimum match b/t read and vector to apply trim
	int min_rlen; // minimum remaining read length
	int min_qual;
};

/* task: paired read merging */
struct prm_t {
	strpair_t op; // output unmerged pairs in fastq format
	std::string os; // output merged seq in single-end fastq format
};

/* input parameters */
class Parameter{

public:
    enum TaskList { DupRm, Trim, PairedReadMerge, SFrqEst };
    strset_t task_set_ { "DupRm", "Trim", "PairedReadMerge", "SFrqEst"};

	strvec_t ipfq; // input paired-read comma separated fastq files;
				  // pair1 and pair2 are adjacent to each other
	strvec_t isfq; // comma separated single end fastq files
	strvec_t ifa;  //  comma separated fasta files
	std::vector<unsigned char> tasks; // task list

	strvec_t opfq;	// final output fastq paired files
	std::string osfq; // final output single fastq file

	// ----- general parameters to be used in multiple tasks ----------
	int lc_n;   // low complexity seq max percentage of ambiguous bases
	int lc_mono; // low complexity seq max percentage of mono bases
	int lc_di; // low complexity seq max percentage of di bases
	int w, w2; // sketching window 1 and 2

	//-------- general performance tuning parameters ------------------
	int batch; // number of of single read/contigs; paired reads to be
			   // read into the memory
	int pthreads;
	bool silent, noclean;

	drm_t drm; // duplicate removal
	trm_t trm; // trim
	prm_t prm; // paired read merging

	int fe_k; // sequence frequency estimate
	/* task: high stringent clustering */
	//bool	 t_hsc; // high stringent read clustering


	Parameter (int argc, char** argv): argnum(argc), arg(argv){
		init ();
		for (int i = 1; i < argnum; i += 2) {
			std::string option = argv[i];
			if (option.compare("-h") == 0 ||	option.compare("--h") == 0) {
				printUsage (argv[0]);
			} else if (option.compare("-ipfq") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string tmp = argv[i+1];
				split(',', tmp, std::back_inserter(ipfq));
				if (ipfq.size() % 2 != 0) {
					abording("pfq should contain even number of files\n type -h to get options");
				}
				//print_1dvec (pfq, &std::cout); // debug
			} else if (option.compare("-isfq") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string tmp = argv[i+1];
				split(',', tmp, std::back_inserter(isfq));
				//print_1dvec (sfq, &std::cout); // debug
			} else if (option.compare("-fa") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string tmp = argv[i+1];
				split(',', tmp, std::back_inserter(ifa));
			} else if (option.compare("-opfq") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string tmp = argv[i+1];
				split(',', tmp, std::back_inserter(opfq));
				if (opfq.size() != 2) {
					abording("-opfq should be 2 paired fastq files\n type -h to get options");
				}
			} else if (option.compare("-osfq") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				osfq = argv[i+1];
			} else if (option.compare("-tasks") == 0) {
				strvec_t tasklist;
				if (argc < i + 2) {
					printUsage (argv[0]);
				}
				std::string tmp = argv[i+1];
				split(',', tmp, std::back_inserter(tasklist));
				tasks.resize(tasklist.size());
				int idx = 0;
				for (auto &x : tasklist) {
					if (x.compare("Trim") == 0) {
						tasks[idx] = TaskList::Trim;
					} else if (x.compare ("DupRm") == 0) {
						tasks[idx] = TaskList::DupRm;
					} else if (x.compare ("PairedReadMerge") == 0) {
						tasks[idx] = TaskList::PairedReadMerge;
					} else if (x.compare ("SFrqEst") == 0){
						tasks[idx] = TaskList::SFrqEst;
					} else {
						std::cout << "task: " << x << "unrecognized\n";
						exit(1);
					}
					++ idx;
				} // for
			} else if (option.compare("-batch") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				batch = atoi (argv[i+1]);
				if (batch < 10000) {
					warning ("-batch is too small, batch is reset to 10k");
					batch = 10000;
				}
			} else if (option.compare("-w") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				w = atoi (argv[i+1]);
			} else if (option.compare("-w2") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				w2 = atoi (argv[i+1]);
			} else if (option.compare("-lc_n") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				lc_n = atoi (argv[i+1]);
				if (lc_n < 1 || lc_n > 100) {
					warning ("-lc_n in [1, 100],  reset to 30");
					lc_n = 30;
				}
			} else if (option.compare("-lc_mono") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				lc_mono = atoi (argv[i+1]);
				if (lc_mono < 1 || lc_mono > 100) {
					warning ("-lc_mono in [1, 100] reset to 50");
					lc_mono = 50;
				}
			} else if (option.compare("-lc_di") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				lc_di = atoi (argv[i+1]);
				if (lc_di < 1 || lc_di > 100) {
					warning ("-lc_di in [1, 100], reset to 80");
					lc_di = 80;
				}
			} else if (option.compare("-pthreads") == 0){
				if (argc < i + 2) printUsage (argv[0]);
				pthreads = atoi (argv[i+1]);
			} else if (option.compare("-silent") == 0) {
				silent = true;
				-- i;
			} else if (option.compare("-noclean") == 0) {
				noclean = true;
				-- i;
			}

			//---------------- duplicate removal task ------------------
			else if (option.compare ("-drm_op") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string tmp = argv[i + 1];
				split (',', tmp, std::back_inserter(drm.op));
		    } else if (option.compare("-drm_perc_sim") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
		    		drm.perc_sim = atoi (argv[i+1]);
		    		if (drm.perc_sim < 95) {
		    			warning ("-drm_perc_sim is too small, it is reset to 95");
		    			drm.perc_sim = 95;
		    		}
		    }/* else if (option.compare ("-drm_max_mismatch") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
		    		drm.max_mismatch = atoi (argv[i + 1]);
		    }*/

			//-------------- paired read removal task ------------------
			else if (option.compare("-prm_op") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string tmp = argv[i+1];
				strvec_t tmp_splits;
				split(',', tmp, std::back_inserter(tmp_splits));
				if (tmp_splits.size() != 2) {
					abording ("-prm_op requires 2 output fq files\n type -h to get options");
				} else prm.op = strpair_t(tmp_splits[0], tmp_splits[1]);
			} else if (option.compare("-prm_os") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				prm.os = argv[i+1];
			}

			//-------------- Trim task ------------------
			else if (option.compare("-trm_vecfa") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				trm.vecfa = argv[i+1];
			} else if (option.compare("-trm_op") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string tmp = argv[i+1];
				split(',', tmp, std::back_inserter(trm.op));
			} else if (option.compare("-trm_os") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				std::string tmp = argv[i+1];
				split(',', tmp, std::back_inserter(trm.os));
			} else if (option.compare("-trm_min_match") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				trm.min_match = atoi (argv[i + 1]);
			} else if (option.compare("-trm_min_rlen") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				trm.min_rlen = atoi (argv[i + 1]);
			} else if (option.compare ("-trm_q") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				trm.min_qual = atoi (argv[i + 1]);
			}

			//-------------- SFrqEst task ------------------
			else if (option.compare("-fe_k") == 0) {
				if (argc < i + 2) printUsage (argv[0]);
				fe_k = atoi (argv[i+1]);
			} else {
				std::cout << "flag: " << option << " is unknown\n"
						<< " type -h to get options\n";
				exit(1);
			}
		} // for (int i = 1; i < argnum; i += 2) {

		printSpec(argv[0]);
	} // Parameter

private:
	int argnum;
	char** arg;

	void init () {
		silent = false;
		noclean = false;
		batch = 500000;
		pthreads = 8;
		lc_n = 30;
		lc_mono = 50;
		lc_di = 80;
		w = 17;
		w2 = 5;
		tasks = {0,1,2,3};

		drm.perc_sim = 98;
		//drm.max_mismatch = 5;

		trm.min_match = 13;
		trm.min_rlen = 70;
		trm.min_qual = 2;

		fe_k = 14;
	}

	void printUsage(char* exe) {
		std::cout << "\n--------------------------------------------------------\n";
		std::cout << "Parameters\n";
		std::cout << "-ipfq: comma separated input paired fastq files; the ith and (i+1)th files form a pair (i is an odd number)\n";
		std::cout << "-isfq: comma separated input single end fastq files\n";
		std::cout << "-fa: comma separated input single end fasta files\n";
		std::cout << "-opfq: comma separated final 2 output fastq paired files\n";
		std::cout << "-osfq: final output singleton fastq file\n";
		std::cout << "-batch: default 500000; number of sequence (pairs) to be loaded in the memory (>=10000)";
		std::cout << "-pthreads: default 8; number of cores to use\n";
		std::cout << "-w, -w2: default 17, 5; sketching window sizes\n";
		std::cout << "-lc_n, -lc_mono, -lc_di: default 30, 50, 80; defining low complexity sequence"
				  << "\n\tmax percentage of ambiguous bases, mono nucleotides, and dinucleotides\n";
		std::cout << "-tasks: default DupRm,Trim,PairedReadMerge,SFrqEst; "
				"\n\ta list of comma separated tasks {DupRm, Trim, PairedReadMerge, SFrqEst}\n";
		std::cout << "-silent: default false; no screen print-out\n";
		std::cout << "-noclean: default false; do not remove intermediate files\n";
		std::cout << std::endl;

		/* duplicate & low complexity fragment removal */
		std::cout << "TASK: DupRm\n";
		std::cout << "-drm_op: comma separated output paired fq files post dup rm\n";
		std::cout << "-drm_perc_sim: default 98; percent similarity\n";
		std::cout << "-drm_max_mismatch: default 5; max mismatches allowed\n";
		std::cout << std::endl;

		/* paired read merging */
		std::cout << "TASK: PairedReadMerge\n";
		std::cout << "-prm_op: 2 comma separated output unmerged fq files\n";
		std::cout << "-prm_os: merged single-end fq file\n";
		std::cout << std::endl;

		/* trimming */
		std::cout << "TASK: Trim\n";
		std::cout << "-trm_vecfa: input fasta file storing vectors/primers\n";
		std::cout << "-trm_op: comma separated output fq paired files\n";
		std::cout << "-trm_os: merged single-end fq files\n";
		std::cout << "-trm_min_match: default 13; min match b/t vector and a read to be trimmed\n";
		std::cout << "-trm_min_rlen: default 70; min read length post-trimming\n";
		std::cout << "-trm_q: default 2 (ASCII 35 #); min phred score (ASCII >= 33)\n";
		std::cout << std::endl;

		/* sequence frequency estimate */
		std::cout << "TASK: SFrqEst -- sequence frequency estimation\n";
		std::cout << "-fe_k: default 14 (<= 16); substring length to calibrate\n";
		std::cout << std::endl;

		std::cout << "----------------------------------------------------------\n";

		exit(1);
	}

	void print_file_list (const std::string& flag, const strvec_t& flist) {
		if (flist.size()) {
			std::cout << flag;
			for (unsigned int i = 0; i < flist.size(); ++ i) {
				std::cout << flist[i];
				if (i!= flist.size() - 1) std::cout << ",";
			}
		}
	}

	void printSpec(char* exe) {
		if (!silent) {
			std::string header = "\n\t";

			std::string mytask;
			std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Running command: \n" << exe;

			print_file_list (header + " -ipfq ", ipfq);
			print_file_list (header + " -isfq ", isfq);
			print_file_list (header + " -fa ", ifa);
			if (opfq.size() != 2) abording ("-opfq not not specified\n type -h to get options");
			print_file_list (header + " -opfq ", opfq);
			if (osfq.empty()) abording ("-osfq not not specified\n type -h to get options");
			std::cout << header << " -osfq " << osfq;
			std::cout << header << " -batch " << batch << header << " -pthreads "
					<< pthreads << header << " -w " << w << header << " -w2 " << w2
					<< header << " -lc_n " << lc_n << header << " -lc_mono " << lc_mono
					<< header << " -lc_di " << lc_di;

			for (auto& task: tasks) {

				switch (task) {
				case DupRm:
					if (mytask.empty()) 	mytask += "DupRm";
					else mytask += ",DupRm";
					if (drm.op.size() == 0) {
						abording ("Task DupRm: -drm_op not specified\n type -h to get options");
					}
					print_file_list (header + " -drm_op ", drm.op);
					std::cout << header << " -drm_perc_sim " << drm.perc_sim;
					//std::cout << header << " -drm_max_mismatch " << drm.max_mismatch;
					break;
				case PairedReadMerge:
					if (mytask.empty()) 	mytask += "PairedReadMerge";
					else mytask += ",PairedReadMerge";

					if (prm.op.first.empty() || prm.os.empty()) {
						std::cout << std::endl;
						abording ("Task PairedReadMerge: -prm_op or -prm_os is not specified\n type -h to get options");
					}
					std::cout << header << " -prm_op " << prm.op.first << "," << prm.op.second;
					std::cout << header << " -prm_os " << prm.os;
					break;
				case Trim:
					if (mytask.empty()) 	mytask += "Trim";
					else mytask += ",Trim";
					if (trm.op.size() == 0 || trm.os.size() == 0) {
						std::cout << std::endl;
						abording ("Task Trim: -trm_op or -trm_os is not specified\n type -h to get options");
					}
					print_file_list (header + " -trm_op ", trm.op);
					print_file_list (header + " -trm_os ", trm.os);
					if (!trm.vecfa.empty()) std::cout << "\n\t -trm_vecfa " << trm.vecfa;
					std::cout << header << " -trm_min_match " << trm.min_match;
					std::cout << header << " -trm_min_rlen " << trm.min_rlen;
					std::cout << header << " -trm_q " << trm.min_qual;
					break;
				case SFrqEst:
					if (mytask.empty()) 	mytask += "SFrqEst";
					else mytask += ",SFrqEst";
					std::cout << header << " -fe_k " << fe_k;
					break;
				}
			}
			std::cout << header << " -tasks " << mytask;
			std::cout << "\n--------------------------------------------------------\n\n";
		}
	} // printSpec
};


#endif /* PARAMETER_H_ */

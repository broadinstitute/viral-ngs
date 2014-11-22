//========================================================================
// Project     : M-Vicuna
// Name        : Trim.h
// Author      : Xiao Yang
// Created on  : Jul 10, 2013
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
//========================================================================


#ifndef TRIM_H_
#define TRIM_H_

#include "xutil.h"
#include "Parameter.h"
#include "xny/seq_cmp.hpp"
#include "xny/file_manip.hpp"
#include "jaz/fastx_iterator.hpp"
#include "ReadBioFile.h"


// kmer location on string s (s_ID, s_startpos, is_fwdstrd)
typedef std::tuple<int, int, bool> kloc_t;
// binary_kmer -> list [ kloc_t ]
typedef std::map<uint32_t, std::vector<kloc_t> > kindex_t;

void trimming (const strvec_t& ipfq, const strvec_t& isfq,
	const trm_t& trm, xny::low_complexity& lc, int batch, bool silent);

void process_vector_file (strvec_t& vectors, kindex_t& kindex,
		const std::string& vecfa, int k);

void trim_pfq (const std::string& ifq, const std::string& ifq2,
	const strvec_t& vectors, const kindex_t& kindex, std::ofstream& ofhfq,
	std::ofstream& ofhfq2, std::ofstream& ofhfa, const trm_t& trm,
	xny::low_complexity& lc, int batch);

int apply_trimming (std::vector<fqtuple_t>& seq, const strvec_t& vectors,
		const kindex_t& kindex, xny::low_complexity& lc, const trm_t& trm);

bool trim_primer (std::string& rSeq, std::string& qual,
		const std::vector<kloc_t>& kloclist,
		const strvec_t& vectors,	int rPos, const trm_t& trm);

bool trim_lc_lq (std::string& rSeq, std::string& qual, int min_qual,
		int min_len, xny::low_complexity& lc);

void trim_sfq (std::ofstream& ofhsfq, const std::string& ifq,
		const trm_t& trm, xny::low_complexity& lc, int batch);


#endif /* TRIM_H_ */

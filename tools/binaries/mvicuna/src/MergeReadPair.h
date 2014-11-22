//========================================================================
// Project     : M-Vicuna
// Name        : MergeReadPair.h
// Author      : Xiao Yang
// Created on  : Jun 3, 2013
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


#ifndef MERGEREADPAIR_H_
#define MERGEREADPAIR_H_

#include "xutil.h"
#include "jaz/fastx_iterator.hpp"
#include "xny/file_manip.hpp"
#include "xny/seq_cmp.hpp"

void merge_paired_read (const strvec_t& ifq, const strpair_t& ofq,
		const std::string& ofa, int batch);


int apply_merging (std::ofstream& fa, std::ofstream& fq, std::ofstream& fq2,
	std::vector<fqtuple_t>& seq);

#endif /* MERGEREADPAIR_H_ */

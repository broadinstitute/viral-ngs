//========================================================================
// Project     : M-Vicuna
// Name        : SeqFrqEst.h
// Author      : Xiao Yang
// Created on  : Aug 6, 2013
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


#ifndef SEQFRQEST_H_
#define SEQFRQEST_H_

#include "xutil.h"
#include "ReadBioFile.h"
#include "jaz/fastx_iterator.hpp"
#include "xny/seq_manip.hpp"
#include "xny/file_manip.hpp"

void estSeqFrq (const strvec_t& ipfqlist, const strvec_t& isfqlist,
		int k, int batch, bool silent);

void obtain_kfrq (umap_t& kcnt, int k, int batch, const strvec_t& files, bool silent);
void kmer_cnt_in_seqs (umap_t& kcnt,	const strvec_t& seqs, int k);


void obtain_seqfrq (ivec_t& frq, const umap_t& kcnt,	int k,
		int batch, const strvec_t& files, bool silent) ;
void seq_freq (ivec_t& frq, const umap_t& kcnt, int k, const strvec_t& seqs);

#endif /* SEQFRQEST_H_ */

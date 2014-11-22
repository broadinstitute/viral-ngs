//========================================================================
// Project     : M-Vicuna
// Name        : ReadBioFile.cpp
// Author      : Xiao Yang
// Created on  : Jun 4, 2013
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

#include "ReadBioFile.h"

/**	Retrieve num records (header, read, qual) from fastq file pointed by
 * iter
 */
void add_fq_reads (std::vector<fqtuple_t>& seq, int num,
	bio::fastq_input_iterator<>& iter, bio::fastq_input_iterator<> end){
	int cnt = 0;
	for (; iter != end; ++ iter) {
		seq.push_back (*iter);
		++ cnt;
		if (cnt >= num) {
			++ iter;
			return;
		}
	}
} // add_fq_reads

void add_fq_reads_only (strvec_t& seq, int num,
	bio::fastq_input_iterator<>& iter, bio::fastq_input_iterator<> end){
	int cnt = 0;
	for (; iter != end; ++ iter) {
		seq.push_back (std::get<1>(*iter));
		++ cnt;
		if (cnt >= num) {
			++ iter;
			return;
		}
	}
}

void add_fa_reads (std::vector<strpair_t>& seq, int num,
	bio::fasta_input_iterator<>& iter, bio::fasta_input_iterator<> end){
	int cnt = 0;
	for (; iter != end; ++ iter) {
		seq.push_back (*iter);
		++ cnt;
		if (cnt >= num) {
			++ iter;
			return;
		}
	}
} //add_fa_reads

void add_fa_reads_only (strvec_t& seq, int num,
	bio::fasta_input_iterator<>& iter, bio::fasta_input_iterator<> end){
	int cnt = 0;
	for (; iter != end; ++ iter) {
		seq.push_back (iter->second);
		++ cnt;
		if (cnt >= num) {
			++ iter;
			return;
		}
	}
} //add_fa_reads



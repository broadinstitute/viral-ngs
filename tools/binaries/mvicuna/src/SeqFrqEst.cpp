//========================================================================
// Project     : M-Vicuna
// Name        : SeqFrqEst.cpp
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


#include "SeqFrqEst.h"

/* @brief	Given a list of paired fastq files [ipfq] and unpaired fastq
 * 	files [isfq], the kmer size [k] estimate the frequency of each input
 * 	read
 */
void estSeqFrq (const strvec_t& ipfq, const strvec_t& isfq,
		int k, int batch, bool silent){

	strvec_t myfiles = ipfq;
	myfiles.insert(myfiles.end(), isfq.begin(), isfq.end());

	if (! silent) std::cout << "\tobtain " << k << "-mer frequency\n";
	umap_t kcnt;
	obtain_kfrq (kcnt, k, batch, myfiles, silent);

	if (! silent) std::cout << "\tnumer of kmers: " << kcnt.size() << "\n";

	ivec_t frq;
	obtain_seqfrq (frq, kcnt, k, batch, myfiles, silent);


	{ // debug
		std::cout << kcnt.size() << " kmers observed\n";
		uvec_t kfrq;
		for (auto & x: kcnt) { kfrq.push_back(x.second); }
		std::sort (kfrq.begin(), kfrq.end());
		for (uvec_t::reverse_iterator it = kfrq.rbegin(); it != kfrq.rend(); ++ it) {
			std::cout << *it << "\t";
		}
		std::cout<< "\n\n";
		std::cout << "# reads " << frq.size() << "\n";
		for (auto & x: frq) std::cout << x << "\n";
	}
} // estSeqFrq

/* @brief	Estimate sequence frequency for all reads in fastq files
 */
void obtain_seqfrq (ivec_t& frq, const umap_t& kcnt,
		int k, int batch, const strvec_t& files, bool silent) {
	strvec_t seqs;
	int cnt = batch, total_reads = 0;
	for (int i = 0; i < (int) files.size(); ++ i) {

		if (!silent) std::cout << "\t\tprocess file: " << files[i] << "\n";

		std::ifstream fh;
		xny::openfile<std::ifstream>(fh, files[i]);

		bio::fastq_input_iterator<> fq(fh), end;

		while (fq != end) {

			add_fq_reads_only (seqs, cnt, fq, end);

			if ((int) seqs.size() >= batch) {
				total_reads += seqs.size();
				seq_freq (frq, kcnt, k, seqs);
				cnt = batch;
				seqs.clear();
			} else { // not enough reads to fill in batch for current file pair
				cnt = batch - seqs.size();
				xny::closefile(fh);
			}
		} // while
	} // for
	total_reads += seqs.size();
	seq_freq (frq, kcnt, k, seqs);
	if (!silent) std::cout << "\t\t" << total_reads << " reads analyzed\n";

} // obtain_seqfrq

/* @brief	Adding frequency of sequences [seqs] to vector [frq]
 */
void seq_freq (ivec_t& frq, const umap_t& kcnt,
		int k, const strvec_t& seqs) {

	int sz = seqs.size();
	ivec_t local_frq (sz, 0);
//	int total_debug = 0;
	#pragma omp parallel for
	for (int i = 0; i < sz; ++ i) {
		int frq_sum = 0;
		uvec_t kIDs;
		xny::get_bitkmer<std::back_insert_iterator<uvec_t>, uint32_t>
				(seqs[i], std::back_inserter(kIDs), k, 0);
		uset_t kID_set (kIDs.begin(), kIDs.end());

		for (auto& id: kID_set) {
			umap_t::const_iterator it = kcnt.find(id);
			if (it != kcnt.end()) frq_sum += it->second;
		}
		if (kID_set.size() != 0) {
			local_frq[i] = frq_sum / kID_set.size();

			{ // debug
				if (local_frq[i] > 1000) {
					std::cout << ">" << i << "\n" <<  seqs[i] << "\n";
					/*
					uvec_t cnts;
					for (auto& id: kIDs) {
						int cnt = kcnt.find(id)->second;
						std::cout << cnt << ", ";
						cnts.push_back(cnt);
					}
					std::cout << "\n";
					std::sort (cnts.begin(), cnts.end());
					for (auto& x : cnts) std::cout << x << ", ";
					std::cout << "\n\n";
					*/
				}
			}
		}
	}
	frq.insert(frq.end(), local_frq.begin(), local_frq.end());
//	std::cout << "total debug = " << total_debug << "\n";
} // seq_freq

/* @brief	Obtain kmer frequency in a given set of fastq files
 */
void obtain_kfrq (umap_t& kcnt, int k, int batch,
		const strvec_t& files, bool silent) {

	strvec_t seqs;
	int cnt = batch, total_reads = 0;
	for (int i = 0; i < (int) files.size(); ++ i) {

		if (!silent) std::cout << "\t\tprocess file: " << files[i] << "\n";

		std::ifstream fh;
		xny::openfile<std::ifstream>(fh, files[i]);

		bio::fastq_input_iterator<> fq(fh), end;

		while (fq != end) {

			add_fq_reads_only (seqs, cnt, fq, end);

			if ((int) seqs.size() >= batch) {
				total_reads += seqs.size();
				kmer_cnt_in_seqs (kcnt, seqs, k);
				cnt = batch;
				seqs.clear();
			} else { // not enough reads to fill in batch for current file
				cnt = batch - seqs.size();
				xny::closefile(fh);
			}
		} // while
	} // for
	total_reads += seqs.size();
	kmer_cnt_in_seqs (kcnt, seqs, k);

	if (!silent) std::cout << "\t\t" << total_reads << " reads analyzed\n";
} // obtain_kfrq

/* @brief	Count kmers in a given set of sequences [seqs]
 */
void kmer_cnt_in_seqs (umap_t& kcnt, const strvec_t& seqs, int k){
	int sz = seqs.size();
	#pragma omp parallel
	{
		umap_t local_kcnt; // private to each thread

		#pragma omp for
		for (int i = 0; i < sz; ++ i) {
			uvec_t kIDs;
			xny::get_bitkmer<std::back_insert_iterator<uvec_t>, uint32_t>
				(seqs[i], std::back_inserter(kIDs), k, 0);
			uset_t kID_set (kIDs.begin(), kIDs.end());

			for (auto& id: kID_set) {
				umap_t::iterator it = local_kcnt.find(id);
				if (it != local_kcnt.end()) ++ it->second;
				else local_kcnt[id] = 1;
			}
		}
		#pragma omp critical
		{
			for (auto& kc: local_kcnt) {
				umap_t::iterator it = kcnt.find(kc.first);
				if (it != kcnt.end()) it->second += kc.second;
				else kcnt.insert(kc);
			}
		}
	} // #pragma omp parallel

} // kmer_cnt_in_seqs

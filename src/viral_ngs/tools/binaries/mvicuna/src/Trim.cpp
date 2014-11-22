//========================================================================
// Project     : M-Vicuna
// Name        : Trim.cpp
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

#include "Trim.h"
/**	Function trimming ()
 *
 */
void trimming (const strvec_t& ipfq, const strvec_t& isfq,
	const trm_t& trm, xny::low_complexity& lc, int batch, bool silent){

	// sanity check
	if ((trm.op.size() == ipfq.size() && trm.os.size() == ipfq.size()/2) ||
		(trm.op.size() == 2 && trm.os.size() == 1)) {	}
	else abording ("In trimming() SC failed");

	// ---------------  read input vector --------------------------
	std::cout << "\tRead in vectors ...";
	strvec_t vectors;
	kindex_t kindex; //[kmerID -> list (vecID, vecPos, dir)]
	process_vector_file (vectors, kindex, trm.vecfa, std::min (trm.min_match, 16));

	if (vectors.size() == 0 || kindex.empty()) {
		std::cout << "\t\tno vector trimming applied\n\n";
		//return;
	}
    // --------------- prepare output when applicable -----------------
	std::ofstream ofhfq, ofhfq2, ofhs;
	if (trm.op.size() == 2 && trm.os.size() == 1) {
		xny::openfile<std::ofstream> (ofhfq, trm.op[0]);
		xny::openfile<std::ofstream> (ofhfq2, trm.op[1]);
		xny::openfile<std::ofstream> (ofhs, trm.os[0]);
	}

	// process every pair of fastq files
	int num_file_pairs = ipfq.size()/2;

	for (int i = 0; i < num_file_pairs; ++ i) {

		int fID = 2*i;

		if (! silent) {
			std::cout << "\tprocess files: " << ipfq[fID] << " and "
					<< ipfq[fID + 1] << "\n\n";
		}

		// ----- output non-redundant read-pairs ---------
		if (trm.op.size() > 2) {
			xny::openfile<std::ofstream> (ofhfq, trm.op[fID]);
			xny::openfile<std::ofstream> (ofhfq2, trm.op[fID + 1]);
			xny::openfile<std::ofstream> (ofhs, trm.os[i]);
		}

		trim_pfq (ipfq[fID], ipfq[fID+1], vectors, kindex, ofhfq,
				ofhfq2, ofhs, trm, lc, batch);

		if (trm.op.size() > 2) {
			xny::closefile(ofhfq);
			xny::closefile(ofhfq2);
			xny::closefile (ofhs);
		}
	} // for (int i = 0

	// process singleton fastq files, the result is stored in the
	// last specified output fastq file
	int num_sfiles = isfq.size();
	std::ofstream ofhsfq;
	xny::openfile<std::ofstream> (ofhsfq, trm.os.back(), std::ios::app);
	for (int fID = 0; fID < num_sfiles; ++ fID) {
		if (! silent) std::cout << "\tprocess files: " << isfq[fID] << "\n\n";
		trim_sfq (ofhsfq, isfq[fID], trm, lc, batch);
	}
	xny::closefile(ofhsfq);
} // trimming


/**	Function process_vector_file ()
 *
 * Given an input fasta file [vecfa], record all vectors in [vectors] and
 * process each vector kmer then store information in
 * [kindex]: [kmerID -> list (vecID, vecPos, dir)]
 */
void process_vector_file (strvec_t& vectors, kindex_t& kindex,
		const std::string& vecfa, int k) {

	/* ------------ read input vectors ----------- */
	std::ifstream fh;
	xny::openfile<std::ifstream>(fh, vecfa);
	bio::fasta_input_iterator<> iter_fa(fh), end;
	add_fa_reads_only (vectors, INT_MAX, iter_fa, end);
	std::cout << "\t\t" << vectors.size() << " vector sequences read\n";
	if (vectors.size() == 0) return;
	xny::closefile (fh);

	/* generate kindex information [kmerID -> list (vecID, vecPos, dir)] */
	int vectorID = 0;
	for (auto& v: vectors) {

		// make sure vector doesn't contain 'n' or 'N'
		if (std::isupper(v.at(0))) std::replace (v.begin(), v.end(), 'N', 'A');
		else std::replace (v.begin(), v.end(), 'n', 'a');

		// generate (k)-spectrum for fwd & rvc strands of vector[i]
		uvec_t kIDs;
		xny::get_bitkmer<std::back_insert_iterator<uvec_t>, uint32_t>
			(v, std::back_inserter(kIDs), k, 2);

		for (unsigned int j = 0; j < kIDs.size(); ++ j) {
			kindex_t::iterator it  = kindex.find(kIDs[j]);
			if (it != kindex.end()) {
				if (j % 2 == 0) { // forward strand
					it->second.push_back(kloc_t(vectorID, j/2, true));
				} else it->second.push_back(kloc_t(vectorID, j/2, false));
			} else {
				if (j % 2 == 0) kindex[kIDs[j]] = { kloc_t(vectorID, j/2, true) };
				else kindex[kIDs[j]] = { kloc_t (vectorID, j/2, false) };
			}
		} // for (unsigned int j

		vectorID ++;
	} // for (auto & v: vectors )

}// process_vector_file


/* @brief Function trim_sfq ()
 *
 * Apply trimming to single end fastq files
 * -- low quality score criterion
 * -- low complexity criterion
 */
void trim_sfq (std::ofstream& ofhsfq, const std::string& ifq,
		const trm_t& trm, xny::low_complexity& lc, int batch) {

	std::ifstream ifhfq;
	xny::openfile<std::ifstream>(ifhfq, ifq);
	bio::fastq_input_iterator<> iter_fq (ifhfq), iter_end;

	int total_reads = 0;
	int num_trimmed = 0;
	std::vector<fqtuple_t> reads;

 	while (iter_fq != iter_end) {
 		add_fq_reads (reads, batch, iter_fq, iter_end);

 		int read_cnt = reads.size();
		#pragma omp parallel for
 		for (int i = 0; i < read_cnt; ++ i) {
 			if (trim_lc_lq (std::get<1>(reads[i]), std::get<2>(reads[i]),
 					trm.min_qual, trm.min_rlen, lc)) {
 				++ num_trimmed;
 			}
 		}
 		// write out
 		for (int i = 0; i < read_cnt; ++ i) {
 			if(! (std::get<1>(reads[i])).empty()) {
 				ofhsfq << "@" << std::get<0>(reads[i]) << "\n";
 				ofhsfq << std::get<1>(reads[i]) << "\n";
 				ofhsfq << "+\n";
 				ofhsfq << std::get<2>(reads[i]) << "\n";
 			}
 		}
		total_reads += read_cnt;

		reads.clear();

	} // while

 	xny::closefile(ifhfq);

	std::cout << "\t\ttotal reads: " << total_reads << ", " << num_trimmed << " trimmed\n";
} // trim_sfq

/**	Function trim_paired_fqfiles ()
 *
 * Apply trimming to two paired input fastq files
 * -- primer trimming
 * -- low quality score criterion
 * -- low complexity criterion
 */
void trim_pfq (const std::string& ifq, const std::string& ifq2,
	const strvec_t& vectors, const kindex_t& kindex, std::ofstream& ofhfq,
	std::ofstream& ofhfq2, std::ofstream& ofhs, const trm_t& trm,
	 xny::low_complexity& lc, int batch){

	std::ifstream ifhfq, ifhfq2;
	xny::openfile<std::ifstream>(ifhfq, ifq);
	xny::openfile<std::ifstream>(ifhfq2, ifq2);
	bio::fastq_input_iterator<> iter_fq (ifhfq), iter_end, iter_fq2 (ifhfq2);

	int total_read_pairs = 0;
	int num_trimmed = 0;
	std::vector<fqtuple_t> pairs;

 	while (iter_fq != iter_end && iter_fq2 != iter_end) {

 		add_fq_reads (pairs, batch/2, iter_fq, iter_end);
		add_fq_reads (pairs, batch/2, iter_fq2, iter_end);

		num_trimmed += apply_trimming (pairs, vectors, kindex, lc, trm);

		// output to file
		int fragnum = pairs.size()/2;
		for (int i = 0; i < fragnum; ++ i) {
			if ((std::get<1>(pairs[i])).empty()) { // first pair is empty
				if (! std::get<1>(pairs[i + fragnum]).empty()) { // 2nd not empty
					ofhs << "@" << std::get<0>(pairs[i + fragnum]) << "\n";
					ofhs << std::get<1>(pairs[i + fragnum]) << "\n";
					ofhs << "+\n";
					ofhs << std::get<2>(pairs[i + fragnum]) << "\n";
				}
			} else if (std::get<1>(pairs[i + fragnum]).empty()) { // 2nd empty
				ofhs << "@" << std::get<0>(pairs[i]) << "\n";
				ofhs << std::get<1>(pairs[i]) << "\n";
				ofhs << "+\n";
				ofhs << std::get<2>(pairs[i]) << "\n";
			} else {
				ofhfq << "@" << std::get<0>(pairs[i]) << "\n";
				ofhfq << std::get<1>(pairs[i]) << "\n";
				ofhfq << "+\n";
				ofhfq << std::get<2>(pairs[i]) << "\n";

				ofhfq2 << "@" << std::get<0>(pairs[i + fragnum]) << "\n";
				ofhfq2 << std::get<1>(pairs[i + fragnum]) << "\n";
				ofhfq2 << "+\n";
				ofhfq2 << std::get<2>(pairs[i + fragnum]) << "\n";
			}
		}

		total_read_pairs += pairs.size()/2;

		pairs.clear();

	} // while

 	xny::closefile(ifhfq);
	xny::closefile(ifhfq2);

	std::cout << "\t\ttotal reads: " << total_read_pairs * 2 << ", "
	 	 << num_trimmed << " trimmed\n";
} //trim_pfq

/** Function apply_trimming ()
 *
 * Apply trimming to each read in [seq] return number of trimmed reads,
 * no paired information is used
 */
int apply_trimming (std::vector<fqtuple_t>& seq, const strvec_t& vectors,
		const kindex_t& kindex, xny::low_complexity& lc, const trm_t& trm) {

	int sz = seq.size();
	int k = std::min (trm.min_match, 16);

	bvec_t trimmed (sz, false); // record which reads are trimmed
	#pragma omp parallel for
	for (int i = 0; i < sz; ++ i) {
		std::string rseq = std::get<1>(seq[i]);
		std::string qual = std::get<2>(seq[i]);
		// make sure vector doesn't contain 'n' or 'N'
		if (std::isupper(rseq.at(0))) {
			std::replace (rseq.begin(), rseq.end(), 'N', 'A');
		} else std::replace (rseq.begin(), rseq.end(), 'n', 'a');

		// Primer Trimming: checking each read against all primers
		while (true) {
			bool trim_applied = false;

			uvec_t rkmers;
			xny::get_bitkmer<std::back_insert_iterator<uvec_t>, uint32_t>
				(rseq, std::back_inserter(rkmers), k, 3);

			for (unsigned int rPos = 0; rPos < rkmers.size(); ++ rPos) {

				kindex_t::const_iterator it  = kindex.find(rkmers[rPos]);

				// check if we can do trimming using current kmer: kmers[j]
				if (it != kindex.end()) {
					if (trim_primer (rseq, qual, it->second,
							vectors, rPos, trm)) {
						trimmed[i] = true;
						trim_applied = true;
						break;
					}
				}
			} // for rPos

			// break when no more trimming can be applied to the read
			if (!trim_applied) break;
		}
		// Low quality score and low complexity trimming
		if (trim_lc_lq (rseq, qual, trm.min_qual, trm.min_rlen, lc)) {
			trimmed[i] = true;
		}

		std::get<1>(seq[i]) = rseq;
		std::get<2>(seq[i]) = qual;
	} // for i

	// counting the trimmed reads
	int num_trimmed = 0;
	for (int i = 0; i < sz; ++ i) { if (trimmed[i]) ++ num_trimmed; }
	return num_trimmed;

} // apply_trimming

/** Function trim_primer ()
 *
 * Apply primer trimming to a read sequence [rSeq] and quality string [qual]
 * For every kmer shared between [rSeq] and an input vector sequence
 * return once the first trimming event happens when
 * the minimum match length is identified; trim is either applied to prefix,
 * suffix, or the complete read by choosing whichever the remaining
 * fragment is longer post-trimming.
 */
bool trim_primer (std::string& rSeq, std::string& qual,
		const std::vector<kloc_t>& kloclist,
		const strvec_t& vectors,	int rPos, const trm_t& trm) {

	int k = std::min (16, trm.min_match);

	int rBegin = rPos, rEnd = rPos + k - 1;
 	int rLen = rSeq.length();

 	for (auto& kloc: kloclist) { // for each matched kmer locus

		std::string myVec = vectors[std::get<0> (kloc)];// vector sequence
		int vBegin = std::get<1> (kloc), vEnd = vBegin + k - 1;
		bool vDir = std::get<2> (kloc);	// fwd or rv direction

		std::string rPrf = rSeq.substr(0, rBegin),
				    rSuf = rSeq.substr(rEnd+1, rLen - rEnd),
					vPrf = myVec.substr(0, vBegin),
					vSuf = myVec.substr(vEnd+1, myVec.length() - vEnd);

		if (vDir) { // fwd
			rEnd += xny::common_pref_len (rSuf, vSuf);
			rBegin -= xny::common_suf_len (rPrf, vPrf);
		} else { // rv
			rEnd += xny::common_pref_len (rSuf, xny::get_rvc_str (vPrf));
			rBegin -= xny::common_suf_len (rPrf, xny::get_rvc_str(vSuf));
		}

		if (rEnd - rBegin + 1 >= trm.min_match) { // apply trimming
			int remain_suffix_len = rLen - rEnd - 1,
				remain_prefix_len = rBegin + 1;
			if (remain_suffix_len >= remain_prefix_len &&
				remain_suffix_len >= trm.min_rlen) { //trim prefix
				rSeq = rSeq.substr(rEnd + 1, remain_suffix_len);
				qual = qual.substr(rEnd + 1, remain_suffix_len);
			} else if (remain_prefix_len >= remain_suffix_len &&
					remain_prefix_len >= trm.min_rlen) { // trim suffix
				rSeq = rSeq.substr(0, remain_prefix_len);
				qual = qual.substr(0, remain_prefix_len);
			} else {
				rSeq = ""; // rm whole read
				qual = "";
			}

			return true;
		}
 	} // for

	return false;
} // trim_primer

/* @brief	Trim read and quality string by min quality score and
 *	low complexity criteria
 */
bool trim_lc_lq (std::string& rSeq, std::string& qual, int min_qual,
		int min_len, xny::low_complexity& lc) {
	/*	Trim max prefix and suffix where quality score of each base <= min_qual
	 */
	int rlen = qual.length();
	int first_hq_idx = 0, last_hq_idx = rlen - 1;
	for (auto& c : qual) { // prefix
		if ((int) c > min_qual) break;
		++ first_hq_idx;
	}
	bool last_pos_found = false;
	for (; last_hq_idx >= min_len + first_hq_idx - 1; -- last_hq_idx) {
		if ((int) qual.at(last_hq_idx) - 33 > min_qual) {
			last_pos_found = true;
			break;
		}
	}
	if (last_pos_found) {
		rSeq = rSeq.substr(first_hq_idx, last_hq_idx - first_hq_idx + 1);
		if (! lc (rSeq)) {
			qual = qual.substr(first_hq_idx, last_hq_idx - first_hq_idx + 1);
			if (rSeq.length() != rlen) return true;
			else return false;
		}
	}
	rSeq.clear();
	qual.clear();
	return true;
} // trim_lc_lq


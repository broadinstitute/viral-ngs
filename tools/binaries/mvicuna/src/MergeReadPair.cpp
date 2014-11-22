//========================================================================
// Project     : M-Vicuna
// Name        : MergeReadPair.cpp
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

#include "MergeReadPair.h"
#include "ReadBioFile.h"

void debug_output_to_file (std::ofstream& ofh, std::vector<fqtuple_t>& seq){
	for (int i = 0; i < seq.size(); ++ i) {
		ofh << "@" << std::get<0> (seq[i]) << "\n"
				<< std::get<1> (seq[i]) << "\n+\n"
				<< std::get<2> (seq[i]) << "\n";
	}
}

void merge_paired_read (const strvec_t& ifq, const strpair_t& ofq,
		const std::string& os, int batch) {

	// sanity check
	if (os.empty()) abording ("merge_paired_read: ofa is empty");

	std::ofstream ofhfq, ofhfq2, ofhs;
	xny::openfile<std::ofstream> (ofhfq, ofq.first);
	xny::openfile<std::ofstream> (ofhfq2, ofq.second);
	xny::openfile<std::ofstream> (ofhs, os);

	int num_merged_pairs = 0, total_read_pairs = 0;
	std::vector<fqtuple_t> seq;
	int cnt = batch;
	for (unsigned int i = 0; i < ifq.size(); i += 2) {

		std::cout << "\tprocess files: " <<  ifq[i] << "\t" << ifq[i + 1] << "\n\n";

		std::ifstream fh, fh2;
		xny::openfile<std::ifstream>(fh, ifq[i]);
		xny::openfile<std::ifstream>(fh2, ifq[i+1]);

		/*{
		std::ofstream debug_output;
		xny::openfile<std::ofstream>(debug_output, "output/test.fq");
		}*/
		bio::fastq_input_iterator<> fq(fh), end;
		bio::fastq_input_iterator<> fq2(fh2);

		//int debug = 1;
		while ((fq != end) && (fq2 != end)) {

			add_fq_reads (seq, cnt/2, fq, end);

			/*{ // check if input is properly read
			//debug_output_to_file (debug_output, seq);
			}*/

			add_fq_reads (seq, cnt/2, fq2, end);

			/*{
				std::cout << debug << "\t" << seq.size() << "\t" << cnt << "\n";
				++ debug;
			}*/

			if ((int) seq.size() >= batch) {
				total_read_pairs += seq.size()/2;
				num_merged_pairs += apply_merging (ofhs, ofhfq, ofhfq2, seq);
				cnt = batch;
				seq.clear();
			} else { // not enough reads to fill in batch for current file pair
				cnt = batch - seq.size();
				xny::closefile(fh);
				xny::closefile(fh2);
			}
		} // while

	}
	total_read_pairs += seq.size()/2;
	num_merged_pairs += apply_merging (ofhs, ofhfq, ofhfq2, seq);

	std::cout << "\tnumber of merged pairs vs total: " << num_merged_pairs
			<< " vs " << total_read_pairs << " ("
			<< 100.0*num_merged_pairs /total_read_pairs << "%) \n";

	xny::closefile (ofhfq);
	xny::closefile (ofhfq2);
	xny::closefile (ofhs);

} // merge_paired_read



/* @brief	Given n read-pairs (1, ...n, n+1, ..., 2n) stored in [seq],
 * where the reads (i, n+i) for 1 <= i <= n form a pair. Every pair is
 * attempted for merging, the ones failed merging are added to [fq] [fq2]
 * and the merged ones are added to [fhs] in single end fastq format
 *
 * NOTE: assume the paired reads are from different strands either:
 * 	----->
 * 		<-----
 * 		OR abnormal case
 * 	<-----
 * 	    ----->
 */
int apply_merging (std::ofstream& fhs, std::ofstream& fhfq,
		std::ofstream& fhfq2, std::vector<fqtuple_t>& seq) {

	bool debug = false;
	int debug_cnter = 0;
	// reverse complementary of second pair
	int num = seq.size();

	#pragma omp parallel for
	for (int i = num/2; i < num; ++ i) {
		xny::rvc_str(std::get<1> (seq[i]));
	}

	xny::suffix_prefix_gap_free_aln aligner (7, 90, 1);

	std::vector<strpair_t> list_merged_fragments (num/2);

	if (debug){ // debug case for merging
		omp_set_num_threads(1);
	}

	#pragma omp parallel for
	for (int i = num/2; i < num; ++ i) {

		coord_t coord =	aligner (std::get<1> (seq[i-num/2]),
				std::get<1> (seq[i]));

		if (std::get<0> (coord) != -1) {
			int start0 = std::get<0> (coord), end0 = std::get<1> (coord),
				start1 = std::get<2> (coord), end1 = std::get<3> (coord),
				l0 = (std::get<1> (seq[i-num/2])).length(),
				l1 = (std::get<1> (seq[i])).length();

			/*if (debug) { // debug print
				std::cout << start0 << "\t" << end0 << "\t"
						<< start1 << "\t" << end1 << "\n";
				int num_dash = start0 - start1;
				std::cout << std::get<1> (seq[i-num/2]) << "\n";
				for (int g = 0; g < num_dash; ++ g) std::cout << "-";
				std::cout << std::get<1> (seq[i]) << "\n\n";
			}*/

			// to merge
			int overlap_sz = end0 - start0 + 1;
			std::string merged_seq, merged_qual;
			for (int idx = 0; idx < overlap_sz; ++ idx) {
				char c0 = (std::get<1> (seq[i-num/2])).at(start0 + idx),
					 c1 = (std::get<1> (seq[i])).at(start1 + idx),
					 q0 = (std::get<2> (seq[i-num/2])).at(start0 + idx),
					 q1 = (std::get<2> (seq[i])).at(start1 + idx);

				if (c0 == c1) {
					merged_seq += c0;
					if (q0 >= q1) merged_qual += q0;
					else merged_qual += q1;
				} else {
					if (q0 >= q1) {
						merged_seq += c0;
						merged_qual += q0;
					} else {
						merged_seq += c1;
						merged_qual += q1;
					}
				}
			} // for (int idx = 0;

			//if (debug) std::cout << "Merged part: \n" << merged_seq << "\n";

			// generate the full string
			if (start0 + 1 >= l0 - end0 - 1) { // normal direction
				merged_seq = (std::get<1> (seq[i-num/2])).substr(0, start0) + merged_seq;
				merged_qual = (std::get<2> (seq[i-num/2])).substr(0, start0) + merged_qual;
				if (end1 < l1 - 1) {
					merged_seq += (std::get<1> (seq[i])).substr(end1 + 1, l1 - end1);
					merged_qual += (std::get<2> (seq[i])).substr(end1 + 1, l1 - end1);
				}
			} else { // abnormal direction
				merged_seq = (std::get<1> (seq[i])).substr(0, start1) + merged_seq;
				merged_qual = (std::get<2> (seq[i])).substr(0, start1) + merged_qual;
				if (end0 < l0 - 1) {
					merged_seq += (std::get<1> (seq[i-num/2])).substr(end0 + 1, l0 - end0);
					merged_qual += (std::get<2> (seq[i-num/2])).substr(end0 + 1, l0 - end0);
				}

				if (debug) {
					std::cout << "> " << debug_cnter << "\n";
					std::cout << merged_seq << "\n\n";

					std::cout << std::get<0> (seq[i-num/2]) << "\n" << std::get<1> (seq[i-num/2]) << "\n";
					std::cout << std::get<0> (seq[i]) << "\n" << std::get<1> (seq[i]) << "\n\n";

					++ debug_cnter;
				}
			}

			/*
			if (debug){
				std::cout << "post-merging, we have\n" << merged_seq
						<< "\n" << merged_qual << "\n";
			}*/

			list_merged_fragments[i - num/2] = strpair_t (merged_seq,
					merged_qual);


		}  // if (std::get<0> (coord) != -1)

	} // 	for (int i = num/2; i < num; ++ i) {

	/* write to file */
	int merged_cnt = 0;
	for (int i = 0; i < num/2; ++ i) {
		if (!list_merged_fragments[i].first.empty()) {
			fhs << "@" << std::get<0> (seq[i + num/2]) << "\n"
					<< list_merged_fragments[i].first << "\n+\n"
					<< list_merged_fragments[i].second << "\n";
			++ merged_cnt;
		} else {
			fhfq << "@" << std::get<0> (seq[i]) << "\n"
					<< std::get<1> (seq[i]) << "\n+\n"
					<< std::get<2> (seq[i]) << "\n";
			fhfq2 << "@" << std::get<0> (seq[i + num/2]) << "\n"
					<< std::get<1> (seq[i + num/2])
					<< "\n+\n" << std::get<2> (seq[i + num/2]) << "\n";
		}
	}

	return merged_cnt;
	//std::cout << "Number of merged read pairs = " << merged << "\n";
	//
	// suffix_prefix_gap_free_aln ();

} // apply_merging


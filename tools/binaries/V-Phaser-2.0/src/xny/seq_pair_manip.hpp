//============================================================================
// Project     : Diversifier
// Name        : seq_pair_manip.hpp
// Author      : Xiao Yang
// Created on  : Sep 26, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#ifndef SEQ_PAIR_MANIP_HPP_
#define SEQ_PAIR_MANIP_HPP_

/*
 * This header file is used along with NCBI c++ tool kit, when compiling
 * needs to specify the path -I/NCBI_tool_kit_path/
 */

/* NCBI library */
#include <algo/align/nw/nw_aligner.hpp>

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <map>
#include "seq_manip.hpp"
//#include "Util.h" // for debugging purprose


namespace xny{

	typedef std::vector<unsigned long> lvec_t;
	typedef std::vector<uint32_t> uvec_t;
	typedef std::vector<int> ivec_t;
	typedef std::pair<uint32_t, uint32_t> upair_t;
	typedef std::pair<int, int> ipair_t;

	// for comparing upair_t by first element
	struct cmp_upair {
		bool operator () (const upair_t& lhs, const upair_t& rhs) const {
			return lhs.first < rhs.first;
		}
	};
	struct cmp_pair_ge {
		bool operator () (const std::pair<int, lvec_t>& lhs,
				const std::pair<int, lvec_t>& rhs) const {
			return rhs.first < lhs.first;
		}
	};
	/* the structure to store alignment coordinates between two strings
	 * s0[s0_l, s0_r] is aligned with s1[s1_l, s1_r]
	 */
	typedef struct ALIGN {
		ALIGN(int f0, int t0, int f1, int t1, const std::string& s):
			s0_l (f0), s0_r (t0), s1_l (f1), s1_r (t1), seq (s) {};
		ALIGN() {};
		int s0_l;
		int s0_r;
		int s1_l;
		int s1_r;
		std::string seq;
	}align_t;


	inline int hamming_dist (const std::string& lhs, const std::string& rhs){
		int len = lhs.length(), hd = 0;
		for (int i = 0; i < len; ++ i)
			if (std::toupper(lhs.at(i)) != std::toupper(rhs.at(i))) ++ hd;
		return hd;
	}


	inline void nwalign(std::string& result, ncbi::CNWAligner& aligner,
		const std::string& s0, const std::string& s1, const lvec_t& hits) {

		/* result consists of (M)atch, (R)eplace, (I)nsert, and (D)elete
		 * frame shift is not penalized, e.g., n*three consecutive gaps */

		aligner.SetSequences(s0.c_str(), s0.length(), s1.c_str(), s1.length());
		aligner.SetPattern(hits);
		aligner.SetEndSpaceFree(true, true, true, true);
	    aligner.SetWg(-3); // set gap opening score
	    aligner.SetWs(-1); // gap ext score
	    aligner.SetWm(1); // set match score
	    aligner.SetWms(-1);
	    aligner.SetScoreMatrix(0);
		//int score = aligner.Run();
	    aligner.Run();
		result = aligner.GetTranscriptString();
	}


	/* @brief	Given an alignment sequence [aln], which consists of
	 * 			alphabet {M, R, D, I}, return a maximum prefix [par_aln]
	 *			with percent identity >= [min_perc_iden].
	 * @note 	Any block of gaps (Is or Ds) are not penalized if
	 * 			1) the block length is 3x
	 * 			2) the block length is <= [max_gap_len]
	 * 			3) both sides of this block is supported by the same
	 * 			number of Ms as the block length
	 */
	inline void get_max_constrained_aln_prefix (std::string& par_aln,
			std::string& aln, int min_perc_iden, int max_gap_len) {
		int max_pos = -1;
		int M_cnt = 0, D_cnt = 0, I_cnt = 0;
		int M_b4_gap = 0, M_after_gap = 0;
		int valid_gaps = 0;
		for (int i = 0; i < (int) aln.length(); ++ i) {
			switch (aln.at(i)){
			case 'M':
				++ M_cnt;
				if (I_cnt || D_cnt) {
					++ M_after_gap;
					if (I_cnt % 3 != 0 || D_cnt % 3 != 0
						|| (I_cnt + D_cnt > max_gap_len)) {
							I_cnt = 0;
							D_cnt = 0;
							M_b4_gap = M_after_gap;
					} else {
						if (M_b4_gap >= I_cnt + D_cnt) {
							if (M_after_gap >= I_cnt + D_cnt) {
								valid_gaps += I_cnt + D_cnt;
								I_cnt = 0;
								D_cnt = 0;
								M_b4_gap = M_after_gap;
							}
						} else {
							I_cnt = 0;
							D_cnt = 0;
							M_b4_gap = M_after_gap;
						}
					}
				} else ++ M_b4_gap;

				/* test if the current position meets min_perc_iden */
				if (100.0 * M_cnt/ (i + 1 - valid_gaps) >= min_perc_iden) {
					max_pos = std::max(i, max_pos);
				}
				break;
			case 'D':
				++ D_cnt;
				I_cnt = 0;
				break;
			case 'I':
				++ I_cnt;
				D_cnt = 0;
				break;
			default:
				break;
			}
		}
		if (max_pos >= 0) par_aln = aln.substr(0, max_pos + 1);
		else par_aln = "";
	}

	/* @brief	Given an alignment result [alignment], which consists of
	 * 			alphabet {M, R, D, I}, return percent identity.
	 * 			3x gaps will not be penalized as long as on both sides of
	 * 			this gap, there are at least [min_gap_support] number of
	 * 			matches e.g., -MMMRIIIMMM-, -MMMDRDDMMM-
	 * @return (#M)/(length of [alignment])
	 */
	inline int align_identity (const std::string& alignment,
			int min_gap_support) {

		int M_cnt = 0, I_counter = 0, D_counter = 0;
		int M_counter_b4_gap = 0, M_counter_after_gap = 0;
		int adj_align_len = alignment.length();

		int total_indels = 0;
		/* R_flag to avoid 3x gap in the form -IIRIM- but allow -IIIRM- */
		//bool R_flag = false;

		for (int j = 0; j < adj_align_len; ++ j) {

			switch (alignment.at(j)){
			case 'I':
				/* if (R_flag) { I_counter = 1;	R_flag = false; } else ++ I_counter;	*/
				++ total_indels;
				++ I_counter;
				D_counter = 0;
				break;
			case 'D':
				/* if (R_flag) { D_counter = 1; R_flag = false; } else ++ D_counter;	*/
				++ total_indels;
				++ D_counter;
				I_counter = 0;
				break;
			case 'R':
				//R_flag = true;
				break;
			case 'M':
				++ M_cnt;
				if (I_counter || D_counter) {
					++ M_counter_after_gap;
					if (I_counter % 3 != 0 || D_counter % 3 != 0
							|| (I_counter + D_counter > 15)) {
						I_counter = 0;
						D_counter = 0;
						M_counter_b4_gap = M_counter_after_gap;
					} else {
						if (M_counter_b4_gap >= I_counter + D_counter) {
							if (M_counter_after_gap >= I_counter + D_counter) {
								adj_align_len -= I_counter + D_counter;
								I_counter = 0;
								D_counter = 0;
								M_counter_b4_gap = M_counter_after_gap;
							}
						} else {
							I_counter = 0;
							D_counter = 0;
							M_counter_b4_gap = M_counter_after_gap;
						}
					}
				} else ++ M_counter_b4_gap;

				//R_flag = false;
				break;
			default:
				std::cout << "[WARNING] align_identity, found "
						"unknown character:" << alignment.at(j) << "\n";
				break;
			}
		}

		/*
		{// debug print
			if (adj_align_len != alignment.length() && total_indels >= 9) {
				std::cout << "M_cnt = " << M_cnt << "\tadj" << adj_align_len
						<< "\t #indels=" << total_indels << "\n";
			}
		}
		*/

		return (100*M_cnt/adj_align_len);

	} // align_identity

	/* @brief 	Check if the current k-frame satisfies the following:
	 * 			a) (k-frame length / contig overlap length) >= [min_perc]
	 * 			b) if the k-frame is within [bd_dist] away from ONE end of
	 * 				either contig under comparison (count for dove-tail
	 * 				alignment)
	 */
	inline bool sat_dist_cri (int l0, int l1, const lvec_t& hits,
			int min_perc, int bd_dist) {

		int hs = hits.size();
		int k_frame_len = std::max (	hits[hs - 1] - hits[2]  + 1,
								hits[hs - 3] - hits[0] + 1);
		int l_min_hang = std::min(hits[0], hits[2]);
		int r_min_hang =	 std::min(l0 - hits[hs - 3] - 1,
								  l1 - hits[hs - 1] - 1);
		int overlap = l_min_hang + r_min_hang + k_frame_len;
		if ((100 * k_frame_len / overlap >= min_perc)
				&&
			(hits[0] < bd_dist || (l0 - hits[hs - 3]) <= bd_dist + 1||
			 hits[2] < bd_dist || (l1 - hits[hs - 3]) <= bd_dist + 1)) {

			return true;
		}
		return false;
	}

	/*
	 * @brief	This aligner is meant for identifying prefix-suffix,
	 * 			containment, dove-tail alignment b/t two input DNA
	 * 			sequences [s0] and [s1], using kmers as seeds.
	 *
	 * 			Only the longest prefix-suffix alignment (if any) will be
	 * 			reported, which satisfies the following criteria:
	 * 			1) the similarity b/t aligned region is >= [min_iden]
	 * 			2) the length of the aligned region is >= [min_overlap]
	 * 			3) exists an overhang (out of four) is within [max_overhang]
	 * 			distance to the end of the sequence.
	 *
	 * @method	First identify all kmers that occur in both s0, s1
	 * 			(only fwd strds are considered). For each common	kmer,
	 * 			generate a pair <p0, p1> specifying its positions
	 * 			in s0 and s1 respectively, this forms a position vector
	 * 			P. Next, extend kmers into frames using func "sat_dist_cri",
	 * 			based on kmer distance and distance to end of sequence.
	 *			Then, generate alignment for each kmer-frame using
	 *			frame anchoring and nw alignment, which further validated
	 *			by func "align_identity". Once the validation is true,
	 *			try extend this alignment towards left and right sides,
	 *			when applicable.
	 *
	 * @para
	 * 	[out]: align_t type, output
	 * 	[s0] & [s1]: input strings to be aligned
	 *	[K]: 	kmer length for seeding the alignment
	 *
	 * @note		When counting identity in the alignment result, gaps with
	 * 			the multiple of 3 are not penalized.
	 * @assure	s0 and s1 are both either lower case or upper case
	 *
	 */
	inline void kmer_anchor_ps_aligner (align_t& out, const std::string& s0,
		const std::string& s1, int K, int min_overlap, int min_iden,
		int max_overhang) {

		bool debug = false;

		/* generate kmers in uint32 bits form for both s0 and s1 */
		uvec_t s0_kmers, s1_kmers;
		xny::get_bitkmer<std::back_insert_iterator<uvec_t>, uint32_t> (
				s0, std::back_inserter(s0_kmers), K, 3);
		xny::get_bitkmer<std::back_insert_iterator<uvec_t>, uint32_t> (
				s1, std::back_inserter(s1_kmers), K, 3);

		/* generate (kmer, pos) pair in s1 */
		std::set<upair_t, cmp_upair> s1_KP;
		for (int i = 0; i < (int) s1_kmers.size(); ++ i)
			s1_KP.insert(upair_t(s1_kmers[i], i));

		/* identify common kmer positions in s0 and s1 */
		std::vector<ipair_t> kp_s0_s1;
		for (int i = 0; i < (int) s0_kmers.size(); ++ i) {
			std::set<upair_t, cmp_upair>::iterator
				it_s1kp = s1_KP.find(upair_t(s0_kmers[i], -1));
			if (it_s1kp != s1_KP.end()) { // found !
				kp_s0_s1.push_back(ipair_t(i, it_s1kp->second));
			}
		}

		if (debug) { // debug
			std::cout << "ALIGN:\n" << ">s0\n" << s0 << "\n"
									<< ">s1\n" << s1 << "\n";
			if (kp_s0_s1.size() == 0) std::cout << "\t no common kmers\n";
		}

		/* no common kmers b/t s0 & s1 */
		if (!kp_s0_s1.size()) return;

		/* merge kmers to form kmer-frames
		/* [block_dist]: an empirical value to determine whether two
		 * 				 blocks are close enough to be merged, and if
		 * 				 the kframe is close enough to contig boundaries
		 * [k_frame_min_perc]: an empirical min_value for the ratio b/t
		 * 			     the kframe length and the length of contig overlap
		 *  				 anchored by this kframe
		 */
		int block_dist = std::max (80, 5 * K);
		int k_frame_min_perc = 50;
		int l0 = s0.length(), l1 = s1.length();
		std::vector<std::pair<int, lvec_t> > frames;
		lvec_t hits;

		std::vector<bool> considered (kp_s0_s1.size(), false);
		for (int i = 0; i < (int) kp_s0_s1.size() - 1; ++ i) {
			if (considered[i]) continue;
			hits.clear();
			hits.push_back (kp_s0_s1[i].first);
			hits.push_back (kp_s0_s1[i].first + K - 1);
			hits.push_back (kp_s0_s1[i].second);
			hits.push_back (kp_s0_s1[i].second + K - 1);

			for (int j = i + 1; j < (int) kp_s0_s1.size(); ++ j) {
				if (considered[j]) continue;

				int s0_step = kp_s0_s1[j].first - hits[hits.size() - 3];
				int s1_step = kp_s0_s1[j].second - hits[hits.size() - 1];
				if (s0_step <= 0) {	/* s0 overlap */
					if (s0_step == s1_step) { /* extend */
						hits[hits.size() - 3] = kp_s0_s1[j].first + K-1;
						hits[hits.size() - 1] = kp_s0_s1[j].second + K-1;
						considered[j] = true;
					} else {
						if (s1_step <= 0 &&
							kp_s0_s1[j].second >= hits[hits.size() - 2]){
							considered[j] = true;
						}
					}
				} else { /* s0_step > 0 */
					if (s0_step <= block_dist) {
						if (s1_step > 0) {
							if (s1_step <= block_dist) {
								hits.push_back(kp_s0_s1[j].first);
								hits.push_back(kp_s0_s1[j].first + K - 1);
								hits.push_back(kp_s0_s1[j].second);
								hits.push_back(kp_s0_s1[j].second + K - 1);
								considered[j] = true;
							} else {
								if (sat_dist_cri (l0, l1, hits,
										k_frame_min_perc, block_dist)) {
									frames.push_back(std::pair<int, lvec_t>
									(hits[hits.size()-1] - hits[2], hits));
								}
								hits.clear();
								break;
							}
						} // else if (s1_step <= 0) is captured in main loop
					} else {
						if (sat_dist_cri (l0, l1, hits,
								k_frame_min_perc, block_dist)) {
							frames.push_back(std::pair<int, lvec_t>
								(hits[hits.size()-1] - hits[2], hits));
						}
						hits.clear();
						break;
					}
				}
			} // for j
			if (hits.size() &&
				sat_dist_cri (l0, l1, hits, k_frame_min_perc, block_dist)){
				frames.push_back(std::pair<int, lvec_t>
					(hits[hits.size()-1] - hits[2], hits));
				hits.clear();
			}
		}// for i

		/* sort [len_hit] w.r.t. len */
		std::sort(frames.begin(), frames.end(), cmp_pair_ge());

		if (debug){ // debug print hits
			std::cout << "\n\tkframe coordinates: \n";
			for (int i = 0; i < frames.size(); ++ i) {
				std::cout << frames[i].first << "\n";
				for (int j = 0; j < frames[i].second.size(); ++ j){
					std::cout << frames[i].second[j] << "\t";
				}
				std::cout << "\n";
			}
			std::cout << "\n";
		}


		/* generate alignment for each kmer-frame, also for left & right
		 * sides of this frame. Concatenate three alignment to yield the
		 * longest alignment that meets the constraints. */

		ncbi::CNWAligner aligner;

		for (int i = 0; i < frames.size(); ++ i) {

			hits = frames[i].second;
			int hs = hits.size();

			/* identify boundaries for three alignments: kmer-frame,
			 *  left, right wrt the contig
			 */

			int l_min_hang = std::min(hits[0], hits[2]);
			int s0_l_end = hits[0] - l_min_hang,
				s1_l_end = hits[2] - l_min_hang;
			int r_min_hang =
				std::min(l0 - hits[hs - 3] - 1, l1 - hits[hs - 1] - 1);
			int s0_r_end = hits[hs - 3] + r_min_hang,
				s1_r_end = hits[hs - 1] + r_min_hang;

			/* adjust s0/1_left/right to account for ins/del outside kframe */
			int l_ext = (100 - min_iden) * l_min_hang/100;
			int r_ext = (100 - min_iden) * r_min_hang/100;
			if (s0_l_end == 0) s1_l_end = std::max(s1_l_end - l_ext, 0);
			else s0_l_end = std::max(s0_l_end - l_ext, 0);
			if (s0_r_end == l0 - 1)
				s1_r_end = std::min(l1 - 1, s1_r_end + r_ext);
			else s0_r_end = std::min(l0 - 1, s0_r_end + r_ext);

			/* calculate coordinates and hit vector for left */
			int l_from_0 = s0_l_end, l_to_0 = hits[1],
				l_from_1 = s1_l_end, l_to_1 = hits[3];
			lvec_t l_hits (4, 0);
			l_hits[0] = hits[0] - s0_l_end;
			l_hits[1] = hits[1] - s0_l_end;
			l_hits[2] = hits[2] - s1_l_end;
			l_hits[3] = hits[3] - s1_l_end;

			/* right */
			int r_from_0 = hits[hs - 4],
				r_to_0 = s0_r_end,
				r_from_1 = hits[hs - 2],
				r_to_1 = s1_r_end;
			lvec_t r_hits (4, 0);
			r_hits[0] = 0;
			r_hits[1] = hits[hs - 3] -	hits[hs - 4];
			r_hits[2] = 0;
			r_hits[3] = hits[hs - 1] - hits[hs - 2];

			/* kmer-frame */
			int k_from_0 = hits[0],
				k_to_0 = hits[hs - 3],
				k_from_1 = hits[2],
				k_to_1 = hits[hs - 1];
			for (int i = 0; i < hs; ++ i) {
				if (i % 4 == 0 || i % 4 == 1)
					hits[i] -= k_from_0;
				else hits[i] -= k_from_1;
			}

			/*** ALIGNMENTS ***/

			/* final alignment coordinates for s0, s1*/
			int s0_lb, s0_rb, s1_lb, s1_rb;

			/* align kframe */
			std::string k_align;
			nwalign(k_align, aligner, s0.substr(k_from_0, k_to_0 - k_from_0 + 1),
					s1.substr(k_from_1, k_to_1 - k_from_1 + 1), hits);

			if (debug) { // debug print
				std::cout << ">1\n" << s0.substr(k_from_0, k_to_0 - k_from_0 + 1) << "\n";
				std::cout << ">2\n" << s1.substr(k_from_1, k_to_1 - k_from_1 + 1) << "\n";
				std::cout << "> k_align \n" << k_align << "\n\n";
			}

			if (debug) {// debug
				if (align_identity (k_align, 3) < min_iden) {
					std::cout << "\t failed\n\n";
				}
			}

			/* if kmer-frame form good alignment */
			if (align_identity (k_align, 3) >= min_iden) {

				s0_lb = k_from_0; s0_rb = k_to_0;
				s1_lb = k_from_1; s1_rb = k_to_1;

				std::string l_align, r_align;

				if (k_from_0 != 0 && k_from_1 != 0) { // align "left"

					nwalign(l_align, aligner, s0.substr(l_from_0, l_to_0 - l_from_0 + 1),
						s1.substr(l_from_1, l_to_1 - l_from_1 + 1), l_hits);

					if (debug){ // debug print
						std::cout << ">1\n" << s0.substr(l_from_0, l_to_0 - l_from_0 + 1) <<"\n";
						std::cout << ">2\n" << s1.substr(l_from_1, l_to_1 - l_from_1 + 1) <<"\n";
						std::cout << "> l_align\n" << l_align << "\n";
					}

					int align_start = l_align.find_first_of('M');
					int l_hit_len = l_hits[1] - l_hits[0] + 1;
					if (align_start < (int) l_align.length() - l_hit_len) {
						int l_align_len = l_align.length() -
											l_hit_len - align_start;
						std::string l_align_region =
								l_align.substr(align_start, l_align_len);

						/* identify left max applicable region */
						std::reverse(l_align_region.begin(), l_align_region.end());
						get_max_constrained_aln_prefix (l_align_region,
								l_align_region, min_iden, 15);
						l_align_len = l_align_region.length ();
						if (l_align_len) {

							std::reverse(l_align_region.begin(), l_align_region.end());
							k_align = l_align_region + k_align;

							if (debug) {
								std::cout << "> l_align_region\n" << l_align_region << "\n\n";
							}
							// count # of chars involved in the l_align_region for
							// s0 and s1, then adjust left boundary
							int num_char_0 = 0, num_char_1 = 0;
							for (int l = 0; l < (int) l_align_len; ++ l){
								switch (l_align_region.at(l)){
								case 'M':
								case 'R':
									++ num_char_0;
									++ num_char_1;
									break;
								case 'I':
									++ num_char_1;
									break;
								case 'D':
									++ num_char_0;
									break;
								}
							}
							s0_lb -= num_char_0;
							s1_lb -= num_char_1;
						}

					} // if (align_start < (int) l_align.length() - l_hit_len)
				}

				if (k_to_0 != (int) s0.length() - 1  &&
						k_to_1 != (int) s1.length() - 1) { // align "right"

					nwalign(r_align, aligner,
							s0.substr(r_from_0, r_to_0 - r_from_0 + 1),
						s1.substr(r_from_1, r_to_1 - r_from_1 + 1), r_hits);

					if (debug) { // debug print
						std::cout << ">1\n" << s0.substr(r_from_0, r_to_0 - r_from_0 + 1) << "\n";
						std::cout << ">2\n" << s1.substr(r_from_1, r_to_1 - r_from_1 + 1) << "\n";
						std::cout << "> r_align\n" << r_align << "\n";
					}

					int align_end = r_align.find_last_of ('M');
					int r_hit_len = r_hits[1] - r_hits[0] + 1;
					if (align_end >= r_hit_len) {
						int r_align_len = align_end - r_hit_len + 1;
						std::string r_align_region =
								r_align.substr(r_hit_len, r_align_len);

						get_max_constrained_aln_prefix (r_align_region,
										r_align_region, min_iden, 15);

						r_align_len = r_align_region.length ();
						if (r_align_len) {

							if (debug) {
								std::cout << "> r_align_region\n" << r_align_region << "\n\n";
							}

							k_align += r_align_region;

							/* adjust right boundary */
							int num_char_0 = 0, num_char_1 = 0;
							for (int r = 0; r < r_align_len; ++ r){
								switch (r_align_region.at(r)){
								case 'M':
								case 'R':
									++ num_char_0;
									++ num_char_1;
									break;
								case 'I':
									++ num_char_1;
									break;
								case 'D':
									++ num_char_0;
									break;
								}
							}
							s0_rb += num_char_0;
							s1_rb += num_char_1;
						}

					} // if (align_end >= K)
				}

				if (k_align.length() >= min_overlap) {
					hits.clear();
					hits.push_back(s0_lb);
					hits.push_back(s0_rb);
					hits.push_back(s1_lb);
					hits.push_back(s1_rb);
					if (sat_dist_cri (l0, l1, hits, 90, max_overhang)){
						out = align_t (s0_lb, s0_rb, s1_lb, s1_rb, k_align);

						if (debug){
							std::cout << ">merged_alignment\n";
							std::cout << k_align << "\n\n";
						}
						return;
					} else {
						if (debug) {
							std::cout << "no valid alignment \n\n";
						}
					}
				}
			} // if (align_identity (k_align, 3) >= min_iden)
		} // for (int i = 0; i < ranges.size(); ++ i)
	} // kmer_anchor_ps_aligner

	/*
	 * return size of common suffix between two strings
	 */
	inline int common_pref_len (const std::string& s0, const std::string& s1) {
		int len = std::min(s0.length(), s1.length());
		int idx = 0;
		for (; idx < len; ++ idx) {
			if (toupper(s0.at(idx)) != toupper(s1.at(idx))) return idx;
		}
		return idx;
	}


	inline int common_suf_len (const std::string& s0, const std::string& s1) {
		int len0 = s0.length(), len1 = s1.length();
		int mlen = std::min(len0, len1);
		int idx = 0;
		for (; idx < mlen; ++ idx) {
			if (toupper (s0.at(len0 - idx - 1)) != toupper(s1.at (len1 - idx -1 )) )
				return idx;
		}
		return idx;
	}

	/* compare strings hamming distance less equal than */
	inline bool hdlet (const char* s1, const char* s2, int length, int d) {
		int hd = 0;
		for (int i = 0; i < length; ++ i) {
			if (s1[i] != s2[i]) ++ hd;
			if (hd > d) return false;
		}
		return true;
	}
	// compare IDs
	template <typename int_t>
	bool hdlet (int_t e1, int_t e2, int k, int d) {
		int_t diff = e1 ^ e2;
		int cnt = 0;
		while (diff) {
			if (diff & 0x3) ++ cnt;
			if (cnt > d) return false;
			diff >>= 2;
		}
		return true;
	}

}
#endif /* SEQ_PAIR_MANIP_HPP_ */

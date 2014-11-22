//========================================================================
// Project     : M-Vicuna
// Name        : seq_cmp.hpp
// Author      : Xiao Yang
// Created on  : Jun 5, 2013
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


#ifndef SEQ_CMP_HPP_
#define SEQ_CMP_HPP_

#include <tuple>
#include <map>
#include <vector>
#include <set>
#include <stdint.h>
#include <string>
#include <algorithm>
#include "seq_manip.hpp"
#include "../jaz/sequence_compare.hpp"

typedef std::pair<int, int> ipair_t;
typedef std::vector<uint32_t> uvec_t;
typedef std::set<uint32_t> uset_t;
typedef std::tuple<int, int, int, int> coord_t;

namespace xny{

	/** Function common_pref_len ()
	 * return length of common prefix between two strings
	 */
	inline int common_pref_len (const std::string& s0, const std::string& s1) {
		int len = std::min(s0.length(), s1.length());
		int idx = 0;
		for (; idx < len; ++ idx) {
			if (toupper(s0.at(idx)) != toupper(s1.at(idx))) return idx;
		}
		return idx;
	} // common_pref_len

	/** Function common_pref_len ()
	 * return length of common suffix between two strings
	 */
	inline int common_suf_len (const std::string& s0, const std::string& s1) {
		int len0 = s0.length(), len1 = s1.length();
		int mlen = std::min(len0, len1);
		int idx = 0;
		for (; idx < mlen; ++ idx) {
			if (toupper (s0.at(len0 - idx - 1)) != toupper(s1.at (len1 - idx -1 )) )
				return idx;
		}
		return idx;
	} // common_suf_len

	/** Function hd ()
	 *
	 * Calculate Hamming distance between two DNA strings
	 * N/n is considered wild char that doesnt contribute to distance incre
	 */
	inline int hd (const std::string& s0, const std::string& s1){
		int l = std::min ((int) s0.length(), (int) s1.length()), dist = 0;
		for (int i = 0; i < l; ++ i) {
			char c0 = std::toupper(s0.at(i)), c1 = std::toupper(s1.at(i));
			if (c0 != c1 && c0 != 'N' && c1 != 'N') {
				++ dist;
			}
		}
		return dist;
	}

	/** Function hdlet()
	 *
	 * Check if the Hamming distance between two DNA strings
	 * is <= maxhd; N/n doesn't contribute to distance increase;
	 *
	 * return -1 if >= maxhd otherwise return actual hd
	 */
	inline int hdlet (const std::string& s0, const std::string& s1, int maxhd){
		int l = std::min ((int) s0.length(), (int) s1.length()), dist = 0;
		for (int i = 0; i < l; ++ i) {
			char c0 = std::toupper(s0.at(i)), c1 = std::toupper(s1.at(i));
			if (c0 != c1 && c0 != 'N' && c1 != 'N') {
				++ dist;
				if (dist > maxhd) return -1;
			}
		}
		return dist;
	}

	/** Function hdlet ()
	 *
	 * Check if the Hamming distance between two DNA in bit format
	 * is <= maxhd; N/n doesn't contribute to distance increase;
	 *
	 * return -1 if >= maxhd otherwise return actual hd
	 */
	template <typename int_t>
	int hdlet (int_t e1, int_t e2, int maxhd) {
		int_t diff = e1 ^ e2;
		int dist = 0;
		while (diff) {
			if (diff & 0x3) ++ dist;
			if (dist > maxhd) return -1;
			diff >>= 2;
		}
		return dist;
	} // hdlet

	/**
	 * Functor implementation of suffix prefix gap free alignment between
	 * DNA strings s0 and s1
	 */
	class suffix_prefix_gap_free_aln {
	public:
		suffix_prefix_gap_free_aln (int cv, int pi, int oh):
			min_len_ (cv),	min_perc_identi_(pi), max_overhang_ (oh) { }

		/** Function: operator()
		 *  Compute alignment between s0 and s1.
		 *
		 *  Returns:
      	 *  4-tuple (aln start, end positions on s0 then on s2).
      	 */
		coord_t operator() (const std::string& s0,
				const std::string& s1) {

			int seed = 7;

			int l0 = s0.length(), l1 = s1.length();
			std::map<int, ipair_t> kID2pos;
			// dividing s0 into non-overlapping seed-mers
			// while neglecting the overhang in the suffix
			int num_sevn_mer = (l0 - max_overhang_)/7;
			int start = l0 - num_sevn_mer * seed - max_overhang_;
			char* addr = const_cast<char*> (s0.c_str());
			for (int i = 0; i < num_sevn_mer; ++ i) {
				uint32_t ID;
				if (xny::str2ID<uint32_t> (ID, addr + start, seed)) {
					kID2pos[ID] = ipair_t (start, -1);
				}
				start += seed;
			}
			// go through s1, identify common seed-mers shared with s0
			addr = const_cast<char*> (s1.c_str());
			for (int l = 0; l < l1 - seed + 1; ++ l) {
				uint32_t ID;
				if (xny::str2ID<uint32_t> (ID, addr + l, seed)) {
					std::map<int, ipair_t>::iterator it = kID2pos.find(ID);
					if (it != kID2pos.end()) {
						it->second.second = l;
					}
				}
			}
			// now anchored by each common seed, check validity of alignment
			//        /
			// -------
			//    ---------
			//   /
			std::map<int, ipair_t>::iterator it = kID2pos.begin();
			for (; it != kID2pos.end(); ++ it) {
				int p0_seed = it->second.first, p1_seed = it->second.second;
				if (p1_seed != -1) { // common seed found
					int lcmplen = std::min(p1_seed - max_overhang_,
							p0_seed - max_overhang_),
					rcmplen = std::min(l1 - p1_seed - seed - max_overhang_,
								l0 - p0_seed - seed - max_overhang_),
					maxhd = (lcmplen + rcmplen + seed)*(100 - min_perc_identi_)/100;

					// min_len_ criterion
					if (lcmplen + rcmplen + seed < min_len_) continue;
					// left side of seed
					int l_hd = 0, r_hd = 0;
					if (lcmplen > 0) {
						l_hd = hdlet (s0.substr(p0_seed - lcmplen, lcmplen),
								s1.substr(p1_seed - lcmplen, lcmplen), maxhd);
						if (l_hd == -1) continue;
					}
					if (rcmplen > 0) {
						r_hd = hdlet (s0.substr(p0_seed + seed, rcmplen),
								s1.substr(p1_seed + seed, rcmplen), maxhd);
						if (r_hd == -1) continue;
					}
					if (l_hd + r_hd <= maxhd) { // valid alignment
						return std::make_tuple (p0_seed - lcmplen,
								p0_seed + seed + rcmplen - 1, p1_seed - lcmplen,
								p1_seed + seed + rcmplen - 1);
					}
				}
			} // for
			return std::make_tuple(-1, -1, -1, -1);
		} // std::tuple<int, int, int, int> operator()

	private:
		int min_len_;
		int min_perc_identi_;
		int max_overhang_;
	}; // class suffix_prefix_gap_free_aln


	/**
	 * Functor implementation of alignment between DNA strings s0 and s1
	 * using the default 7-mers as seeds to anchor the alignment.
	 */
	class kmer_anchor_aln {
	// match block type --common substring between s0 and s1
	// s0_mb_start, s0_mb_end, s1_mb_start, s1_mb_end
	public:
		kmer_anchor_aln (int minl, int pi, int oh):
			min_len_(minl),	min_perc_identi_(pi), max_overhang_ (oh)
			{ seed_ = 7; }

		void set_seed (int kmer) {
			if (kmer < 3) {
				std::cout << "Warning kmer_anchor_aln: set_seed () the seed is"
						"ridiculously small, I'm resetting it to 7... \n";
				seed_ = 7;
			} else seed_ = kmer;
		} // set_seed

		/** Function: operator()
		 *  Compute alignment between s0 and s1 seeded by unique 7 mers
		 *
		 *  Returns: 5-tuple (aln start positions on s0, s0_end, start s1,
		 *  					 end_s1, alignment string wrt s0).
		 *  alignment string notations
      	 *  - M: match
      	 *  - R: mismatch
      	 *  - I: insertion
      	 *  - D: deletion
      	 *
      	 *  If no valid alignment found, return (-1, -1, -1, -1, "")
      	 */
		std::tuple<coord_t, std::string> operator () (const std::string& s0,
				const std::string& s1) {

			bool debug = true;


			int l0 = s0.length(), l1 = s1.length();

			std::map<int, ipair_t> seed2pos;

			// ----------------- s0 seeds   -----------------------------
			char* addr = const_cast<char*> (s0.c_str());
			uset_t dupl;
			for (int i = 0; i < l0 - seed_ + 1; ++ i) {
				uint32_t seedID;
				if (xny::str2ID<uint32_t> (seedID, addr + i, seed_)) {
					if (seed2pos.count (seedID)) dupl.insert (seedID);
					else seed2pos[seedID] = ipair_t (i, -1);
				}
			}
			// remove duplicates
			for (uset_t::iterator it = dupl.begin(); it != dupl.end(); ++ it){
				seed2pos.erase(*it);
			}

			// ----------------- s1 seeds   -----------------------------
			// go through s1, identify common seed-mers shared with s0
			dupl.clear();
			addr = const_cast<char*> (s1.c_str());
			for (int i = 0; i < l1 - seed_ + 1; ++ i) {
				uint32_t seedID;
				if (xny::str2ID<uint32_t> (seedID, addr + i, seed_)) {
					std::map<int, ipair_t>::iterator it = seed2pos.find(seedID);
					if (it != seed2pos.end()) {
						if (it->second.second != -1) dupl.insert(seedID);
						else it->second.second = i;
					}
				}
			}
			// remove duplicates
			for (uset_t::iterator it = dupl.begin(); it != dupl.end(); ++ it){
				seed2pos.erase(*it);
			}

			// debug print out seed2pos
			if (debug) {
				std::cout << "reads: \n" << s0 << "\n" << s1 << "\n\n";
				std::cout << "seed2pos:\n";
				std::map<int, ipair_t>::iterator it = seed2pos.begin();
				for (; it != seed2pos.end(); ++ it) {
					std::cout << xny::ID2Str(it->first, seed_) << "\t"
							<< it->second.first << "," << it->second.second <<"\n";
				}
				std::cout << "\n";
			}

			// ------- identify seed-frames and generate alignment ------

			std::vector<ipair_t> match_pos; // match positions in s0 and s1 sorted in
											// increasing order on s0
			std::map<int, ipair_t>::iterator it_s2p = seed2pos.begin();
			for (; it_s2p != seed2pos.end(); ++ it_s2p) {
				if (it_s2p->second.second != -1) {
					match_pos.push_back(it_s2p->second);
				}
			}
			std::sort(match_pos.begin(), match_pos.end());

			if (debug) {
				std::cout << "match_pos in s0 and s1:\n";
				for (int i = 0; i < (int) match_pos.size(); ++ i) {
					std::cout << match_pos[i].first << ","
							<< match_pos[i].second << "\t";
 				}
				std::cout << "\n";
			}

			if (match_pos.size() == 0) return std::make_tuple(
					std::make_tuple(-1,-1,-1,-1), "");

			// ---------------------------------------------------------

			// an empirical value to define the distance between two mbs
			// in a frame and an empirical value to define the min span for
			// a frame over the entire potential alignment region for it
			// to be further considered for alignment
			int mb_dist = std::min (80, 5 * seed_);
			int min_frame_span = 50;

			std::vector<bool> considered (match_pos.size(), false);

			for (int i = 0; i < (int) match_pos.size() - 1; ++ i) {

				if (considered[i]) continue;

				// frame consists of a list of mbs
				std::vector<coord_t> frame;

				frame.push_back(coord_t (match_pos[i].first,
					match_pos[i].first + seed_ - 1, match_pos[i].second,
					match_pos[i].second + seed_ - 1));

				for (int j = i + 1; j < (int) match_pos.size(); ++ j) {
					int prv_mb0_end = std::get<1> (frame.back());
					int step0 = match_pos[j].first + seed_ - 1
							- std::get<1> (frame.back()),
					    step1 = match_pos[j].second + seed_ - 1
					    		- std::get<3> (frame.back());
					if (match_pos[j].first <= prv_mb0_end) { // overlap on s0
						//extending the kmer
						if (step0 == step1) {
							std::get<1> (frame.back()) += step0;
							std::get<3> (frame.back()) += step0;
							considered[j] = true;
						}
					} else { // non-overlap on s0
						if (step0 <= mb_dist) {
							int prv_mb1_end = std::get<3> (frame.back());
							if (match_pos[j].second > prv_mb1_end) { // non-overlap on s1
								if (step1 <= mb_dist) {
									frame.push_back(
										std::make_tuple(match_pos[j].first,
										match_pos[j].first + seed_ - 1,
										match_pos[j].second,
										match_pos[j].second + seed_ - 1));
								} else {
									if (is_valid_frame (l0, l1, frame,
											min_frame_span, mb_dist)) {

										std::tuple<coord_t, std::string>
											result = align (	s0, s1, frame);

										if (is_valid_alignment(l0, l1, result)) {
											return result;
										}

										break;
									}
								}
							}
						} else {
							if (is_valid_frame (l0, l1, frame,
									min_frame_span, mb_dist)) {

								std::tuple<coord_t, std::string>
									result =  align (s0, s1, frame);
								if (is_valid_alignment(l0, l1, result)) {
									return result;
								}
								break;
							}
						}
					}

				} // for (int j = 1

				if (is_valid_frame(l0, l1, frame, min_frame_span, mb_dist)){

					std::tuple<coord_t, std::string>
						result =  align (s0, s1, frame);
					if (is_valid_alignment(l0, l1, result)) return result;
					break;
				}
			} // for (int i = 0

			return std::make_tuple(std::make_tuple(-1, -1, -1, -1), "");

		} // function: operator ()

	private:
		int min_len_;
		int min_perc_identi_;
		int max_overhang_;
		int seed_;

		bool is_valid_frame (int l0, int l1, const std::vector<coord_t>& frame,
				int min_frame_span, int frame_dist_to_boundary) {
			int sz = frame.size();
			if (!sz) return 0;
			int frame_length = std::max(
				std::get<1> (frame.back()) - std::get<0>(frame.front()) + 1,
				std::get<3> (frame.back()) - std::get<2>(frame.front()) + 1);
			int l_min_hang = std::min(std::get<0>(frame.front()),
					std::get<2>(frame.front())),
				r_min_hang =	 std::min(l0 - std::get<1> (frame.back()) - 1,
						l1 - std::get<3> (frame.back()) - 1);

			int span = 100 * frame_length / (frame_length + l_min_hang + r_min_hang);
			if (span >= min_frame_span && l_min_hang <= frame_dist_to_boundary
					&& r_min_hang <= frame_dist_to_boundary) {
				return true;
			}
			return false;
		} // is_valid_frame

		std::tuple<coord_t, std::string> align (const std::string& s0,
				const std::string& s1, const std::vector<coord_t>& frame) {

			int l0 = s0.length(), l1 = s1.length();
			// parameter from bwa-sw (http://bio-bwa.sourceforge.net/bwa.shtml)
			bio::global_alignment aligner (1, -3, -5, -2);

			//------------------ frame alignment -------------------------
			if (!frame.size()) {
				return std::make_tuple (std::make_tuple(-1,-1,-1,-1), "");
			}

			std::string alnstr (
				std::get<1>(frame[0]) - std::get<0>(frame[0]) + 1, 'M');

			for (unsigned int i = 1; i < frame.size(); ++ i) {
				// standard global alignment of strings b/t two mbs
				std::string sub0 = s0.substr(std::get<1>(frame[i-1]) + 1,
					std::get<0>(frame[i]) - std::get<1>(frame[i-1]) - 1),
				    sub1 = s1.substr(std::get<3>(frame[i-1]) + 1,
					std::get<2>(frame[i]) - std::get<3>(frame[i-1]) - 1);
				aligner (sub0, sub1);
				alnstr += aligner.path();

				alnstr += std::string (
					std::get<1>(frame[i]) - std::get<0>(frame[i]) + 1, 'M');
			}

			// ---------------- set alignment type -----------------------
			aligner.set_alignment_type(2);

			//------------------ left overhang alignment -----------------
			int l_min_hang = std::min(std::get<0>(frame[0]),
									  std::get<2> (frame[0]));
			int start0 = std::get<0> (frame[0]) - l_min_hang,
				start1 = std::get<2> (frame[0]) - l_min_hang,
				len0 = l_min_hang, len1 = l_min_hang;

			if (l_min_hang != 0) {
				// append extra max l_min_hang/2 length substring for aln
				if (start0 > 0) {
					start0 = std::max (0, start0 - l_min_hang/2);
					len0 = std::get<0> (frame[0]) - start0;
				} else if (start1 > 0) {
					start1 = std::max (0, start1 - l_min_hang/2);
					len1 = std::get<2>(frame[0]) - start1;
				}
				std::string substr0 = s0.substr(start0, len0),
						substr1 = s1.substr(start1, len1);
				std::reverse(substr0.begin(), substr0.end());
				std::reverse(substr1.begin(), substr1.end());

				aligner (substr0, substr1);

				std::string lalnstr = aligner.path();

				{ // debug print
					std::cout << "starts, lengths: " << start0 << ", " << len0 << ", "
							<< start1 << ", " << len1 << "\n";
					std::cout << "substrings: \n" << substr0 << "\n" << substr1 << "\n";
					std::cout << "lalnstr = " << lalnstr << "\n";
				}

				std::reverse(lalnstr.begin(), lalnstr.end());

				// stripping of gaps and mismatches M/I/R from the beginning
				bool stop_strip = false;
				int strip_len = 0, llen = lalnstr.length();
				for (int i = 0; i < llen; ++ i) {
					switch (lalnstr.at(i)) {
					case 'M':
						stop_strip = true;
						break;
					case 'I':
						++ start1;
						-- len1;
						++ strip_len;
						break;
					case 'R':
						++ start0;
						++ start1;
						-- len0;
						-- len1;
						++ strip_len;
						break;
					case 'D':
						++ start0;
						-- len0;
						++ strip_len;
						break;
					}
					if (stop_strip) break;
				}

				lalnstr = lalnstr.substr(strip_len, llen - strip_len);

				{ // debug print
					std::cout << "starts, lengths: " << start0 << ", " << len0 << ", "
							<< start1 << ", " << len1 << "\n";
					std::cout << "substrings: \n" << substr0 << "\n" << substr1 << "\n";
					std::cout << "lalnstr = " << lalnstr << "\n";
				}

				alnstr = lalnstr + alnstr;
			}

			//------------------ right overhang alignment ----------------
			int 	r_min_hang = std::min(l0 - std::get<1> (frame.back()) - 1,
						l1 - std::get<3> (frame.back()) - 1);

			int end0 = std::get<1> (frame.back()) + r_min_hang,
				end1 = std::get<3> (frame.back()) + r_min_hang;
			if (r_min_hang > 0) {

				len0 = r_min_hang;
				len1 = r_min_hang;

				// append extra max r_min_hang/2 length substring for aln
				if (end0 + 1 < l0) {
					end0 = std::min(l0 - 1, end0 + r_min_hang/2);
					len0 = end0 - std::get<1> (frame.back());
				} else if (end1 + 1 < l1) {
					end1 = std::min(l1 - 1, end1 + r_min_hang/2);
					len1 = end1 - std::get<3> (frame.back());
				}
				aligner (s0.substr(end0 - len0 + 1, len0),
						s1.substr(end1 - len1 + 1, len1));

				std::string ralnstr = aligner.path();

				{ // debug print
					std::cout << "ends, lengths: " << end0 << ", " << len0 << ", "
							<< end1 << ", " << len1 << "\n";
					std::cout << "substrings: \n" << s0.substr(end0 - len0 + 1, len0)
							<< "\n" << s1.substr(end1 - len1 + 1, len1) << "\n";
					std::cout << "ralnstr = " << ralnstr << "\n";
				}

				// stripping of gaps and mismatches M/I/R from the end
				int rlen = ralnstr.length(), strip_len = 0;
				bool stop_strip = false;
				for (int i = rlen - 1; i >= 0; -- i) {
					switch (ralnstr.at(i)) {
					case 'M':
						stop_strip = true;
						break;
					case 'I':
						-- len1;
						++ strip_len;
						-- end1;
						break;
					case 'R':
						-- len0;
						-- len1;
						++ strip_len;
						-- end0;
						-- end1;
						break;
					case 'D':
						-- len0;
						++ strip_len;
						-- end0;
						break;
					}
					if (stop_strip) break;
				} // for (int i

				ralnstr = ralnstr.substr(0, rlen - strip_len);
				{ // debug print
					std::cout << "starts, lengths: " << end0 << ", " << len0 << ", "
							<< end1 << ", " << len1 << "\n";
					std::cout << "substrings: \n" << s0.substr(end0, len0)
							<< "\n" << s1.substr(end1, len1) << "\n";
					std::cout << "ralnstr = " << ralnstr << "\n";
				}

				alnstr += ralnstr;
			}

			{
				std::cout << "final alignment:\n";
				std::cout << "start0, end0 = " << start0 << ", " << end0 << "\n";
				std::cout << "start1, end1 = " << start1 << ", " << end1 << "\n";
				std::cout << alnstr << "\n\n";
			}

			return std::make_tuple (
					std::make_tuple(start0, end0, start1, end1), alnstr);
		} // align

		/**	Given an initial suffix-prefix alignment between s0 and s1,
		 * check if any sub-alignment conforming the input constraint
		 */
		bool is_valid_alignment (int l0, int l1,
				std::tuple<coord_t, std::string>& result) {

			int start0 = std::get<0> (std::get<0> (result)),
				end0 = std::get<1> (std::get<0> (result)),
				start1 = std::get<2> (std::get<0> (result)),
				end1 = std::get<3> (std::get<0> (result));
			std::string aln = std::get<1> (result);

			int alnlen = aln.length(), Ms = 0, Rs = 0, Is = 0, Ds = 0;
			for (int i = 0; i < alnlen; ++ i) {
				switch (aln.at(i)) {
				case 'M':	++ Ms;	break;
				case 'R':	++ Rs;	break;
				case 'I':	++ Is;	break;
				case 'D':	++ Ds;	break;
				default:	 std::cout << "is_valid_alignment SC failed\n";
					exit(1);
				}
			}

			// identify alignment conforming criteria
			bool to_check_prefix = true, to_check_suffix = true;
			while (true) {
				if (100 * Ms/alnlen >= min_perc_identi_) {
					if (std::min (start0, start1) <= max_overhang_ &&
						std::min(l0 - end0 - 1, l1 - end1 - 1) <= max_overhang_ &&
						alnlen >= min_len_){

						std::get<0> (std::get<0> (result)) = start0;
						std::get<1> (std::get<0> (result)) = end0;
						std::get<2> (std::get<0> (result)) = start1;
						std::get<3> (std::get<0> (result)) = end1;
						std::get<1> (result) = aln;
						return true;

					} else break;
				} else {
					if (!to_check_prefix && !to_check_suffix) break;
					// strip off low similarity prefix or suffix until
					// unable to satisfy max_overhang_

					int global_sim = 100 * Ms /alnlen;

					//check prefix or suffix based on flag
					int prfx_M, prfx_R, prfx_I, prfx_D; // length, num of non Ms

					if (to_check_prefix) {
						if (!get_non_m_prefix (prfx_M, prfx_R,
								prfx_I, prfx_D, aln)) break;

						int prfx_len = prfx_M + prfx_R + prfx_I + prfx_D,
							prfx_sim = 100 * prfx_M / prfx_len;

						// if the similarity is lower compared to the global sim
						if (prfx_sim < global_sim) {

							// satisfy overhang criterion
							if (std::min (start0 + prfx_len - prfx_I,
								start1 + prfx_len - prfx_D) <= max_overhang_
								&&
								(alnlen - prfx_len) >= min_len_) {

								start0 += prfx_len - prfx_I;
								start1 += prfx_len - prfx_D;
								aln = aln.substr(prfx_len, alnlen - prfx_len);
								alnlen = aln.length();
								Ms -= prfx_M;

							} else to_check_prefix = false;
						} else to_check_prefix = false;
					} // if (to_check_prefix)

					if (to_check_suffix) {

						std::string rvs_aln = aln;
						std::reverse(rvs_aln.begin(), rvs_aln.end());

						if (!get_non_m_prefix(prfx_M, prfx_R,
								prfx_I, prfx_D, rvs_aln)) break;

						int prfx_len = prfx_M + prfx_R + prfx_I + prfx_D,
							prfx_sim = 100 * prfx_M / prfx_len;

						// if similarity is lower compared to the global
						if (prfx_sim < global_sim) {

							// satisfy overhang criterion
							if (std::min(l0 - end0 - 1 + (prfx_len - prfx_I),
								l1 - end1 - 1 + (prfx_len - prfx_D)) <= max_overhang_
								&& (alnlen - prfx_len) >= min_len_ ) {
								end0 -= (prfx_len - prfx_I);
								end1 -= (prfx_len - prfx_D);
								aln = aln.substr(0, alnlen - prfx_len);
								alnlen = aln.length();
								Ms -= prfx_M;
							} else to_check_suffix = false;
						} else to_check_suffix = false;
					} // if (to_check_suffix)


				} // if ... else

			} // while (true)
			return false;
		} // is_valid_alignment

		/**	Obtain the minimal prefix of alignment string that contains
		 * the first substring of the non-Ms, its length of the prefix and
		 * number of non-Ms within this prefix.
		 * e.g. MIIIRMM, returns 5, 1 as such a prefix = MIIIR
		 *
		 * Return true if existing such a prefix; false, otherwise
		 */
		bool get_non_m_prefix (int& Ms, int& Rs,
				int& Is, int& Ds, const std::string& alnstr){
			Ms = Rs = Is = Ds = 0;
			unsigned int first_non_M = alnstr.find_first_not_of('M');
			if (first_non_M == std::string::npos) return false;

			Ms += first_non_M;

			bool is_next_M_found = false;
			for (unsigned int i = first_non_M; i < alnstr.length(); ++ i){
				switch (alnstr.at(i)) {
				case 'M': is_next_M_found = true; break;
				case 'R': ++ Rs; break;
				case 'I': ++ Is;	break;
				case 'D': ++ Ds;	break;
				}
				if (is_next_M_found) break;
			}

			return is_next_M_found;
		} // get_non_m_prefix

	}; //class kmer_anchor_aln

} // namespace xny


#endif /* SEQ_CMP_HPP_ */

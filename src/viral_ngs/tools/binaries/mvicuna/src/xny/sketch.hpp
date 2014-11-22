//========================================================================
// Project     : M-Vicuna
// Name        : sketch.hpp
// Author      : Xiao Yang
// Created on  : Jul 2, 2013
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


#ifndef SKETCH_HPP_
#define SKETCH_HPP_

#include "seq_manip.hpp"

typedef std::pair<uint64_t, int> sketch_t;

namespace xny{
	/**
	 * Functor implementation of generating sketch list of a given string
	 */
	class sketch_list {
	public:
		/** constructor sketch_list ()
		 *
		 * is_fwd_only = true for protein sequence
		 */
		sketch_list (int k, bool is_fwd_only):
			k_ (k),	is_fwd_only_ (is_fwd_only) { }

		/** Function: operator() sketch_list
		 *
		 *  Compute sketch list for input [seq]
		 *	Return vector of (val, pos), where val is the hash value of
		 *	the alphabetically larger kmer starting at pos. If rvc strand
		 *	is used to obtain val, pos is set to be negative
		 *
      	 */
		template <typename T>
		std::vector<sketch_t> operator() (const std::string& seq, const T& hfunc) {

			std::vector<sketch_t> list;

			int num = seq.length() - k_ + 1;

			for (int i = 0; i < num; ++ i) {
				int pos = i;
				std::string s = seq.substr(i, k_);
				if (! is_fwd_only_) {
					std::string s_rvc = xny::get_rvc_str(s);
					if (s < s_rvc) {
						s = s_rvc;
						pos *= -1 ;
					}
				}
				list.push_back({hfunc (s), pos});
			}
			return list;
		} // operator()

	private:
		int k_;
		bool is_fwd_only_;
	}; // class sketch_list

	/**	Functor to compare sketch_t wrt the val field
	 */
	struct cmp_sketch {
		bool operator () (const sketch_t& lhs, const sketch_t& rhs) const {
			return lhs.first < rhs.first;
		}
	};

	/**
	 * Functor implementation of generating the super sketch given
	 * the list of sketch of a given string
	 */
	class super_sketch {
		public:
			/** constructor super_sketch ()
			 */
			super_sketch (int w): w_ (w) { }

			/** Function: operator() super_sketch
			 *
			 * return false if super_sketch doesn't exist
			 * super_sketch.first: ss value
			 * super_sketch.second: the position of the sketch on the string
			 *  that is the first sketch of the super_sketch
	      	 */
			template <typename T>
			bool operator() (sketch_t& super_sketch,
					const std::vector<sketch_t>& list, const T& hfunc) {

				int num = list.size() - w_ + 1;
				if (num <= 0) return false;

				std::vector<uint64_t> sketch_val_list;
				for (auto& x: list) sketch_val_list.push_back(x.first);

				std::vector<sketch_t> super_sketch_list (num);
			    unsigned int v_sz = w_ * sizeof(uint64_t);
			    for (unsigned int i = 0; i < num; ++ i) {
			    		const char* v = reinterpret_cast<const char*>(&sketch_val_list[i]);
			    		super_sketch_list[i].first = hfunc (v, v_sz);
			    		super_sketch_list[i].second = i;
			    }

			    std::sort(super_sketch_list.begin(), super_sketch_list.end(),
			    		cmp_sketch());

			    super_sketch.first = super_sketch_list.front().first;
			    int idx_min_ss = super_sketch_list.front().second;
			    super_sketch.second = list[idx_min_ss].second;
				return true;
			} // operator()

		private:
			int w_; // number of sketches to form a window
		}; // class super_sketch

} // namespace xny

#endif /* SKETCH_HPP_ */

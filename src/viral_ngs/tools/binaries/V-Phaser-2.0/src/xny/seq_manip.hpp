//============================================================================
// Name        : seq_manip.hpp
// Author      : Xiao Yang
// Created on  : Aug 3, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================
/*
std::vector<int> x(100);
std::vector<int> y(1000);
jaz::random_sample_n(y.begin(), y.end(), 100, x.begin());
x.clear();
jaz::random_sample_n(y.begin(), y.end(), 100, std::back_inserter(x));
template <Iter> void f(Iter first, Iter last) {
}
now, if Iter is RandomAccess Iterator you can write:

first[10] = ...;
unsigned int x = last - first;
http://www.sgi.com/tech/stl/stl_introduction.html
*/

#ifndef SEQ_MANIP_HPP_
#define SEQ_MANIP_HPP_

#include <iostream>
#include <string>
#include <ctype.h>
#include <iterator>
#include <cmath>
#include <vector>
#include <algorithm>

namespace xny{

	inline char bits2char(int value) {
		switch (value) {
        case 0:
            return 'a';
        case 1:
            return 'c';
        case 2:
            return 'g';
        case 3:
            return 't';
        default:
            return 'n';
		}
	}

	inline int char2bits (char c) {
		int cvalue = -1;
		switch (std::toupper(c)) {
			case 'A':
				cvalue = 0;
				break;
			case 'C':
				cvalue = 1;
				break;
			case 'G':
				cvalue = 2;
				break;
			case 'T':
				cvalue = 3;
				break;
		}
		if (cvalue != -1) return cvalue;
		else return -1;
	} // char2bits

	inline void remove_gap (std::string& seq) {
		std::string::iterator it_end =
				std::remove(seq.begin(), seq.end(), '-');
		seq.resize(it_end - seq.begin());
	}
	/*
	 *	DNA string to type T \in {uint32_t, uint64_t}
	 *	Consider alphabet ACGTacgt only
	 */
	template <typename T>
	bool str2ID (T& ID, char* addr, int len) {
		ID = 0;
		for (int i = 0; i < len; ++ i){
			int c = char2bits(addr[i]);
			if (c == -1) return false;
		ID  = (ID << 2 | c);
		}
		return true;
	} // str2ID

	/* convert uint32_t or uint64_t type to DNA sequence
	 */
	template<typename in_t>
	std::string ID2Str(in_t ID, int len){
	    std::string kmer = "";
	    for (int i = 0; i < len; ++ i){
	        int last = (ID & 0x3);
	        char c;
	        switch (last){
	            case 0: c = 'a';
	            break;
	            case 1: c = 'c';
	            break;
	            case 2: c = 'g';
	            break;
	            case 3: c = 't';
	            break;
	        }
	        kmer += c;
	        ID = ID >> 2;
	    }
	    std::reverse(kmer.begin(), kmer.end());
	    return kmer;
	}

	/* @brief	Convert a DNA string in bit format to its reverse
	 *			complementary bit format
	 */
	template <typename sz_t>
	sz_t get_rvc_bits (sz_t bitSeq, int len){
	    sz_t rv = 0;
	    for (int i = 0; i < len; ++ i){
	        int cvalue = bitSeq & 0x3;
	        switch (cvalue){
	            case 0:
	                rv = (rv << 2) | 0x3;
	                break;
	            case 1:
	                rv = (rv << 2) | 0x2;
	                break;
	            case 2:
	                rv = (rv << 2) | 0x1;
	                break;
	            case 3:
	                rv = rv << 2;
	                break;
	        }
	        bitSeq >>= 2;
	    }
	    return rv;
	}

	/* check if the char is a proper nucleotide (acgtACGT) or not */
	inline bool is_nt (char c) {
		switch (std::toupper(c)) {
		case 'A':
		case 'C':
		case 'G':
		case 'T':
			return true;
		default: return false;
		}
	}// is_nt

	/* check if the oligo consists of acgt/ACGT only */
	inline bool is_nt_string (const std::string& oligo) {
		for (int i = 0; i < (int) oligo.length(); ++ i) {
			if (!is_nt(oligo.at(i))) return false;
		}
		return true;
	}

	inline char get_rvc_base (char in) {
		switch (in){
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'a':
			return 't';
		case 't':
			return 'a';
		case 'c':
			return 'g';
		case 'g':
			return 'c';
		default:
			return 'N';
		}
	}
	/* Convert an input DNA string (will not be modified)
	 * to its reverse complementary; the result is stored in another string.
	 * All non-ACGTacgt characters will be returned as 'n'
	 */
	inline std::string get_rvc_str (std::string fwd) {
		std::reverse(fwd.begin(), fwd.end());
		for (int i = 0; i < (int) fwd.length(); ++ i) {
			fwd[i] = get_rvc_base (fwd[i]);
		}
		return fwd;
	} // get_rvc_str

	/* @brief 	The input DNA string [fwd] is converted to its rvc, the
	 * 			result is stored in [fwd].
	 * @note		All non-ACGTacgt characters will be returned as 'n'
	 */
	inline void rvc_str (std::string& fwd) {
		std::reverse(fwd.begin(), fwd.end());
		for (int i = 0; i < (int) fwd.length(); ++ i){
			switch (fwd.at(i)){
			case 'A':
				fwd[i] = 'T';
				break;
			case 'T':
				fwd[i] = 'A';
				break;
			case 'C':
				fwd[i] = 'G';
				break;
			case 'G':
				fwd[i] = 'C';
				break;
			case 'a':
				fwd[i] = 't';
				break;
			case 't':
				fwd[i] = 'a';
				break;
			case 'c':
				fwd[i] = 'g';
				break;
			case 'g':
				fwd[i] = 'c';
				break;
			default:
				fwd[i] = 'N';
				break;
			}
		}
	} // get_rvc_str
	/*
	 * extract all kmers (including duplicated ones) in bit form from a
	 * DNA sequence, including forward and/or reverse complementary
	 * strands according to the [code]:
	 * code =
	 * 		 0: alphabetically smaller kmers are retained
	 * 		 1: alphabetically larger kmers are retained
	 * 		 2: kmers are retained in the order of
	 * 		 	fwd0, rv0, fwd1, rv1 ...
	 * 		 3: only fwd strand fw0, fw1 ... fw_i
	 * 		 4: only rv strand rv0 rv1 ...rv_i, i is w.r.t the fwd strand
	 *
	 * Usage: xny::get_bitkmer<std::back_insert_iterator<Container_t>, out_t>()
	 */
	template <typename outputIter, typename out_t>
	void get_bitkmer (const std::string& seq,
			outputIter out, int k, int code) {

		if (k > (int) seq.length()) return;
		// lowest 2k bits = 1
		out_t MaskLowerKbits = (out_t) pow(2, 2*k) - 1;
		std::vector<out_t> MaskTop2bits(4, 0);
		for (out_t i = 1; i < 4; ++ i){
			MaskTop2bits[i] = (i << (2*(k - 1)));
		}
		std::reverse(MaskTop2bits.begin(), MaskTop2bits.end()); // t g c a

		char* addr = const_cast<char*> (seq.c_str());
		int len = (int) seq.length() - k + 1;

		out_t ID = 0, rcID = 0;
		int j = 0;
		bool flag = false; //true: the previous kmer contains no Ns
		while (j < len){
			if (flag){
				int c = char2bits (*(addr + j + k -1));
				if (c == -1) {
					j += k;
					flag = false;
				}
				else {
					ID  = ((ID << 2 | c) & MaskLowerKbits);
					rcID = (rcID >> 2);
					rcID |= MaskTop2bits[c];

					switch (code) {
					case 0:
						if (rcID < ID) { *out = rcID; ++ out; }
						else  {	*out = ID; ++ out; }
						break;
					case 1:
						if (rcID < ID) { *out = ID; ++ out; }
						else  {	*out = rcID; ++ out; }
						break;
					case 2:
						*out = ID; ++ out;
						*out = rcID; ++ out;
						break;
					case 3:
						*out = ID; ++ out;
						break;
					case 4:
						*out = rcID; ++ out;
						break;
					default:
						std::cout << "In func get_bitkmer, [code] out of range\n";
						exit(1);
					}

					++ j;
				}
			}
			else{
				if (str2ID (ID, addr + j, k)) {
					rcID = get_rvc_bits <out_t> (ID, k);
					flag = true;

					switch (code) {
					case 0:
						if (rcID < ID) { *out = rcID; ++ out; }
						else  {	*out = ID; ++ out; }
						break;
					case 1:
						if (rcID < ID) { *out = ID; ++ out; }
						else  {	*out = rcID; ++ out; }
						break;
					case 2:
						*out = ID; ++ out;
						*out = rcID; ++ out;
						break;
					case 3:
						*out = ID; ++ out;
						break;
					case 4:
						*out = rcID; ++ out;
						break;
					default:
						std::cout << "In func get_bitkmer, [code] out of range\n";
						exit(1);
					}

				}
				++ j;
			}
		}
	} // get_bitkmer

	/********************************************************************
	 * extract all kmers (including duplicated ones) in string form from a
	 * DNA sequence, including both forward and reverse complementary
	 * strands, depends on the [code], different kmers are retained:
	 * code =
	 * 		 0: alphabetically smaller kmers are retained
	 * 		 1: alphabetically larger kmers are retained
	 * 		 2: both kmers are retained, in the order of
	 * 		 	fwd0, rv0, fwd1, rv1 ...
	 * 		 3: only fwd strand fw0, fw1 ... fw_i
	 * 		 4: only rv strand rv0 rv1 ...rv_i, i is w.r.t the fwd strand
	 *
	 * 	Usage:
	 * 	xny::get_strkmer(std::string, std::back_inserter(Container_t), int, int);
	 ********************************************************************/
	template <typename outputIter>
	void get_strkmer (const std::string& seq, outputIter out, int k, int code) {

		if (k > seq.length()) return; // kmer is longer than the read

		std::string kmer, rcKmer;
		kmer = seq.substr(0, k);
		rcKmer = get_rvc_str (kmer);
		std::transform(kmer.begin(), kmer.end(), kmer.begin(), toupper);
		std::transform(rcKmer.begin(), rcKmer.end(), rcKmer.begin(), toupper);

		switch (code) {
			case 0:
				if (kmer.compare(rcKmer) <= 0) { *out = kmer; ++ out; }
				else  {	*out = rcKmer; ++ out; }
				break;
			case 1:
				if (kmer.compare(rcKmer) > 0) { *out = kmer; ++ out; }
				else  {	*out = rcKmer; ++ out; }
				break;
			case 2:
				*out = kmer; ++ out;
				*out = rcKmer; ++ out;
				break;
			case 3:
				*out = kmer; ++ out;
				break;
			case 4:
				*out = rcKmer; ++ out;
				break;
			default:
				std::cout << "In func get_strkmer, [code] out of range\n";
				exit(1);
		}
		int j = k;
		while (j < seq.length()){
			kmer.erase(0, 1);
			kmer = kmer + (char) toupper(seq.at(j));

			std::string tmp (1, std::toupper(get_rvc_base (seq.at(j))));
			rcKmer = tmp + rcKmer.substr(0, kmer.size() - 1);

			//rcKmer[k - 1] = '\0'; // remove the last position
			//rcKmer = (char) toupper( get_rvc_str(seq.substr(j, 1)).at(0) ) + rcKmer;

			switch (code) {
				case 0:
					if (kmer.compare(rcKmer) <= 0) { *out = kmer; ++ out; }
					else  {	*out = rcKmer; ++ out; }
					break;
				case 1:
					if (kmer.compare(rcKmer) > 0) { *out = kmer; ++ out; }
					else  {	*out = rcKmer; ++ out; }
					break;
				case 2:
					*out = kmer; ++ out;
					*out = rcKmer; ++ out;
					break;
				case 3:
					*out = kmer; ++ out;
					break;
				case 4:
					*out = rcKmer; ++ out;
					break;
				default:
					std::cout << "In func get_strkmer, [code] out of range\n";
					exit(1);
			}

			++ j;
		}
	} //get_stringkmer


	/********************************************************************
	 * @brief	Display a uint64 int in an arbitrary-base format
	 * @example	convBase (10, 2) ==> 1010
	 *******************************************************************/
	inline std::string convBase(unsigned long v, long base)
	{
		std::string digits = "0123456789abcdef";
		std::string result;
		if((base < 2) || (base > 16)) {
			result = "convBase: error -- base out of range [2, 16].\n";
		}
		else {
			do {
				result = digits[v % base] + result;
				v /= base;
			}
			while(v);
		}
		return result;
	}

	/********************************************************************
	 * @brief	Given a DNA sequence and a seed, denoted by a sequence
	 * 			of '0's and '1's, e.g., 100101, extract all seed window
	 * 			(including duplicated ones) in bit form, where positions
	 * 			of '0's are masked.
	 * 			Seed windows are retained according to [code]
	 * 	code =
	 * 		 0: alphabetically smaller seeds are retained
	 * 		 1: alphabetically larger seeds are retained
	 * 		 3: both seeds are retained, in the order of
	 * 		 	fwd0, rv0, fwd1, rv1 ...
	 * 		 	e.g., S = acggt, and seed = 1001, two windows are considered:
	 * 		 	--->    --->
	 * 		 	acgg		cggt
	 *          tgcc		gcca
	 *          <---    <---
	 *			after masking, we have a--g, c--t, c--t, a--g will be
	 *			retained. (in fact, since the mask is '0', the actual
	 *			4 resulting seeds are: aaag, caat, caat, aaag, where '-'
	 *			is replaced by 'a')
	 * 		 4: only fwd strand fw0, fw1 ... fw_i
	 * 		 5: only rv strand rv0 rv1 ...rv_i, i is w.r.t the fwd strand
	 *
	 * @ example usage:
	 * 			std::vector<std::pair<uint64_t, int> > seededkmers;
	 * 			std::vector<bool> seed (6, true);
	 * 			seed[1] = seed[4] = 0;
	 * 			xny::get_bit_gapped_seeds (std::back_inserter(seededkmers),
	 * 									   seed, "ACGTTTATG", 2);
	 *
	 ********************************************************************/
	template <typename outputIter>
	void get_bit_gapped_seeds (outputIter out, const std::vector<bool>& seed,
			const std::string& seq, int code) {

		int seedsz = seed.size();
		if (seedsz > 32) {
			std::cout << "[WARNING] get_str_seeded_kmer: seed size "
					<< seedsz << " is > 32, seeding not applicable\n";
			return;
		} else if (seedsz > seq.length()) {
			/*
			std::cout << "[WARNING] get_str_seeded_kmer: seed size "
					<< seedsz << " is > read length " << seq.length()
					<< ", seeding not applicable\n";
			*/
			return;
		}

		/* convert the seed into 64 bit format and
		 * generate a mask with lower 2*seedsz bits to be 1*/
		uint64_t bit_seed = 0;
		uint64_t MaskLowerKbits = 3;
		for (int i = 0; i < seedsz; ++ i) {
			if (seed[i]) bit_seed = (bit_seed << 2 | 3);
			else bit_seed <<= 2;
			MaskLowerKbits = (MaskLowerKbits << 2 | 3);
		}

		std::vector<uint64_t> MaskTop2bits(4, 0);
		for (uint64_t i = 1; i < 4; ++ i){
			MaskTop2bits[i] = (i << (2*(seedsz - 1)));
		}
		std::reverse(MaskTop2bits.begin(), MaskTop2bits.end()); // t g c a

		char* addr = const_cast<char*> (seq.c_str());
		int len = (int) seq.length() - seedsz + 1;

		uint64_t ID = 0, rcID = 0;
		int j = 0;
		bool flag = false; //true: the previous seed does not contain 'N'
		while (j < len){
			if (flag){
				int c = char2bits (*(addr + j + seedsz - 1));
				if (c == -1) {
					j += seedsz;
					flag = false;
				} else {
					ID  = ((ID << 2 | c) & MaskLowerKbits);
					rcID = (rcID >> 2);
					rcID |= MaskTop2bits[c];

					switch (code) {
					case 0:
						if (rcID < ID) {
							*out = std::pair<uint64_t, int>(
									(rcID & bit_seed), j);
							++ out;
						}
						else  {
							*out = std::pair<uint64_t, int> (
									(ID & bit_seed), j);
							++ out;
						}
						break;
					case 1:
						if (rcID < ID) {
							*out =  std::pair<uint64_t, int>(
									(ID & bit_seed), j);
							++ out;
						}
						else  {
							*out =  std::pair<uint64_t, int>(
									(rcID & bit_seed), j);
							++ out;
						}
						break;
					case 2:
						*out = std::pair<uint64_t, int>(
								(ID & bit_seed), j);
						++ out;
						*out = std::pair<uint64_t, int>(
								(rcID & bit_seed), j);
						++ out;
						break;
					case 3:
						*out = std::pair<uint64_t, int>(
								(ID & bit_seed), j);
						++ out;
						break;
					case 4:
						*out = std::pair<uint64_t, int>(
								(rcID & bit_seed), j);
						++ out;
						break;
					default:
						std::cout << "get_bit_gapped_seeds: [code] out of range\n";
						exit(1);
					}

					++ j;
				}
			}
			else{
				if (str2ID (ID, addr + j, seedsz)) {

					rcID = get_rvc_bits <uint64_t> (ID, seedsz);
					flag = true;

					switch (code) {
					case 0:
						if (rcID < ID) {
							*out = std::pair<uint64_t, int>(
									(rcID & bit_seed), j);
							++ out;
						}
						else  {
							*out = std::pair<uint64_t, int> (
									(ID & bit_seed), j);
							++ out;
						}
						break;
					case 1:
						if (rcID < ID) {
							*out =  std::pair<uint64_t, int>(
									(ID & bit_seed), j);
							++ out;
						}
						else  {
							*out =  std::pair<uint64_t, int>(
									(rcID & bit_seed), j);
							++ out;
						}
						break;
					case 2:
						*out = std::pair<uint64_t, int>(
								(ID & bit_seed), j);
						++ out;
						*out = std::pair<uint64_t, int>(
								(rcID & bit_seed), j);
						++ out;
						break;
					case 3:
						*out = std::pair<uint64_t, int>(
								(ID & bit_seed), j);
						++ out;
						break;
					case 4:
						*out = std::pair<uint64_t, int>(
								(rcID & bit_seed), j);
						++ out;
						break;
					default:
						std::cout << "get_bit_gapped_seeds: [code] out of range\n";
						exit(1);
					} // switch
				} // if
				++ j;
			} // else
		}
	} // get_bit_gapped_seeds

} // namespace xny


#endif /* DNASEQ_MANIP_HPP_ */

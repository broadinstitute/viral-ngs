//============================================================================
// Project     : Diversifier
// Name        : Util.h
// Author      : Xiao Yang
// Created on  : Jul 14, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#ifndef XUTIL_H_
#define XUTIL_H_


#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <omp.h>
#include <iomanip>
#include <tuple>
#include <stdlib.h> // rand
#include <limits.h>
#include "jaz/fasta_file.hpp"

#if defined (_MSC_VER)
  #include "xny/getWinTime.hpp"
#else
  #include <sys/time.h>
#endif

typedef std::vector<int> ivec_t;
typedef std::vector<int64_t> i64vec_t;
typedef std::vector<i64vec_t> ii64vec_t;
typedef std::vector<char> cvec_t;
typedef std::vector<bool> bvec_t;
typedef std::vector<ivec_t> iivec_t;
typedef std::vector<uint32_t> uvec_t;
typedef std::vector<uvec_t> uuvec_t;
typedef std::vector<std::string> strvec_t;
typedef std::pair<int, int> ipair_t;
typedef std::pair<bool, bool> bpair_t;
typedef std::pair<uint32_t, uint32_t> upair_t;
typedef std::pair<uint8_t, uint8_t> u8pair_t;
typedef std::pair<std::string, std::string> strpair_t;
typedef std::map<int, int> imap_t;
typedef std::map<int64_t, int64_t> i64map_t;
typedef std::map<uint32_t, uint32_t> umap_t;
typedef std::set<int> iset_t;
typedef std::set<uint32_t> uset_t;
typedef std::set<std::string> strset_t;
typedef std::tuple<std::string, std::string, std::string> fqtuple_t;

// matching coordinates of subseqs between two strings s0, s1 (s0_start, s0_end, s1_start, s1_end)
typedef std::tuple<int, int, int, int> coord_t;



inline double get_time() {
      timeval t;
      gettimeofday(&t, 0);
      return t.tv_sec + (0.000001 * t.tv_usec);
}

inline void print_time (const std::string& msg, double& timing){
    double cur_time = get_time();
    std::cout << "\n" << msg << " *** " << (cur_time - timing)/60 << " mins ***\n\n";
    timing = cur_time;
}

inline void abording (const std::string& msg) {
	std::cout << "\n[EXIT]: " << msg << "\n";
	exit(1);
}

inline void warning (const std::string& msg) {
	std::cout << "\t[WARNING]: " << msg << "\n";
}

template <typename T>
void print_1dvec (const std::vector<T>& myVec, std::ostream* out) {
	std::cout << "vector content:\n";
	for (int i = 0; i < (int) myVec.size(); ++ i) {
		(*out) << "\t(" << i << "," << myVec[i] << ")";
		if ((i!=0) && (i % 6 == 0)) (*out) << "\n";
	}
	(*out) << "\n";
}

// print to file
inline void print_2dvec (const iivec_t& myVec, std::ostream* out){
	iivec_t::const_iterator it1 (myVec.begin());
	std::vector<int>::const_iterator it2;
	int index = 0;
	for (; it1 != myVec.end(); ++ it1){
		(*out) << (index ++) << ": ";
		for (it2 = it1->begin(); it2 != it1->end(); ++ it2){
			(*out) << (*it2) << "\t";
		}
		(*out) << "\n";
	}
} // print_2dvec

/** Split an input string, s, into a list of strings, out,
 * 	according to the separator, c.
 */
template <typename Iter>
void split(char c, const std::string& s, Iter out) {
    unsigned int pos = 0;
    for (unsigned int i = 0; i < s.length(); ++i) {
	  if (s[i] == c) {
	      if (i - pos > 0) *(out++) = s.substr(pos, i - pos);
	      pos = i + 1;
	  }
    }
    *(out++) = s.substr(pos, s.length() - pos);
} // split


// ------------------ union find structure --------------------------
/** Function find(): part of union find algorithm
 *
 * identify the root of leaf [elem] in the union find structure [uf],
 * where the [uf] is updated during the process
 */
template <typename T>
T uf_find(T elem, std::vector<T>& uf){

    if (elem == uf[elem])
        return elem;
    else {
    		uf[elem] = uf_find (uf[elem], uf);
        return uf[elem];
    }
}

/** Function clsfind (): part of union find algorithm
 *
 * identify the root of leaf [elem] in the union find structure [uf] is
 * not modified
 */
template <typename T>
T uf_clsfind (T elem, const std::vector<T>& uf){
    while (elem != uf[elem]) elem = uf[elem];
    return elem;
} // clsfind

//----- generate final cluster { clusterID --> fragment IDs } ------
/** Funtion uf_generate_cls ();
 *	Given the union find structure uf, generate the final
 *	cluster in 2d vector, where each 1-d vector stores elements belonging
 *	to the same cluster. The elems are sorted in the increasing order
 */

template <typename T>
void uf_generate_cls (iivec_t& clusters, const std::vector<T>& uf) {
	std::map<int, ivec_t> root_elems;
	std::map<int, ivec_t>::iterator it;
	int sz = uf.size();
	for (int i = 0; i < sz; ++ i) {
		int root_i = uf_clsfind<T> (i, uf);
		it = root_elems.find (root_i);
		if (it != root_elems.end()) it->second.push_back(i);
		else root_elems[root_i] = ivec_t (1, i);
	}

	// go through the map and produce clusters in sorted vector format
	for (it = root_elems.begin(); it != root_elems.end(); ++ it) {
		if (it->second.size() > 1) {
			std::sort (it->second.begin(), it->second.end());
			clusters.push_back(it->second);
		}
	} // for (it
} // uf_generate_cls
#endif /* XUTIL_H_ */

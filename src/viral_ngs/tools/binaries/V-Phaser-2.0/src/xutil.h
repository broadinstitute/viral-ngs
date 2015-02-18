//============================================================================
// Project     : VPhaser 2.0
// Name        : xutil.h
// Author      : Xiao Yang
// Created on  : Jul 14, 2011
// Version     :
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//============================================================================


#ifndef XUTIL_H_
#define XUTIL_H_

// #include <sys/resource.h>  Getrusage()

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
#include <stdint.h>
//#if defined (_MSC_VER)
//  #include "xny/getWinTime.hpp"
//#else
  #include <sys/time.h>
//#endif

#include "xny/file_manip.hpp"
#include "jaz/string_add.hpp"

typedef std::vector<int> ivec_t;
typedef std::vector<char> cvec_t;
typedef std::vector<double> dvec_t;
typedef std::vector<ivec_t> iivec_t;
typedef std::vector<dvec_t> ddvec_t;
typedef std::vector<uint32_t> uvec_t;
typedef std::vector<uvec_t> uuvec_t;
typedef std::vector<std::string> strvec_t;
typedef std::pair<int, int> ipair_t;
typedef std::pair<bool, bool> bpair_t;
typedef std::pair<uint32_t, uint32_t> upair_t;
typedef std::pair<uint8_t, uint8_t> u8pair_t;
typedef std::pair<std::string, std::string> strpair_t;
typedef std::map<int, int> imap_t;
typedef std::map<std::string, int> strimap_t;
typedef std::map<uint32_t, uint32_t> umap_t;
typedef std::set<int> iset_t;
typedef std::set<std::string> strset_t;


#include <sstream>



template <typename T>
inline std::string to_string (const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}

template <typename T>
inline T string_to (const std::string& st)
{
	T t;
	std::stringstream (st) >> t;
	return t;
}

inline void abording (const std::string& msg) {
	std::cout << "[EXIT]: " << msg << "\n";
	exit(1);
}

inline void warning(const std::string& msg) {
	std::cout << "[WARNING]: " << msg << "\n\n";
}

inline void print_strvec (const std::string& delim,
		const std::vector<std::string>& list) {
	for (int i = 0; i < (int) list.size(); ++ i) {
		std::cout << delim << list[i] << "\n";
	}
	std::cout << "\n";
} //print_strvec

/* @brief	Get the index of an element in a linear array.
 * 			Return -1 if not found.
 */
template <typename T>
int get_index (const T& elem, const std::vector<T>& array) {
	for (int i = 0; i < (int) array.size(); ++ i) {
		if (array[i] == elem) return i;
	}
	return -1;
} // get_index


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

template <typename T>
void debug_print_vec(const std::vector<T>& vec) {
	for (int i = 0; i < (int) vec.size(); ++ i) {
		std::cout << vec[i] << "\t";
	}
	std::cout << "\n";
}

#endif /* XUTIL_H_ */

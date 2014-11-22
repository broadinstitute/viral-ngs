//========================================================================
// Project     : M-Vicuna
// Name        : DuplRm.h
// Author      : Xiao Yang
// Created on  : Jun 25, 2013
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


#ifndef DUPLRM_H_
#define DUPLRM_H_

#include "xutil.h"
#include "Parameter.h"
#include "ReadBioFile.h"
#include "xny/file_manip.hpp"
#include "xny/seq_cmp.hpp"
#include "xny/sketch.hpp"
#include "jaz/fastx_iterator.hpp"
#include "jaz/hash.hpp"

void debug_print_fragments (const ivec_t& fragIDs, const std::string& fq,
		const std::string& fq2);

void duplicate_removal (const strvec_t& ifqs, const drm_t& drm, int w,
		int w2, xny::low_complexity& lc, int batch, bool silent);

void clustering_via_ss (ivec_t& uf_clst, const ii64vec_t& list_seeds,
	const std::string& f, const std::string& f2, int batch,
	xny::sketch_list& slistgen, xny::super_sketch& ssgen, int max_mismatch,
	bool silent);

void get_super_sketches (std::vector<sketch_t>& super_sketches,
	const std::string& f, const std::string& f2, xny::sketch_list& slistgen,
	xny::super_sketch& ssgen, jaz::murmur264& hashfunc, int batch, bool silent);

void generate_super_sketches (std::vector<sketch_t>& super_sketches,
		const strvec_t& pairs, xny::sketch_list& slistgen,
		xny::super_sketch ssgen, jaz::murmur264& hashfunc);

void get_seeds (ii64vec_t& list_seeds, const std::string& f,
		const std::string& f2, int seed_len, int batch, bool silent);

void generate_seeds (ii64vec_t& list_seeds, const strvec_t& pairs,
	int seed_len);

void get_seeds_per_fragment (i64vec_t& seeds, const std::string& frag,
	 int seed_len);

void clustering_via_seeds (ivec_t& uf_clst, ii64vec_t& list_seeds,
		int num_seed, int max_mismatch, bool silent);

void validate_clusters (ivec_t& uf_clst, const iivec_t& init_clusters,
	 const ii64vec_t& list_seeds, int max_mismatch, int max_cls_sz);

void make_cluster (iivec_t& clusters, const ii64vec_t& list_seeds,
		const ivec_t& init_cluster, int max_mismatch, const ivec_t& uf_clst);

void update_uf (ivec_t& uf_clst, const ivec_t& clusters) ;

bool is_similar (const i64vec_t& s0, const i64vec_t& s1, int max_mismatch);

void clean_dupl_frag (const std::string& ifq, const std::string& ifq2,
	std::ofstream& ofhfq, std::ofstream& ofhfq2, const iset_t& duplIDs,
	xny::low_complexity& lc, int batch);

void check_low_complexity (iset_t& low_complex_frag, int start_fragID,
		const std::vector<fqtuple_t>& pairs, xny::low_complexity& lc);

/* sort seeds wrt the idx-th element */
struct cmp_seed{
public:
	cmp_seed (int index): idx_(index) {}
	bool operator () (const i64vec_t& lhs, const i64vec_t& rhs)
	const {	return lhs[idx_] < rhs[idx_]; }
private:
	int idx_;
};


#endif /* DUPLRM_H_ */

//========================================================================
// Project     : VariantCaller
// Name        : pileup.h
// Author      : Xiao Yang
// Created on  : Sep 25, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================


#ifndef PILEUP_H_
#define PILEUP_H_

#include "stat_func.h"
#include "aln_manip.h"

struct pileup_profile_t {
	double lamda_S, lamda_L; // lamda for SNP and LP column
	dvec_t p_S, p_L; 		 // probability vector
	strimap_t SProfile, LProfile; // actual profile
	pileup_profile_t () {
		lamda_S = 0.0;
		lamda_L = 0.0;
	}
};

struct pileup_list_t {
	iset_t S;
	iset_t L;
};

typedef std::map<int, strset_t> var_t; // refpos --> var_str0, var_str1...

void pileup (var_t& nuvars, pileup_list_t& nulist, pileup_list_t& rmlist,
		const pileup_list_t& pileuplist, const var_t& gvars,
		const col_t& col, const std::vector<jeb_t>& jeb,
		const GlobalParam& gParam, int iter);

void pileup_profiling (pileup_profile_t& profile_info, const col_t& col,
		const std::vector<jeb_t>& jeb, const GlobalParam& gParam);

bool call_pileup_variant (strset_t& vars, int threshold,
		const strimap_t& profile);

void add_var_entry (var_t& lhs, int pos, const strset_t& types);

double get_pe (int bk_numer, int bk_denom, const double& prior_pe);

#endif /* PILEUP_H_ */

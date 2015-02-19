//========================================================================
// Project     : VariantCaller
// Name        : phase.h
// Author      : Xiao Yang
// Created on  : Sep 25, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================


#ifndef PHASE_H_
#define PHASE_H_

#include "pileup.h"

struct phase_profile_t {
	dvec_t p_SS, p_SL, p_LS, p_LL; // probabilities for the column for poisson
	double lamda_SS, lamda_SL, lamda_LS, lamda_LL;
	strimap_t SS, SL, LS, LL;
	phase_profile_t () {
		lamda_SS = 0.0, lamda_SL = 0.0,
	    lamda_LS = 0.0, lamda_LL = 0.0;
	}
};

struct phase_nb_t {
	std::map<int, iset_t> SS;
	std::map<int, iset_t> SL;
	std::map<int, iset_t> LS;
	std::map<int, iset_t> LL;
};

void phase (var_t& nuvars, phase_nb_t& nu_nb, phase_nb_t& rm_nb,
		const phase_nb_t& phase_nb, const pileup_list_t& pileuplist,
		const var_t& gvars, const col_t& lhs, const col_t& rhs,
		const std::vector<jeb_t>& jeb, const GlobalParam& gParam, int iter) ;

void phase_profiling (phase_profile_t& profile, const col_t& lhs,
		const col_t& rhs, const std::vector<jeb_t>& jeb,
		const GlobalParam& gParam);

void add_to_phaselist (std::map<int, iset_t>& list, int p0, int p1);

bool call_phase_variant (strset_t& l_vars, strset_t& r_vars,	int l_pos,
	int r_pos, int threshold, const strimap_t& profile, const var_t& gvars);

bool is_in_phase (const std::map<int, iset_t>& list, int p0, int p1) ;

bool is_var_exist (int pos, const std::string& var, const var_t& gvars);

#endif /* PHASE_H_ */

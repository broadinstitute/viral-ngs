//========================================================================
// Project     : VariantCaller
// Name        : em.h
// Author      : Xiao Yang
// Created on  : Apr 21, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================


#ifndef EM_H_
#define EM_H_

#include "api/BamReader.h"
#include "xutil.h"
#include "Parameter.h"
#include "xny/EnumParser.hpp"
#include "bam_manip.h"
#include "aln_manip.h"
#include "stat_func.h"
#include "pileup.h"
#include "phase.h"

using namespace BamTools;

/* the following two struct is used for gathering variant information for
 * output to file
 */
struct varinfo_t {
	bool is_pass_strd_test;
	bool is_pass_strd_fdr_test;
	bool is_snp;		//variant type snp or lp
	int ref_pos;
	std::string var;
	double p_value;
	std::map<std::string, ipair_t> profile_strd;
	//varinfo_t () {
	//	p_value = 0.0;
	//	is_snp = true;
	//}
	varinfo_t (bool issnp, int pos, const std::string& str,
		const double& pval, const std::map<std::string, ipair_t>& prfl){
			is_snp = issnp;
			ref_pos = pos;
			var = str;
			p_value = pval;
			profile_strd = prfl;
			is_pass_strd_fdr_test = true;
			is_pass_strd_test = true;
	};
};

struct cmp_varinfo_pval {
	bool operator () (const varinfo_t& lhs, const varinfo_t& rhs) {
		return lhs.p_value < rhs.p_value;
	}
};

struct cmp_varinfo_refpos {
	bool operator () (const varinfo_t& lhs, const varinfo_t& rhs) {
		return lhs.ref_pos < rhs.ref_pos;
	}
};

template<typename T>
struct cmp_p {
	bool operator () (const std::pair<T, int>& lhs,
			const std::pair<T, int>& rhs) {
		return lhs.first < rhs.first;
	}
};

/*
struct varprofile_t {
	bool is_snp;
	int ref_pos;
	std::map<std::string, ipair_t> profile_strd; //polymor--> (#fwd, #rv)
	varprofile_t (bool issnp, int pos, const std::map<std::string, ipair_t>& p):
		is_snp (issnp), ref_pos (pos), profile_strd (p) {}
}; */

void EM (GlobalParam& gParam, const Parameter& myPara);

bool infer_variant (var_t& gvars, pileup_list_t& pileuplist,
		phase_nb_t& phase_nb, const std::vector<col_t>& cols,
		int maxPileupIdx, const std::vector<jeb_t>& jeb,
		const Parameter& myPara, const GlobalParam& gParam, int iter);

void merge_phase_nb (std::map<int, iset_t>& to,
		const std::map<int, iset_t>& from);

void subtract_phase_nb (std::map<int, iset_t>& lhs,
		const std::map<int, iset_t>& rhs);

// ------------------         ------------------------
void update_variants (std::vector<jeb_t>& jeb,
		std::map<int, bpair_t>& called,	Bucket& bkinfo, 	const var_t& gvars,
		const strvec_t& alnfiles, const Parameter& myPara);

bool is_contain_vars (const strset_t& strvars);

void clean_jeb (std::vector<jeb_t>& jeb, std::map<int, bpair_t>& called,
		Bucket& bkinfo,	const var_t& gvars,	const std::vector<col_t>& cols);

void get_consensus (std::string& cons_p, std::string& cons_s, const col_t& col);

ipair_t decre_jeb (std::vector<jeb_t>& jeb, const std::string& cons,
		const std::vector<entry_t>& entries, bool is_poly);

// ----------- output variants ----------------------
void finalize_variant (const var_t& gvars, const std::string& refName,
		const strvec_t& alnfiles, const Parameter& myPara);

void analyze_vars (std::vector<varinfo_t>& list_varInfo, const double& alpha,
	const std::vector<col_t>& cols, const std::vector<strset_t>& list_vars);
void strdbias_pval (std::vector<varinfo_t>& list_varInfo);
void benjamini_hchberg_fdr (std::vector<varinfo_t>& list_varInfo,
		const double& alpha);

void write_var (int& polycnt, int& scnt, std::ofstream& oHandle,
		const col_t& col, var_t::const_iterator it_v);
void output_vars (const std::vector<varinfo_t>& list_varInfo,
		const std::string& dir, const std::string& refName);
// ----------- debug print ----------------------
void debug_print_phasing_info (const phase_profile_t& profile,
		int l_pos, const strset_t& l_vars, int r_pos,
		const strset_t& r_vars, const var_t& local_vars);
void debug_print_var (const var_t& vars);
void debug_print_profile (const strimap_t& profile);
#endif /* EM_H_ */

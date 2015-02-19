//========================================================================
// Project     : VariantCaller
// Name        : pileup.cpp
// Author      : Xiao Yang
// Created on  : Sep 25, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================

#include "pileup.h"


/* @brief	For bases, differentiate by (base, #occur); For deletions,
 * 			differentiate by (D, #); for insertions, differentiate by
 * 			(I, ins_str)
 *
 * @param	[nuvars]: vars that were not previously called in [gvars]
 */
void pileup (var_t& nuvars, pileup_list_t& nulist, pileup_list_t& rmlist,
		const pileup_list_t& pileuplist, const var_t& gvars,
		const col_t& col, const std::vector<jeb_t>& jeb,
		const GlobalParam& gParam, int iter) {

	//{
	//	std::cout << "pos = " << col.ref_pos << "\n";
 	//}
	int threshold = INT_MAX;

 	// ---------- generate profile and calculate lamda -------------
	pileup_profile_t profile_info;
	strset_t vars;
	if (iter == 0) { // first iteration pop up the pic

		pileup_profiling (profile_info, col, jeb, gParam);

		if (profile_info.lamda_L) {
			threshold = getThreshold (profile_info.lamda_L,
					profile_info.p_L, gParam.bonferroni);
			// ----------------- analyze profile and call vars ------------------
			if (call_pileup_variant (vars, threshold, profile_info.LProfile))	{
				nulist.L.insert(col.ref_pos);
			}
		}

		if (profile_info.lamda_S) {
			threshold = getThreshold (profile_info.lamda_S,
					profile_info.p_S, gParam.bonferroni);
			if (call_pileup_variant (vars, threshold, profile_info.SProfile)) {
				nulist.S.insert(col.ref_pos);
			}
		}

		add_var_entry (nuvars, col.ref_pos, vars);

		{
		//	debug_print_var (nuvars);
		}
	} else { // not the first iteration
		if (pileuplist.L.count(col.ref_pos) || pileuplist.S.count(col.ref_pos)) {

			pileup_profiling (profile_info, col, jeb, gParam);

			if (profile_info.lamda_L && pileuplist.L.count(col.ref_pos)) {
				threshold = getThreshold (profile_info.lamda_L,
						profile_info.p_L, gParam.bonferroni);
				if (! call_pileup_variant (vars, threshold, profile_info.LProfile)){
					rmlist.L.insert(col.ref_pos);
				}
			}

			if (profile_info.lamda_S && pileuplist.S.count(col.ref_pos)) {
				threshold = getThreshold (profile_info.lamda_S,
						profile_info.p_S, gParam.bonferroni);
				if (! call_pileup_variant (vars, threshold, profile_info.SProfile)) {
					rmlist.S.insert(col.ref_pos);
				}
			}

			// identify variants that were not yet in [gvars]
			var_t::const_iterator it_gvar = gvars.find(col.ref_pos);
			if (it_gvar != gvars.end()) {
				strset_t diffs;
				strset_t::const_iterator it_s = vars.begin();
				for (; it_s != vars.end(); ++ it_s){
					if (!it_gvar->second.count(*it_s)) diffs.insert(*it_s);
				}
				if (diffs.size()) add_var_entry(nuvars, col.ref_pos, diffs);
			} else add_var_entry (nuvars, col.ref_pos, vars);

			{
				//debug_print_var (nuvars);
			}

		} // if

	} // else
} // pileup


/* @brief	For an alignment column
 * 		- generate cons_type [SProfile] for SNP entries;
 * 		- generate [LProfile] for LP entries i.e. Dx or Ix.
 * 		- Compute lamda values [lamda_S/L], and raw probabilities in [p_S/L]
 */
void pileup_profiling (pileup_profile_t& profile_info, const col_t& col,
		const std::vector<jeb_t>& jeb, const GlobalParam& gParam) {

	int num_entry = col.entries.size();
	bool is_ins_col = (col.ref_pos % 2 == 1) ? true : false;

	for (int i = 0; i < num_entry; ++ i) {
		if (col.entries[i].cons_type.at(0) == 'D'
				|| col.entries[i].cons_type.at(0) == 'I') {
			add_to_profile (profile_info.LProfile, col.entries[i].cons_type);
		} else {
			if (is_ins_col)	add_to_profile (profile_info.LProfile, "d");
			else {
				add_to_profile (profile_info.LProfile, "i");
				add_to_profile (profile_info.SProfile, col.entries[i].cons_type);
			}
		}

		if (jeb[col.entries[i].eb_index].lp_sum != 0) { //

			//double tmp = (jeb[col.entries[i].eb_index].lp_err + 1.0)/
			//		(jeb[col.entries[i].eb_index].lp_sum + 1.0);
			//tmp = (tmp == 1.0) ? 0.5 : tmp;
			int idx = col.entries[i].eb_index;
			// constant prior
			//double tmp = (jeb[idx].lp_err + gParam.bkinfo.get_prior_L())/
			//			 (jeb[idx].lp_sum + gParam.bkinfo.average_L());
			// weighted prior
			//double tmp = (jeb[idx].lp_err * gParam.bkinfo.average_L() +
			//		jeb[idx].lp_sum * gParam.bkinfo.get_prior_L())/
			//		(2 * jeb[idx].lp_sum * gParam.bkinfo.average_L());

			//lamda_lp += std::min(tmp, gParam.lenpoly_pe); // over calling
			double pe = get_pe (jeb[idx].lp_err, jeb[idx].lp_sum,
				1.0 * gParam.bkinfo.get_prior_L()/ gParam.bkinfo.average_L());

			profile_info.lamda_L += pe;		// under calling
			profile_info.p_L.push_back(pe);
		}

		if (! is_ins_col && col.entries[i].cons_type.at(0) != 'D'
				&& jeb[col.entries[i].eb_index].snp_sum != 0) { // SC

			//double tmp = (jeb[col.entries[i].eb_index].snp_err + 1.0)/
			//			(jeb[col.entries[i].eb_index].snp_sum + 1.0);
			//tmp = (tmp == 1.0) ? 0.5 : tmp;
			int idx = col.entries[i].eb_index;
			//double tmp = (jeb[idx].snp_err + gParam.bkinfo.get_prior_S())/
			//			(jeb[idx].snp_sum + gParam.bkinfo.average_S());

			// weighted piror
			//double tmp = (jeb[idx].snp_err * gParam.bkinfo.average_S() +
			//		jeb[idx].snp_sum * gParam.bkinfo.get_prior_S())/
			//		(2 * jeb[idx].snp_sum * gParam.bkinfo.average_S());
			//lamda_snp += gParam.pe; // uniform
			//lamda_snp += std::min(gParam.pe, tmp); // over-calling

			double pe = get_pe (jeb[idx].snp_err, jeb[idx].snp_sum,
				1.0 * gParam.bkinfo.get_prior_S()/ gParam.bkinfo.average_S());
			profile_info.lamda_S += pe;
			profile_info.p_S.push_back(pe);
		}

	}// for (int i = 0
} // pileup_profiling

/* @brief	Identify all variants (including the consensus) that pass
 * 			the significant threshold; return true if this column has
 * 			candidate variants that should be checked further.
 */
bool call_pileup_variant (strset_t& vars, int threshold,
		const strimap_t& profile){

	bool is_to_further_check = false;
	strimap_t::const_iterator it_p;
	for (it_p = profile.begin(); it_p != profile.end(); ++ it_p){
		if (it_p->second >= threshold) {
			vars.insert(it_p->first);
		} else if (it_p->second > 1) {

			//std::cout << it_p->first << "\t" << it_p->second << "\n";
			is_to_further_check = true;
		}
	} // for
	return is_to_further_check;
} // call_pileup_variant

/* @brief	Add variants of a particular ref position [pos] to existing
 * 			variant list [lhs]
 */
void add_var_entry (var_t& lhs, int pos, const strset_t& types) {
	if (types.size()) {
		var_t::iterator it = lhs.find(pos);
		if (it != lhs.end()) {
			it->second.insert(types.begin(), types.end());
		} else lhs[pos] = types;
	}
} // add_var_entry

/* @brief	Calculate PE for each jeb entry
 *
 */
double get_pe (int bk_numer, int bk_denom, const double& prior_pe){
	double min_cnt = 20.0;
	if (bk_denom < min_cnt) {
		return (bk_numer + min_cnt * prior_pe) / (bk_denom + min_cnt);
	} else {
		if (bk_numer == bk_denom) return bk_numer/(bk_denom + 1.0);
		else return (bk_numer + 0.5)/(bk_denom + 0.5);
	}
} //get_pe

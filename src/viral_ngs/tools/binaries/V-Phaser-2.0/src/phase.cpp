//========================================================================
// Project     : VariantCaller
// Name        : phase.cpp
// Author      : Xiao Yang
// Created on  : Sep 25, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================

#include "phase.h"


/* @brief
 *
 */
void phase (var_t& nuvars, phase_nb_t& nu_nb, phase_nb_t& rm_nb,
	const phase_nb_t& phase_nb, const pileup_list_t& pileuplist,
	const var_t& gvars, const col_t& lhs, const col_t& rhs,
	const std::vector<jeb_t>& jeb, const GlobalParam& gParam, int iter) {

	int bc = gParam.bonferroni *  (gParam.bonferroni - 1) / 2;
	//------ generate profile, calculate lamda, and stores P vector -----
	phase_profile_t profile_info;

	if (iter == 1) { // calculate phaselist
		phase_profiling (profile_info, lhs, rhs, jeb, gParam);
		strset_t l_nuvars, r_nuvars;

		if (pileuplist.L.count(lhs.ref_pos) &&
			pileuplist.L.count(rhs.ref_pos) && profile_info.lamda_LL) {
			int threshold = getThreshold (profile_info.lamda_LL,
					profile_info.p_LL, bc);
			if (call_phase_variant (l_nuvars, r_nuvars, lhs.ref_pos,
					rhs.ref_pos, threshold, profile_info.LL, gvars)){
				add_to_phaselist (nu_nb.LL, lhs.ref_pos, rhs.ref_pos);
			}
		}

		if (pileuplist.L.count(lhs.ref_pos) &&
			pileuplist.S.count(rhs.ref_pos) && profile_info.lamda_LS){

			int threshold = getThreshold (profile_info.lamda_LS,
					profile_info.p_LS, bc);
			if (call_phase_variant (l_nuvars, r_nuvars, lhs.ref_pos,
					rhs.ref_pos, threshold, profile_info.LS, gvars)){
				add_to_phaselist (nu_nb.LS, lhs.ref_pos, rhs.ref_pos);
			}
		}

		if (pileuplist.S.count(lhs.ref_pos) &&
			pileuplist.S.count(rhs.ref_pos) && profile_info.lamda_SS){

			int threshold = getThreshold (profile_info.lamda_SS,
					profile_info.p_SS, bc);
			if (call_phase_variant (l_nuvars, r_nuvars, lhs.ref_pos,
					rhs.ref_pos, threshold, profile_info.SS, gvars)){
				add_to_phaselist (nu_nb.SS, lhs.ref_pos, rhs.ref_pos);
			}
		}

		if (pileuplist.S.count(lhs.ref_pos) &&
			pileuplist.L.count(rhs.ref_pos) && profile_info.lamda_SL){

			int threshold = getThreshold (profile_info.lamda_SL,
					profile_info.p_SL, bc);
			if (call_phase_variant (l_nuvars, r_nuvars, lhs.ref_pos,
					rhs.ref_pos, threshold, profile_info.SL, gvars)){
				add_to_phaselist (nu_nb.SL, lhs.ref_pos, rhs.ref_pos);
			}
		}

		 add_var_entry (nuvars, lhs.ref_pos, l_nuvars);
		 add_var_entry (nuvars, rhs.ref_pos, r_nuvars);

	} else { //check phaselist first
		bool is_phase_LL = is_in_phase (phase_nb.LL, lhs.ref_pos, rhs.ref_pos),
			is_phase_LS = is_in_phase (phase_nb.LS, lhs.ref_pos, rhs.ref_pos),
			is_phase_SS = is_in_phase (phase_nb.SS, lhs.ref_pos, rhs.ref_pos),
			is_phase_SL = is_in_phase (phase_nb.SL, lhs.ref_pos, rhs.ref_pos);

		if (is_phase_LL || is_phase_LS || is_phase_SS || is_phase_SL) {

			phase_profiling (profile_info, lhs, rhs, jeb, gParam);
			strset_t l_nuvars, r_nuvars;

			if (is_phase_LL && profile_info.lamda_LL) {
				int threshold = getThreshold (profile_info.lamda_LL,
											profile_info.p_LL, bc);
				if (!call_phase_variant (l_nuvars, r_nuvars, lhs.ref_pos,
						rhs.ref_pos, threshold, profile_info.LL, gvars)){
					add_to_phaselist (rm_nb.LL, lhs.ref_pos, rhs.ref_pos);
				}
			}

			if (is_phase_LS && profile_info.lamda_LS){
				int threshold = getThreshold (profile_info.lamda_LS,
											profile_info.p_LS, bc);
				if (!call_phase_variant (l_nuvars, r_nuvars, lhs.ref_pos,
						rhs.ref_pos, threshold, profile_info.LS, gvars)){
					add_to_phaselist (rm_nb.LS, lhs.ref_pos, rhs.ref_pos);
				}
			}

			if (is_phase_SS && profile_info.lamda_SS){
				int threshold = getThreshold (profile_info.lamda_SS,
											profile_info.p_SS, bc);
				if (!call_phase_variant (l_nuvars, r_nuvars, lhs.ref_pos,
						rhs.ref_pos, threshold, profile_info.SS, gvars)){
					add_to_phaselist (rm_nb.SS, lhs.ref_pos, rhs.ref_pos);
				}
			}

			if (is_phase_SL && profile_info.lamda_SL){
				int threshold = getThreshold (profile_info.lamda_SL,
											profile_info.p_SL, bc);
				if (!call_phase_variant (l_nuvars, r_nuvars, lhs.ref_pos,
						rhs.ref_pos, threshold, profile_info.SL, gvars)){
					add_to_phaselist (rm_nb.SL, lhs.ref_pos, rhs.ref_pos);
				}
			}

			 add_var_entry (nuvars, lhs.ref_pos, l_nuvars);
			 add_var_entry (nuvars, rhs.ref_pos, r_nuvars);
		} // if (is_phase_LL ||

	} // else

} // phase

/* @brief	Generate profile for 4 different combination between two
 * 			columns, and calculate lamda values correspondingly.
 */
void phase_profiling (phase_profile_t& prflinfo, const col_t& lhs,
	const col_t& rhs, const std::vector<jeb_t>& jeb, const GlobalParam& gParam) {

	bool l_is_ins = (lhs.ref_pos % 2 == 1) ? true : false,
		 r_is_ins = (rhs.ref_pos % 2 == 1) ? true : false;

	// identify pair of entries on the same read/pair
	int l_sz = lhs.entries.size(), r_sz = rhs.entries.size();
	int l = 0, r = 0;
	while (l < l_sz && r < r_sz) {
		if (lhs.entries[l].rID == rhs.entries[r].rID) {
			std::string ltype = lhs.entries[l].cons_type,
				 rtype = rhs.entries[r].cons_type;
			std::string jtype;
			if (ltype.at(0) == 'D' || ltype.at(0) == 'I') { //
				if (rtype.at(0) == 'D' || rtype.at(0) == 'I') { //PP
					jtype += ltype + "," + rtype;
					add_to_profile (prflinfo.LL, jtype);
				} else {	//
					if (r_is_ins) { // PP(Iins)
						jtype += ltype + ",d";
						add_to_profile (prflinfo.LL, jtype);
					} else { // PS
						jtype +=  ltype + "," + rtype;
						add_to_profile (prflinfo.LS, jtype);
						jtype.clear();
						jtype += ltype + ",i";
						add_to_profile (prflinfo.LL, jtype);
					}
				}
			} else {
				if (rtype.at(0) == 'D' || rtype.at(0) == 'I') { //SP
					if (l_is_ins) {
						jtype += "d," + rtype;
						add_to_profile (prflinfo.LL, jtype);
					} else {
						jtype += "i," + rtype;
						add_to_profile (prflinfo.LL, jtype);
						jtype.clear();
						jtype += ltype + "," + rtype;
						add_to_profile (prflinfo.SL, jtype);
					}
				} else { //SS
					if (l_is_ins) {
						if (r_is_ins) { // Ins Ins
							add_to_profile (prflinfo.LL, "d,d"); // d d
						} else { // Ins Non-ins
							jtype += "d," + rtype;
							add_to_profile (prflinfo.LS, jtype);
							add_to_profile (prflinfo.LL, "d,i"); // d i
						}
					} else { // non-ins ins
						if (r_is_ins) {
							jtype += ltype + ",d";
							add_to_profile (prflinfo.SL, jtype);
							add_to_profile (prflinfo.LL, "i,d"); // i d
						} else { // non-ins non-ins
							jtype += ltype + "," + rtype;
							add_to_profile (prflinfo.SS, jtype);
							add_to_profile (prflinfo.LL, "i,i"); // i i
						}
					}
				} // else SS
			}
			// calculate lamdas
			int idx_l = lhs.entries[l].eb_index,
				idx_r = rhs.entries[r].eb_index;
			/*
			double l_p = (jeb[idx_l].lp_sum == 0) ? 0 :
				(jeb[idx_l].lp_err + 1.0)/(jeb[idx_l].lp_sum + 1.0),
			  l_s = (jeb[idx_l].snp_sum == 0) ? 0 :
				(jeb[idx_l].snp_err + 1.0)/(jeb[idx_l].snp_sum + 1.0),
			  r_p = (jeb[idx_r].lp_sum == 0) ? 0 :
				(jeb[idx_r].lp_err + 1.0)/(jeb[idx_r].lp_sum + 1.0),
			  r_s = (jeb[idx_r].snp_sum == 0) ? 0 :
				(jeb[idx_r].snp_err + 1.0)/(jeb[idx_r].snp_sum + 1.0);
			l_p = (l_p == 1.0) ? 0.5 : l_p;
			l_s = (l_s == 1.0) ? 0.5 : l_s;
			r_p = (r_p == 1.0) ? 0.5 : r_p;
			r_s = (r_s == 1.0) ? 0.5 : r_s;
			*/
			// constant prior
			/*
			double l_p = (jeb[idx_l].lp_sum == 0) ? 0 :
				(jeb[idx_l].lp_err + gParam.bkinfo.get_prior_L())/
				(jeb[idx_l].lp_sum + gParam.bkinfo.average_L()),
			  l_s = (jeb[idx_l].snp_sum == 0) ? 0 :
				(jeb[idx_l].snp_err + gParam.bkinfo.get_prior_S())/
				(jeb[idx_l].snp_sum + gParam.bkinfo.average_S()),
			  r_p = (jeb[idx_r].lp_sum == 0) ? 0 :
				(jeb[idx_r].lp_err + gParam.bkinfo.get_prior_L())/
				(jeb[idx_r].lp_sum + gParam.bkinfo.average_L()),
			  r_s = (jeb[idx_r].snp_sum == 0) ? 0 :
				(jeb[idx_r].snp_err + gParam.bkinfo.get_prior_S())/
				(jeb[idx_r].snp_sum + gParam.bkinfo.average_S());
			*/
			// weighted prior
			/*
			double l_p = (jeb[idx_l].lp_sum == 0) ? 0 :
				(jeb[idx_l].lp_err * gParam.bkinfo.average_L() +
				gParam.bkinfo.get_prior_L() * jeb[idx_l].lp_sum)/
				(2 * jeb[idx_l].lp_sum * gParam.bkinfo.average_L()),
			  l_s = (jeb[idx_l].snp_sum == 0) ? 0 :
				(jeb[idx_l].snp_err * gParam.bkinfo.average_S() +
			  	 gParam.bkinfo.get_prior_S() * jeb[idx_l].snp_sum)/
				(2 * jeb[idx_l].snp_sum * gParam.bkinfo.average_S()),
			  r_p = (jeb[idx_r].lp_sum == 0) ? 0 :
				(jeb[idx_r].lp_err * gParam.bkinfo.average_L() +
				 gParam.bkinfo.get_prior_L() * jeb[idx_r].lp_sum)/
				(2 * jeb[idx_r].lp_sum * gParam.bkinfo.average_L()),
			  r_s = (jeb[idx_r].snp_sum == 0) ? 0 :
				(jeb[idx_r].snp_err * gParam.bkinfo.average_S() +
				gParam.bkinfo.get_prior_S() * jeb[idx_r].snp_sum)/
				(2 * jeb[idx_r].snp_sum * gParam.bkinfo.average_S());
			*/
			double l_p = (jeb[idx_l].lp_sum == 0) ? 0 : get_pe (
				jeb[idx_l].lp_err, jeb[idx_l].lp_sum, 1.0 *
				gParam.bkinfo.get_prior_L()/ gParam.bkinfo.average_L()),
			  l_s = (jeb[idx_l].snp_sum == 0) ? 0 : get_pe (
				jeb[idx_l].snp_err, jeb[idx_l].snp_sum, 1.0 *
				gParam.bkinfo.get_prior_S()/gParam.bkinfo.average_S()),
			  r_p = (jeb[idx_r].lp_sum == 0) ? 0 : get_pe (
				jeb[idx_r].lp_err, jeb[idx_r].lp_sum, 1.0 *
				 gParam.bkinfo.get_prior_L()/gParam.bkinfo.average_L()),
			  r_s = (jeb[idx_r].snp_sum == 0) ? 0 : get_pe (
				jeb[idx_r].snp_err, jeb[idx_r].snp_sum, 1.0 *
				gParam.bkinfo.get_prior_S() /gParam.bkinfo.average_S());

			if (l_p) {
				if (r_p) { // ------ PP ------
					//prflinfo.lamda_pp +=
					//	std::min(l_p, gParam.lenpoly_pe) * std::min(r_p, gParam.lenpoly_pe);
					double tmp = l_p * r_p;
					prflinfo.lamda_LL += tmp;
					prflinfo.p_LL.push_back(tmp);
					//prflinfo.lamda_pp = gParam.lenpoly_pe * gParam.lenpoly_pe;
				}

				if (!r_is_ins && rtype.at(0) != 'D' && r_s) { // ------ PS ------
					//prflinfo.lamda_ps +=
					//	std::min(l_p, gParam.lenpoly_pe) * std::min(r_s, gParam.pe);
					double tmp = l_p * r_s;
					prflinfo.lamda_LS += tmp;
					prflinfo.p_LS.push_back(tmp);
					//prflinfo.lamda_ps += gParam.lenpoly_pe * gParam.pe;
				}
			}

			if (!l_is_ins && ltype.at(0) != 'D' && l_s) {
				if (r_p) { // ------ SP ------
					//prflinfo.lamda_sp +=
					//	std::min(l_s, gParam.pe) * std::min(r_p, gParam.lenpoly_pe);
					double tmp = l_s * r_p;
					prflinfo.lamda_SL += tmp;
					prflinfo.p_SL.push_back(tmp);
					//prflinfo.lamda_sp += gParam.pe * gParam.lenpoly_pe;
				}

				if (!r_is_ins && rtype.at(0) != 'D' && r_s) { //------ SS ------
					//prflinfo.lamda_ss +=
					//	std::min(l_s, gParam.pe) * std::min(r_s, gParam.pe);
					double tmp = l_s * r_s;
					prflinfo.lamda_SS += tmp;
					prflinfo.p_SS.push_back(tmp);
					//prflinfo.lamda_ss += gParam.pe * gParam.pe;
				}
			}

			++ l;
			++ r;
		} else {
			if (lhs.entries[l].rID > rhs.entries[r].rID) ++ r;
			else ++ l;
		}
	} // while (l < l_sz && r < r_sz) {
} // phase_profiling

/* @brief	Add p1 to p0's neighborhood
 */
void add_to_phaselist (std::map<int, iset_t>& list, int p0, int p1) {
	std::map<int, iset_t>::iterator it_l = list.find(p0);
	if (it_l != list.end()) it_l->second.insert(p1);
	else {
		iset_t tmp;
		tmp.insert(p1);
		list[p0] = tmp;
	}
} //add_to_phaselist


/* @brief	Given two ref positions, determine if new phase variants exist
 * 			if so, store the results in [l_vars] and [r_vars], otherwise,
 * 			determine if the phase columns contain any phasing variants
 * 			that needs to be checked further, if so, return true;
 *
 */
bool call_phase_variant (strset_t& l_vars, strset_t& r_vars,	int l_pos,
	int r_pos, int threshold, const strimap_t& profile, const var_t& gvars) {

	bool is_to_further_check = false;

	strimap_t::const_iterator it_p;
	for (it_p = profile.begin(); it_p != profile.end(); ++ it_p){
		if (it_p->second >= threshold) {
			strvec_t splits;
			jaz::split(',', it_p->first, std::back_inserter(splits));
			if (!is_var_exist (l_pos, splits[0], gvars) &&
					!is_var_exist (r_pos, splits[1], gvars)) {

				l_vars.insert(splits[0]);
				r_vars.insert(splits[1]);

				{ //xxxxxxxxxxxx
//					std::cout << "pos: " << l_pos/2 << ", " << r_pos/2 << ":\t";
//					std::cout << splits[0] << ", " << splits[1] << " -- "
//					<< it_p->second << "\n";
				}
			}
		} else if (it_p->second > 1) {
			strvec_t splits;
			jaz::split(',', it_p->first, std::back_inserter(splits));
			if (!is_var_exist (l_pos, splits[0], gvars) &&
					!is_var_exist (r_pos, splits[1], gvars)) {
				is_to_further_check = true;
				{
		//			std::cout << "skip pos: " << l_pos/2 << ", " << r_pos/2 << ":\t";
		//			std::cout << splits[0] << ", " << splits[1] << " -- "
		//			<< it_p->second << "\n";
				}
			}
		}
	} // for (it_p
	return is_to_further_check;
} // call_phase_variant

/* @brief	Check if p1 is in p0's neighborhood
 */
bool is_in_phase (const std::map<int, iset_t>& list, int p0, int p1) {
	std::map<int, iset_t>::const_iterator it_l = list.find(p0);
	if (it_l != list.end()) {
		if (it_l->second.count(p1)) return true;
		else return false;
	} else return false;
} // is_in_phase

/* @brief	return true if a specific var exists in [gvars] already
 */
bool is_var_exist (int pos, const std::string& var, const var_t& gvars){
	var_t::const_iterator it_g = gvars.find(pos);
	if (it_g != gvars.end()) {
		if (it_g->second.count(var)) return true;
		else return false;
	} else return false;
}

//========================================================================
// Project     : VariantCaller
// Name        : em.cpp
// Author      : Xiao Yang
// Created on  : Apr 21, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================

#include "em.h"

/* @brief
 *
 */
void EM (GlobalParam& gParam, const Parameter& myPara){

	std::cout << "\n\tCalibrate errors...\n\n";

	// obtain jeb, and calculate bucket size info
	std::vector<jeb_t> jeb;
	read_jeb (gParam, jeb);
	//debug_print_jeb (jeb);

	int file_idx = 0;
	Ref::iterator it_ref = gParam.reference.begin();

	for (; it_ref != gParam.reference.end(); ++ it_ref) {

		// record refPos with variants called so that the whole column
		// will be removed from jeb calibration
		std::map<int, bpair_t> called; // ref_pos --> (lp, snp)
		var_t gvars; //std::map<int, strset_t> global variants
		pileup_list_t pileuplist;// pileup candidate positions to be checked
		phase_nb_t phase_nb;  // phase position and its neighbors

		int iter = 0;
		while (iter < 10) {

			bool is_var_found = false;

			std::cout << "\n\t\tIteration " << iter << "\n\n";

			int num_file = gParam.alnfiles[file_idx].size();

 			for (int f = 0; f < num_file; ++ f) {

				std::cout << "\t\tProcess " << gParam.alnfiles[file_idx][f] << "\n";

				std::vector<col_t> cols;
				ipair_t region;
				read_aln_file (region, cols, ivec_t(), gParam.alnfiles[file_idx][f]);
 				// ----------------------------------------------------------
				/* process column : calculate errbucket index, pe, type
				 calculate boundary to make sure pos0 <= regions.second */
				int num_col = cols.size();
				int maxPileupIdx = -1;

				//#pragma omp parallel for shared (cols)
				// identify complete aln column boundary
				for (int i = 0; i < num_col; ++ i) {
					if (cols[i].ref_pos > region.second) {
						maxPileupIdx = i;
						break;
					}
				}
				if (maxPileupIdx == -1) maxPileupIdx = num_col;

				// ----------------------------------------------------------
				// actual em
				if (infer_variant (gvars, pileuplist, phase_nb, cols,
						maxPileupIdx, jeb, myPara, gParam, iter)) {
					is_var_found = true;
				}


			} // for (int f = 0;


			{ // print out new_vars
			//	debug_print_var (ivars);
			//	std::cout << "Total Var pos: " << ivars.size() << "\n";
			}

			if (! is_var_found) break;

			std::cout << "\n\t\tTo consider in the new iteration:\n";
			std::cout << "\n\t\t\t" << pileuplist.S.size() << " snp";
			std::cout << "\n\t\t\t" << pileuplist.L.size() << " lv";
			std::cout << "\n\t\t\t" << phase_nb.SS.size() << " snp_snp";
			std::cout << "\n\t\t\t" << phase_nb.LS.size() << " lv_snp";
			std::cout << "\n\t\t\t" << phase_nb.SL.size() << " snp_lv";
			std::cout << "\n\t\t\t" << phase_nb.LL.size() << " lv_lv\n";
			std::cout << "\n\t\tUpdate Variants\n";

			update_variants (jeb, called, gParam.bkinfo, gvars,
					gParam.alnfiles[file_idx], myPara);

			gParam.bkinfo.calculate_prior();
			gParam.bkinfo.print();

			{
			//	debug_print_var (gvars);
			//	std::cout << "\tTotal Var pos: " << gvars.size() << "\n\n";
			//	abording ("requested exit");
			}

			++ iter;
 		//	abording ("requested exit");

		} // while (iter <

		finalize_variant (gvars, it_ref->first, gParam.alnfiles[file_idx],
				myPara);
		++ file_idx;

	} // for (; it_ref;
} // EM

/* @brief
 */
bool infer_variant (var_t& gvars, pileup_list_t& pileuplist,
		phase_nb_t& phase_nb, const std::vector<col_t>& cols,
		int maxPileupIdx, const std::vector<jeb_t>& jeb,
		const Parameter& myPara, const GlobalParam& gParam, int iter) {

	bool is_new_var = false;
	int num_cols = cols.size();

	#pragma omp parallel // --- parallel
	{
		var_t nuvars;
		pileup_list_t nulist, rmlist;
		phase_nb_t nu_nb, rm_nb;

		#pragma omp for schedule (dynamic, 20) // --- parallel
		for (int i = 0; i < maxPileupIdx; ++ i) {

			int p0 = cols[i].ref_pos;

			{
				//if ((p0/2) % 50 == 0) std::cout << p0/2 << "\n";
			}

			if (iter == 0) {
				pileup (nuvars, nulist, rmlist, pileuplist, gvars,
						cols[i], jeb, gParam, iter);
			}

			// phasing from 2nd iter
			if (iter && (myPara.errModel == 1)) {
				// neither position contains vars during pileup
				if (pileuplist.L.count(p0) || pileuplist.S.count(p0)) {

					for (int j = i + 1; j < num_cols; ++ j) {

						int p1 = cols[j].ref_pos;

						if (pileuplist.L.count(p1) || pileuplist.S.count(p1)) {

							phase (nuvars, nu_nb, rm_nb, phase_nb, pileuplist,
								gvars, cols[i], cols[j], jeb, gParam, iter);
						}
					} // for (int j =

				} // if (pileuplist

			} // if (iter)

		} // for (int i = 0;


		#pragma omp critical // --- parallel
		{ // for [nuvars] [nu_list] [nu_nb]
			if (nuvars.size()) is_new_var = true;
			var_t::iterator it = nuvars.begin();
			for (; it != nuvars.end(); ++ it) {
				add_var_entry (gvars, it->first, it->second);
			}

			pileuplist.L.insert(nulist.L.begin(), nulist.L.end());
			pileuplist.S.insert(nulist.S.begin(), nulist.S.end());

			for (iset_t::iterator it = rmlist.L.begin();
				it != rmlist.L.end(); ++ it) pileuplist.L.erase(*it);

			for (iset_t::iterator it = rmlist.S.begin();
				it != rmlist.S.end(); ++ it) pileuplist.S.erase(*it);

			merge_phase_nb (phase_nb.LL, nu_nb.LL);
			merge_phase_nb (phase_nb.LS, nu_nb.LS);
			merge_phase_nb (phase_nb.SL, nu_nb.SL);
			merge_phase_nb (phase_nb.SS, nu_nb.SS);

			subtract_phase_nb (phase_nb.LL, rm_nb.LL);
			subtract_phase_nb (phase_nb.LS, rm_nb.LS);
			subtract_phase_nb (phase_nb.SL, rm_nb.SL);
			subtract_phase_nb (phase_nb.SS, rm_nb.SS);

		}
	}

	return is_new_var;
} // calibrate_pe

/*	@brief	merge [from] to [to]
 */
void merge_phase_nb (std::map<int, iset_t>& to,
		const std::map<int, iset_t>& from) {

	std::map<int, iset_t>::const_iterator it_from;
	std::map<int, iset_t>::iterator it_to;
	for (it_from = from.begin(); it_from != from.end(); ++ it_from) {
		it_to = to.find(it_from->first);
		if (it_to != to.end()) {
			it_to->second.insert(it_from->second.begin(),
					it_from->second.end());
		} else to.insert(*it_from);
	}
} // merge_phase_nb

/* @brief	Substract rhs from lhs
 *
 */
void subtract_phase_nb (std::map<int, iset_t>& lhs,
		const std::map<int, iset_t>& rhs) {
	std::map<int, iset_t>::const_iterator it_r;
	std::map<int, iset_t>::iterator it_l;
	for (it_r = rhs.begin(); it_r != rhs.end(); ++ it_r) {
		it_l = lhs.find(it_r->first);
		if (it_l != lhs.end()) {
			iset_t::const_iterator it_s = it_r->second.begin();
			for (; it_s != it_r->second.end(); ++ it_s){
				if (it_l->second.size()) it_l->second.erase (*it_s);
				else lhs.erase(it_l);
			}
		} else {
			abording ("subtract_phase_nb: nb not found !");
		}
	}
} // merge_phase_nb

/* @brief	Check if a pos contain any non-consensus variant
 */
bool is_contain_vars (const strset_t& strvars) {
	int lvcnt = 0, snpcnt = 0;
	strset_t::const_iterator it_s = strvars.begin();
	for (; it_s != strvars.end(); ++ it_s) {
		char v0 = it_s->at(0);
		if (v0 == 'I' || v0 == 'D' || v0 == 'd' || v0 == 'i') ++ lvcnt;
		else ++ snpcnt;
	} // for
	if (lvcnt > 1 || snpcnt > 1) return true;
	else return false;
} // is_contain_vars

/* @brief
 * 	1. identify which columns contain non-consensus vars (.size() > 1)
 * 		that were not previously updated for jeb
 * 	2. Fetch those columns
 * 	3. Update [jeb] and [called]
 *			return true if new variants are found.
 */
void update_variants (std::vector<jeb_t>& jeb, std::map<int, bpair_t>& called,
		Bucket& bkinfo, 	const var_t& gvars, const strvec_t& alnfiles,
		const Parameter& myPara) {

	ivec_t ref_col_to_fetch;
	var_t::const_iterator it_g = gvars.begin();
	for (; it_g != gvars.end(); ++ it_g) {
		if (is_contain_vars (it_g->second)) { // contain non-consensus variant
			int ref_pos = it_g->first;
			std::map<int, bpair_t>::iterator it_c = called.find(ref_pos);
			if (it_c != called.end ()) {
				if (it_c->second.first == false || // lp or snp for non-ins column
				  ((ref_pos % 2 == 0) && it_c->second.second == false)) {
					 ref_col_to_fetch.push_back(ref_pos);
				}
			} else ref_col_to_fetch.push_back(ref_pos);
		}
	}

	std::cout << "\n\t\t\tFetching " << ref_col_to_fetch.size() << " cols\n";
	// --------- fetch columns and store in [cols] --------------
	std::vector<col_t> cols;
	for (int i = 0; i < (int) alnfiles.size(); ++ i) {
		ipair_t region;
		read_aln_file (region, cols, ref_col_to_fetch, alnfiles[i]);
	}

	// --------- check gvars and update jeb when it hasn't
	// --------- been considered (based on [called]) previously -----
	clean_jeb (jeb, called, bkinfo, gvars, cols);

} // update_variants

/* @brief	Go over all columns in [cols], and check if current column
 * 			was used to update jeb previously, if not, update [jeb] and
 * 			[called]
 */
void clean_jeb (std::vector<jeb_t>& jeb, std::map<int, bpair_t>& called,
	Bucket& bkinfo, const var_t& gvars,	const std::vector<col_t>& cols) {

	int num_update_lv = 0, num_update_snp = 0;
	for (int i = 0; i < (int) cols.size(); ++ i) {

		std::string cons_L, cons_S;
		get_consensus (cons_L, cons_S, cols[i]);

		int ref_pos = cols[i].ref_pos;
		// find this ref_pos in [gvars]
		var_t::const_iterator it_g = gvars.find(ref_pos);
		if (it_g == gvars.end()) { // SC
			std::cout << "\nref_pos = " << ref_pos << "\n";
			abording ("clean_jeb: can't find ref_pos in gvars");
		}

		// for each variant in ref_pos
		strset_t::const_iterator it_str = it_g->second.begin();

		for (; it_str != it_g->second.end(); ++ it_str) {

			if (it_str->compare(cons_L) != 0 && it_str->compare(cons_S) != 0){
				ipair_t	update_cnt (0, 0);
				char v0 = it_str->at(0);
				std::map<int, bpair_t>::iterator it_c = called.find(ref_pos);

				ipair_t cnt;
				if (v0 == 'I' || v0 == 'D' || v0 == 'd' || v0 == 'i') {
					if (it_c == called.end()){
						cnt = decre_jeb (jeb, cons_L, cols[i].entries, true);
						called[ref_pos] = bpair_t (true, false);
						++ num_update_lv;
					} else if (!it_c->second.first) {
						cnt = decre_jeb (jeb, cons_L, cols[i].entries, true);
						it_c->second.first = true;
						++ num_update_lv;
					}
					bkinfo.decre_err_L(cnt.first);
					bkinfo.decre_sum_L(cnt.second);
				} else { // snp variant
					if (it_c == called.end()){
						cnt = decre_jeb (jeb, cons_S, cols[i].entries, false);
						called[ref_pos] = bpair_t (false, true);
						++ num_update_snp;
					} else if (!it_c->second.second) {
						cnt = decre_jeb (jeb, cons_S, cols[i].entries, false);
						it_c->second.second = true;
						++ num_update_snp;
					}
					bkinfo.decre_err_S(cnt.first);
					bkinfo.decre_sum_S(cnt.second);
				}
			}// if
		} // for (;
	} // for (int i = 0; i < (int) cols.size(); ++ i)

	std::cout << "\n\t\t\t" << num_update_snp << " snp updates for jeb\n";
	std::cout << "\n\t\t\t" << num_update_lv << " lv updates for jeb\n";

} // clean_jeb

/* @brief	Return # numerator & demoninator
 */
ipair_t decre_jeb (std::vector<jeb_t>& jeb, const std::string& cons,
		const std::vector<entry_t>& entries, bool is_lp) {

	ipair_t update (0, 0); // counting # updated numerator & denominator
	update.first = update.second = 0;

	if (is_lp) {
		if (cons.at(0) == 'd' || cons.at (0) == 'i') {
			for (int i = 0; i < (int) entries.size(); ++ i) {
				char c0 = entries[i].cons_type.at(0);
				int idx = entries[i].eb_index;

				if (c0 == 'I' || c0 == 'D') {
					if (jeb[idx].lp_err <= 0) { // SC
						std::cout << "jeb index = " << idx << "\n";
						abording ("decre_jeb: d, i case err_cnt <= 0");
					}
					-- jeb[idx].lp_err;
					++ update.first;
				}
				-- jeb[idx].lp_sum;
				++ update.second;
			}
		} else { // cons is I or D
			for (int i = 0; i < (int) entries.size(); ++ i) {
				int idx = entries[i].eb_index;

				if (cons.compare(entries[i].cons_type) != 0) {
					if (jeb[idx].lp_err <= 0) { // SC
						std::cout << "jeb index = " << idx << "\n";
						abording ("decre_jeb: I, D case err_cnt <= 0");
					}
					-- jeb[idx].lp_err;
					++ update.first;
				}
				-- jeb[idx].lp_sum;
				++ update.second;
			}
		}
	} else {
		for (int i = 0; i < (int) entries.size(); ++ i) {
			char c0 = entries[i].cons_type.at(0);
			int idx = entries[i].eb_index;
			if (xny::is_nt(c0)) {
				if (cons.compare(entries[i].cons_type) != 0) {
					if (jeb[idx].snp_err <= 0) { // SC
						std::cout << "jeb index = " << idx << "\n";
						abording ("decre_jeb: subst case cnt <= 0");
					}
					-- jeb[idx].snp_err;
					++ update.first;
				}
				-- jeb[idx].snp_sum;
				++ update.second;
			}
		}
	}
	return update;
}// decre_jeb


/* @brief	Given an aln column, determine length polymorphic consensus,
 * 			and substitution consensus.
 */
void get_consensus (std::string& cons_p, std::string& cons_s, const col_t& col) {

	strimap_t profile, lprofile;

	// obtain both types of profiles
	bool is_ins = (col.ref_pos % 2 == 1) ? true : false;
	for (int i = 0; i < (int) col.entries.size(); ++ i) {
		if (col.entries[i].cons_type.at(0) == 'D'
				|| col.entries[i].cons_type.at(0) == 'I') {
			add_to_profile (lprofile, col.entries[i].cons_type);
		} else {
			if (is_ins)	add_to_profile (lprofile, "d");
			else {
				add_to_profile (lprofile, "i");
				add_to_profile (profile, col.entries[i].cons_type);
			}
		}
	}

	ivec_t cnts;
	analyze_profile (cons_p, cnts, lprofile);
	analyze_profile (cons_s, cnts, profile);

} // get_consensus

/* @brief	Finalize variant: get each var and profile information for
 * the corresponding aln column; strandbias test; then output
 */
void finalize_variant (const var_t& gvars, const std::string& refName,
		const strvec_t& alnfiles, const Parameter& myPara){


	// ref pos to be retrieved from file
	ivec_t ref_pos_to_fetch;
	std::vector<strset_t> list_vars;
	var_t::const_iterator it_g = gvars.begin();
	for (; it_g != gvars.end(); ++ it_g) {
		// contain non-consensus variant
		if (is_contain_vars (it_g->second)) {
		//if (it_g->second.size() > 1) {
			ref_pos_to_fetch.push_back(it_g->first);
			list_vars.push_back(it_g->second);
		}
	}

	// --------- fetch targeted columns in [cols] --------------
	std::vector<col_t> cols;
	for (int i = 0; i < (int) alnfiles.size(); ++ i) {
		ipair_t region;
		read_aln_file (region, cols, ref_pos_to_fetch, alnfiles[i]);
	}

	// --------- analyze vars --------------------
	std::vector<varinfo_t> list_varInfo;
	if (list_vars.size()) { // fixed bug, if no var exists do not run analyze_vars
		analyze_vars (list_varInfo, myPara.alpha, cols, list_vars);
	}

	// --------- counting strand information for each var and cons ----
	if (list_varInfo.size()) { // fixed bug
		output_vars(list_varInfo, myPara.oDirNm, refName);
	}
} // finalize_variant

void output_vars (const std::vector<varinfo_t>& list_varInfo,
		const std::string& dir, const std::string& refName){

	std::string path_nofdr = dir, path_nofdr_detail = dir, path_fdr = dir,
			path_fdr_detail = dir, path_raw = dir, path_raw_detail = dir;

//	path_nofdr += "/" + refName + ".nofdr.var.txt";
	path_nofdr_detail += "/" + refName + ".nofdr.var.txt";
//	path_fdr += "/" + refName + ".fdr.var.txt";
	path_fdr_detail += "/" + refName + ".fdr.var.txt";
//	path_raw += 	"/" + refName + ".var.raw.txt";
	path_raw_detail += "/" + refName + ".var.raw.txt";

	std::stringstream ss;
	std::ofstream oFdr, oFdrd, oRaw, oRawd, oNoFdr, oNoFdrd;

//	create_output_file (oFdr, path_fdr);
	create_output_file (oFdrd, path_fdr_detail);
//	create_output_file (oNoFdr, path_nofdr);
	create_output_file (oNoFdrd, path_nofdr_detail);
//	create_output_file (oRaw, path_raw);
	create_output_file (oRawd, path_raw_detail);

	int sz = list_varInfo.size();
//	oFdr << "# Ref_Pos\tVar\tCons\n";
//	oNoFdr << "# Ref_Pos\tVar\tCons\n";
//	oRaw << "# Ref_Pos\tVar\tCons\n";
	oFdrd << "# Ref_Pos\tVar\tCons\tStrd_bias_pval\tType\tVar_perc\tSNP_or_LP_Profile\n";
	oNoFdrd << "# Ref_Pos\tVar\tCons\tStrd_bias_pval\tType\tVar_perc\tSNP_or_LP_Profile\n";
	oRawd << "# Ref_Pos\tVar\tCons\tStrd_bias_pval\tType\tVar_perc\tSNP_or_LP_Profile\n";
//	oFdr << "# ------------------------------------------------------------\n";
//	oNoFdr << "# ------------------------------------------------------------\n";
//	oRaw << "# ------------------------------------------------------------\n";
	oFdrd << "# ------------------------------------------------------------\n";
	oNoFdrd << "# ------------------------------------------------------------\n";
	oRawd << "# ------------------------------------------------------------\n";

	int cnt_S_fdr = 0, cnt_S_nofdr = 0, cnt_S_raw = 0,
		cnt_L_fdr = 0, cnt_L_nofdr = 0, cnt_L_raw = 0;
	std::map<std::string, ipair_t>::const_iterator it_p;
	for (int i = 0; i < sz; ++ i) {

		// -------------------- identify consensus --------------
		std::string cons;
		int cons_cnt = -1;
		it_p = list_varInfo[i].profile_strd.begin();
		for (; it_p != list_varInfo[i].profile_strd.end(); ++ it_p) {
			int tmp = it_p->second.first +  it_p->second.second;
			if (tmp > cons_cnt) {
				cons = it_p->first;
				cons_cnt = tmp;
			} else if (tmp == cons_cnt && it_p->first.compare(cons) < 0) {
				cons = it_p->first;
			}
		}

/*		oRaw << std::left << std::setw(9) << list_varInfo[i].ref_pos/2
				<< std::setw(8) << list_varInfo[i].var << std::setw(8)
				<< cons;
		oRawd << std::left << std::setw(9) << list_varInfo[i].ref_pos/2
				<< std::setw(8) << list_varInfo[i].var << std::setw(8)
				<< cons	<< std::setprecision(4) << std::fixed
				<< std::setw(10) << list_varInfo[i].p_value;
*/
		int ref_pos = list_varInfo[i].ref_pos/2 + 1;
		oRawd << ref_pos << "\t" << list_varInfo[i].var	<< "\t"	<< cons
				<< "\t" << std::setprecision(4) << list_varInfo[i].p_value;

		// ---------------------------
		if (list_varInfo[i].is_pass_strd_fdr_test) {
			/*
			oFdr << std::left << std::setw(9) << list_varInfo[i].ref_pos/2
				<< std::setw(8) << list_varInfo[i].var << std::setw(8)
				<< cons;
			oFdrd << std::left << std::setw(9) << list_varInfo[i].ref_pos/2
				<< std::setw(8) << list_varInfo[i].var << std::setw(8)
				<< cons	<< std::setprecision(4) << std::fixed
				<< std::setw(10) << list_varInfo[i].p_value;
			*/
			oFdrd << ref_pos << "\t" << list_varInfo[i].var << "\t" << cons
				<< "\t" << std::setprecision(4) << list_varInfo[i].p_value;

		}
		// ---------------------------
		if (list_varInfo[i].is_pass_strd_test) {
			/*
			oNoFdr << std::left << std::setw(9) << list_varInfo[i].ref_pos/2
				<< std::setw(8) << list_varInfo[i].var << std::setw(8)
				<< cons;
			oNoFdrd << std::left << std::setw(9) << list_varInfo[i].ref_pos/2
				<< std::setw(8) << list_varInfo[i].var << std::setw(8)
				<< cons	<< std::setprecision(4) << std::fixed
				<< std::setw(10) << list_varInfo[i].p_value;
			*/
			oNoFdrd << ref_pos << "\t" << list_varInfo[i].var << "\t" << cons
					<< "\t" << std::setprecision(4) << list_varInfo[i].p_value;
		}

		if (list_varInfo[i].is_snp) {
			++ cnt_S_raw;
//			oRawd << std::setw(7) << "snp";
			oRawd << "\t" << "snp";

			if (list_varInfo[i].is_pass_strd_test) {
				++ cnt_S_nofdr;
//				oNoFdrd << std::setw(7) << "snp";
				oNoFdrd << "\t" << "snp";
			}
			if (list_varInfo[i].is_pass_strd_fdr_test) {
				++ cnt_S_fdr;
//				oFdrd << std::setw(7) << "snp";
				oFdrd << "\t" << "snp";
			}
		} else {
			++ cnt_L_raw;
//			oRawd << std::setw(7) << "lp";
			oRawd << "\t" << "lp";
			if (list_varInfo[i].is_pass_strd_test) {
				++ cnt_L_nofdr;
//				oNoFdrd << std::setw(7) << "lp";
				oNoFdrd << "\t" << "lp";
			}
			if (list_varInfo[i].is_pass_strd_fdr_test) {
				++ cnt_L_fdr;
//				oFdrd << std::setw(7) << "lp";
				oFdrd << "\t" << "lp";
			}
		}

		// calculate frequency first
		int sum_freq = 0, var_freq = 0;
		it_p = list_varInfo[i].profile_strd.begin();
		for (; it_p != list_varInfo[i].profile_strd.end(); ++ it_p) {
			if (it_p->first.compare(list_varInfo[i].var) == 0) {
				var_freq += it_p->second.first + it_p->second.second;
			}
			sum_freq += it_p->second.first + it_p->second.second;
		}

		// output frequency
		double frq = 0.0;
		if (sum_freq) frq = var_freq * 100.0/sum_freq;
		oRawd << "\t" << std::setprecision(4) << frq;
		if (list_varInfo[i].is_pass_strd_test) {
			oNoFdrd << "\t" << std::setprecision(4) << frq;
		}
		if (list_varInfo[i].is_pass_strd_fdr_test) {
			oFdrd << "\t" << std::setprecision(4) << frq;
		}

		it_p = list_varInfo[i].profile_strd.begin();
		for (; it_p != list_varInfo[i].profile_strd.end(); ++ it_p) {
//			ss.str("");
//			ss << it_p->first << ":" << it_p->second.first
//					<< ":" << it_p->second.second;
//			oRawd << std::setw(20) << ss.str();
			oRawd << "\t" << it_p->first << ":" << it_p->second.first
										<< ":" << it_p->second.second;
			if (list_varInfo[i].is_pass_strd_test) {
//				oNoFdrd << std::setw(20) << ss.str();
				oNoFdrd << "\t" << it_p->first << ":" << it_p->second.first
											<< ":" << it_p->second.second;
			}
			if (list_varInfo[i].is_pass_strd_fdr_test) {
//				oFdrd << std::setw(20) << ss.str();
				oFdrd << "\t" << it_p->first << ":" << it_p->second.first
											<< ":" << it_p->second.second;
			}
		}

//		oRaw << "\n";
		oRawd << "\n";

		if (list_varInfo[i].is_pass_strd_test) {
//			oNoFdr << "\n";
			oNoFdrd << "\n";
		}
		if (list_varInfo[i].is_pass_strd_fdr_test) {
//			oFdr << "\n";
			oFdrd << "\n";
		}
	}

//	oFdr << "# ------------------------------------------------------------\n";
//	oNoFdr << "# ------------------------------------------------------------\n";
//	oRaw << "# ------------------------------------------------------------\n";
	oFdrd << "# ------------------------------------------------------------\n";
	oNoFdrd << "# ------------------------------------------------------------\n";
	oRawd << "# ------------------------------------------------------------\n";

//	oFdr << "# Summary: SNPV: " << cnt_S_fdr << ";\t LPV: " << cnt_L_fdr << "\n";
//	oNoFdr << "# Summary: SNPV: " << cnt_S_nofdr << ";\t LPV: " << cnt_L_nofdr << "\n";
//	oRaw << "# Summary: SNPV: " << cnt_S_raw << ";\t LPV: " << cnt_L_raw << "\n";
	oFdrd << "# Summary: SNPV: " << cnt_S_fdr << ";\t LPV: " << cnt_L_fdr << "\n";
	oNoFdrd << "# Summary: SNPV: " << cnt_S_nofdr << ";\t LPV: " << cnt_L_nofdr << "\n";
	oRawd << "# Summary: SNPV: " << cnt_S_raw << ";\t LPV: " << cnt_L_raw << "\n";
//	oFdr.close();
	oFdrd.close();
//	oNoFdr.close();
	oNoFdrd.close();
//	oRaw.close();
	oRawd.close();
} // output_vars


void analyze_vars (std::vector<varinfo_t>& list_varInfo, const double& alpha,
	const std::vector<col_t>& cols, const std::vector<strset_t>& list_vars){

	int sz = cols.size();
	for (int i = 0; i < sz; ++ i) {

		std::string cons_p, cons_s;
		get_consensus (cons_p, cons_s, cols[i]);

		std::map<std::string, ipair_t> sprofile, lprofile;

		//----- obtain snp and lp strand profiles -------------------------
		bool is_ins = (cols[i].ref_pos % 2 == 1) ? true : false;
		for (int j = 0; j < (int) cols[i].entries.size(); ++ j) {
			if (cols[i].entries[j].cons_type.at(0) == 'D'
					|| cols[i].entries[j].cons_type.at(0) == 'I') {
				add_to_profile_strand (lprofile, cols[i].entries[j].cons_type,
						cols[i].entries[j].is_rv);
			} else {
				if (is_ins)	add_to_profile_strand (lprofile, "d",
						cols[i].entries[j].is_rv);
				else {
					add_to_profile_strand (lprofile, "i", cols[i].entries[j].is_rv);
					add_to_profile_strand (sprofile, cols[i].entries[j].cons_type,
							cols[i].entries[j].is_rv);
				}
			}
		}

		// ------ update [list_varInfo] --------------------
		//if (list_vars.size() == 0)
		strset_t::const_iterator it_str = list_vars[i].begin();
		for (;	it_str != list_vars[i].end(); ++ it_str) {
			if (it_str->compare(cons_p) != 0 &&
					it_str->compare(cons_s) != 0){
				if (it_str->at(0) == 'D' || it_str->at(0) == 'I' ||
					it_str->at(0) == 'i' || it_str->at(0) == 'd'){

					list_varInfo.push_back(varinfo_t(false, cols[i].ref_pos,
							*it_str, 0.0, lprofile));
				} else {
					list_varInfo.push_back(varinfo_t(true, cols[i].ref_pos,
							*it_str, 0.0, sprofile));
				}
			}
		}// for
	} // for (int i

	strdbias_pval (list_varInfo);
	benjamini_hchberg_fdr (list_varInfo, alpha);
	// pval threshold applies to every single entry
	for (int i = 0; i < (int) list_varInfo.size(); ++ i) {
		if (list_varInfo[i].p_value <= alpha) {
			list_varInfo[i].is_pass_strd_test = false;
		}
	}
} //analyze_vars

/* @brief	FDR control of family err (better than bonferroni)
 */
void benjamini_hchberg_fdr (std::vector<varinfo_t>& list_varInfo,
		const double& alpha) {
	std::sort (list_varInfo.begin(), list_varInfo.end(), cmp_varinfo_pval());
	int m = list_varInfo.size();

	for (int i = 0; i < m; ++ i) {
		if (list_varInfo[i].p_value <= (i + 1.0) * alpha/m) {
			list_varInfo[i].is_pass_strd_fdr_test = false;
		}
		//std::cout << std::setw(5) << list_varInfo[i].p_value << " ,";
	}
	std::cout << "\n";

	std::sort (list_varInfo.begin(), list_varInfo.end(), cmp_varinfo_refpos());
} // benjamini_hchberg_fdr

/* @brief	Apply Chi-sqr or fisher's exact test, calculate p-value
 *			for each entry of list_varInfo
 *				var  not_var
 *			fwd	x00		x01
 *			rv  x10 		x11
 */
void strdbias_pval (std::vector<varinfo_t>& list_varInfo) {
	int n = list_varInfo.size();
	std::map<std::string, ipair_t>::iterator it_p;
	for (int i = 0; i < n; ++ i) {
		int x00 = 0, x01 = 0, x10 = 0, x11 = 0;
		it_p = list_varInfo[i].profile_strd.begin();
		for (; it_p != list_varInfo[i].profile_strd.end(); ++ it_p) {
			if (it_p->first.compare(list_varInfo[i].var) == 0) {
				x00 = it_p->second.first;
				x10 = it_p->second.second;
			} else {
				x01 += it_p->second.first;
				x11 += it_p->second.second;
			}
		}
		/* conservative approach to test a) make sure the larger number
		 * is on the same side, and b) when an entry is 0, let p-val = 0*/
		if ((x00 < x10 && x01 > x11) || (x00 > x10 && x01 < x11)) {
			//swap x01 and x11
			std::swap(x01, x11);
		}
		if (!x00 || !x01 || !x10 || !x11) list_varInfo[i].p_value = 0;
		else {
			/* create table for stat test */
			iivec_t mytable;
			ivec_t row (2);
			row[0] = x00; row[1] = x01;
			mytable.push_back(row);
			row[0] = x10; row[1] = x11;
			mytable.push_back(row);
			//std::cout << x00 << ",\t" << x01 << ",\t" << x10 << ",\t" << x11 << "\n";
			double p_value = chi_sqr_pvalue (mytable);
			if (p_value == -1.0) {
				p_value = tbt_fisher_exact_twotail_pval (x00, x01, x10, x11);
			}
			list_varInfo[i].p_value = p_value;
		}
	}
} // strdbias_test

/* @brief	Write out the variants for a specified ref col
 */
void write_var (int& cnt_L, int& cnt_S, std::ofstream& oHandle,
		const col_t& col, var_t::const_iterator it_v) {

	std::string cons_p, cons_s;
	get_consensus (cons_p, cons_s, col);

	std::map<std::string, ipair_t> profile, lprofile;
	std::map<std::string, ipair_t>::iterator it_p;

	//----- obtain both types of profiles, and strand count  ------------
	bool is_ins = (col.ref_pos % 2 == 1) ? true : false;
	for (int i = 0; i < (int) col.entries.size(); ++ i) {

		if (col.entries[i].cons_type.at(0) == 'D'
				|| col.entries[i].cons_type.at(0) == 'I') {
			add_to_profile_strand (lprofile, col.entries[i].cons_type,
					col.entries[i].is_rv);
		} else {
			if (is_ins)	add_to_profile_strand (lprofile, "d",
					col.entries[i].is_rv);
			else {
				add_to_profile_strand (lprofile, "i", col.entries[i].is_rv);
				add_to_profile_strand (profile, col.entries[i].cons_type,
						col.entries[i].is_rv);
			}
		}
	}

	// ------ write out variants and profiles --------------------
	bool is_poly_var = false, is_subst_var = false;
	strset_t::const_iterator it_str = it_v->second.begin();
	for (;	it_str != it_v->second.end(); ++ it_str) {
		if (it_str->compare(cons_p) != 0 && it_str->compare(cons_s) != 0){
			if (!is_poly_var && !is_subst_var) {
				oHandle << it_v->first << " (" << it_v->first/2 << "):\t\t" ;
			}

			if (it_str->at(0) == 'D' || it_str->at(0) == 'I' ||
				it_str->at(0) == 'i' || it_str->at(0) == 'd'){
				is_poly_var = true;
			} else is_subst_var = true;
			oHandle << *it_str << "\t";
		}
	}

	if (is_subst_var) {
		++ cnt_S;
		oHandle << "\n                     snp\t\t";
		for (it_p = profile.begin(); it_p != profile.end(); ++ it_p) {
			oHandle  << it_p->first << ":" << it_p->second.first
					<< "," << it_p->second.second << "\t\t";
		}
		oHandle << "\n";
	}
	if (is_poly_var) {
		++ cnt_L;
		oHandle << "\n                     len\t\t";
		for (it_p = lprofile.begin(); it_p != lprofile.end(); ++ it_p) {
			oHandle  << it_p->first << ":" << it_p->second.first
					<< "," << it_p->second.second << "\t\t";
		}
		oHandle << "\n";
	}

} // write_var

//------------------------ debug printing facilities ---------------------
void debug_print_var (const var_t& vars) {
	for (var_t::const_iterator it = vars.begin(); it != vars.end();
						++ it) {
		std::cout << it->first << " (" << it->first/2 << " ):\t" ;
		for (strset_t::const_iterator it_str = it->second.begin();
			it_str != it->second.end(); ++ it_str) {
			std::cout << *it_str << "\t";
		}
		std::cout << "\n";
	}
}
void debug_print_profile (const strimap_t& profile) {
	std::cout << "--- profile ---\n";
	strimap_t::const_iterator it_p;
	for (it_p = profile.begin(); it_p != profile.end(); ++ it_p) {
		std::cout << it_p->first << "   " << it_p->second << "\n";
	}
	std::cout << "-------------\n\n";
}

void debug_print_phasing_info (const phase_profile_t& profile,
		int l_pos, const strset_t& l_vars, int r_pos,
		const strset_t& r_vars, const var_t& local_vars){
	std::cout << "LL:\n";
	debug_print_profile(profile.LL);
	std::cout << "LS:\n";
	debug_print_profile(profile.LS);
	std::cout << "SL:\n";
	debug_print_profile(profile.SL);
	std::cout << "SS:\n";
	debug_print_profile(profile.SS);

	// debug
	std::cout << l_pos << "(" << l_pos/2 << "):\t";
	for (strset_t::iterator it = l_vars.begin(); it != l_vars.end(); ++ it) {
		std::cout << "\t" << *it << "\t";
	}
	std::cout << "\n";
	std::cout << "local_vars:\t";
	var_t::const_iterator it_var = local_vars.find(r_pos);
	for (strset_t::const_iterator it = it_var->second.begin();
			it != it_var->second.end(); ++ it) {
		std::cout << *it << "\t";
	}
	std::cout << "\n------------\n";

	std::cout << r_pos << "(" << r_pos/2 << "):\t";
	for (strset_t::iterator it = r_vars.begin(); it != r_vars.end(); ++ it) {
		std::cout << "\t" << *it << "\t";
	}
	std::cout << "\nlocal_vars:\t";
	it_var = local_vars.find(r_pos);
	if (it_var != local_vars.end()) {
		for (strset_t::const_iterator it = it_var->second.begin();
				it != it_var->second.end(); ++ it) {
			std::cout << *it << "\t";
		}
		std::cout << "\n-------------\n";
	} else std::cout << "-\n------------\n";

}



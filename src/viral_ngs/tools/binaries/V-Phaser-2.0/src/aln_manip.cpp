//========================================================================
// Project     : VPhaser 2.0
// Name        : aln_manip.cpp
// Author      : Xiao Yang
// Created on  : Jul 3, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).
// Description :
//========================================================================

#include "aln_manip.h"

void debug_print_eb (std::vector<ipair_t>& errbucket){
	int sz = errbucket.size();
	std::cout << sz << "\n";
	for (int i = 0; i < sz; ++ i) {
		std::cout << errbucket[i].first << "\t" << errbucket[i].second << "\n";
	}
} // debug_print_eb

void debug_print_alnfile (ipair_t& region, std::vector<col_t>& cols) {
	std::cout << region.first << "\t" << region.second << "\n";

	/* print each column*/
	for (int i = 0; i < (int) cols.size(); ++ i) {
		int num_row = cols[i].entries.size();
		std::cout << cols[i].ref_pos << "\t" << num_row << "\n";
		for (int j = 0; j < num_row; ++ j) {
			std::cout << cols[i].entries[j].rID << "\t"
					//<< cols[i].entries[j].is_firstMate << "\t"
					//<< cols[i].entries[j].rcycle << "\t"
					//<< cols[i].entries[j].qt << "\t"
					//<< cols[i].entries[j].dt << "\t"
					<< cols[i].entries[j].eb_index << "\t"
					<< cols[i].entries[j].cons_type << "\n";
			//if (cols[i].entries[j].flanking_base)
			//	std::cout <<"\t" << cols[i].entries[j].flanking_base << "\n";
			//else std::cout << "\n";
		}
	}
}// debug_print

bool get_next_nonempty_line_buf (std::istringstream& buf, std::ifstream& iHandle) {
	std::string line;
	while (std::getline (iHandle, line)) {
		if (line.find_first_not_of(" \t\v\r\n") == std::string::npos) {
			continue;
		}
		buf.clear();
		buf.str(line);
		return true;
	}
	return false;
} // get_next_nonempty_line_buf

//
void read_disjoint_eb (eb_t& disjoint_eb, const std::string& ebfile ) {
	std::ifstream iHandle (ebfile.c_str());
	if (!iHandle.good()) {
		abording ("read_disjoint_eb: can't open file " + ebfile);
	}
    std::istringstream buf;
	read_pairs (disjoint_eb.which_pair, buf, iHandle);
	read_pairs (disjoint_eb.dt, buf, iHandle);
	read_pairs (disjoint_eb.qt, buf, iHandle);
	read_pairs (disjoint_eb.cycle, buf, iHandle);
	iHandle.close();
} // read_disjoint_eb

void read_pairs (std::vector<ipair_t>& lhs, std::istringstream& buf,
		std::ifstream& iHandle) {
	int sz ;
	if (get_next_nonempty_line_buf (buf, iHandle)) buf >> sz;
	lhs.resize(sz);
	for (int i = 0; i < sz; ++ i) {
		if (get_next_nonempty_line_buf (buf, iHandle)) {
			buf >> lhs[i].first >> lhs[i].second;
		} else abording ("read_pairs: reading disjoint eb");
	}
} // read_pairs

/* @brief	Read in [bonferroni] and [jeb], calculate and record in gParam
 * 			the err cnt and sum of snpv and lv;
 */
void read_jeb (GlobalParam& gParam, std::vector<jeb_t>& jeb) {

	std::ifstream iHandle (gParam.ebfile.c_str());
	if (!iHandle.good()) {
		abording ("read_errbucket: can't open file " + gParam.ebfile);
	}
    std::istringstream buf;
    int sz;
    if (get_next_nonempty_line_buf (buf, iHandle)) {
    		buf >> sz >> gParam.bonferroni;
    }
    gParam.bonferroni *= 2;
    jeb.resize(sz);

    imap_t pprofile, profile;
	for (int i = 0; i < sz; ++ i) {
		if (get_next_nonempty_line_buf(buf, iHandle)) {
			buf >> jeb[i].lp_err >> jeb[i].lp_sum
				>> jeb[i].snp_err >> jeb[i].snp_sum;
			gParam.bkinfo.add_sum_L (jeb[i].lp_sum);
			gParam.bkinfo.add_err_L (jeb[i].lp_err);
			gParam.bkinfo.add_sum_S (jeb[i].snp_sum);
			gParam.bkinfo.add_err_S (jeb[i].snp_err);
			if (jeb[i].snp_sum) {
				imap_t::iterator it = profile.find(jeb[i].snp_sum);
				if (it != profile.end()) ++ it->second;
				else profile[jeb[i].snp_sum] = 1;
				gParam.bkinfo.incre_bknum_S();
			}
			if (jeb[i].lp_sum) {
				imap_t::iterator it = pprofile.find(jeb[i].lp_sum);
				if (it != pprofile.end()) ++ it->second;
				else pprofile[jeb[i].lp_sum] = 1;
				gParam.bkinfo.incre_bknum_L();
			}
		} else {
			warning ("read_jeb: read insufficient entries");
			break;
		}
	}
	iHandle.close();

	gParam.bkinfo.calculate_prior();
	gParam.bkinfo.print();
} //read_jeb

/* @brief	[target_pos] specifies a sorted target set of ref positions
 * 			to consider, when empty, consider every position
 */
void read_aln_file (	ipair_t& region, std::vector<col_t>& cols,
		const ivec_t& target_pos, const std::string& alnfile) {

	int num_target = target_pos.size();

	std::ifstream iHandle (alnfile.c_str());
	if (!iHandle.good()) {
		abording ("calibrate_pe: can't open file " + alnfile);
	}

    std::istringstream buf;
    if (get_next_nonempty_line_buf (buf, iHandle)) {
		buf >> region.first >> region.second;
    }

    // identify the range in target_pos that current file contains
    int begin_target_idx = -1, end_target_idx = -1;
    if (num_target) {
    		for (int i = 0; i < num_target; ++ i) {
    			if (begin_target_idx == -1 && target_pos[i] >= region.first) {
    				begin_target_idx = i;
    			}
    			if (target_pos[i] <= region.second) {
    				end_target_idx = i;
    			} else break;
    		}
    }
    if (begin_target_idx > end_target_idx) { // nothing to be read
    		iHandle.close();
    		return;
    } else {

		while (get_next_nonempty_line_buf(buf, iHandle)) {
			col_t cur_col;
			int num;
			buf >> cur_col.ref_pos >> num;
			if (num_target) {
				if (cur_col.ref_pos != target_pos[begin_target_idx]) {
					// skip num rows
					for (int i = 0; i < num; ++ i) {
						if (!get_next_nonempty_line_buf(buf, iHandle)) {
							abording ("read_aln_file: SC failed");
						}
					}
					continue;
				}
			}

			for (int i = 0; i < num; ++ i) {
				if (get_next_nonempty_line_buf(buf, iHandle)){
					entry_t entry;
//					buf >> entry.rID >> entry.is_firstMate >> entry.rcycle
//						>> entry.qt >> entry.dt >> entry.cons_type;
//					if (entry.dt.at(0) == 'D' ||
//							(entry.dt.size() > 1 && entry.dt.at(1) == 'D')) {
//						buf >> entry.flanking_base;
//					}
					buf >> entry.rID >> entry.eb_index >> entry.cons_type >> entry.is_rv;
					if (*entry.cons_type.rbegin() == 'D') {
						entry.cons_type = entry.cons_type.substr(
								0, entry.cons_type.size() - 1);
						entry.is_D_entry = true;
					}
					cur_col.entries.push_back(entry);
				} else break;
			}
			cols.push_back(cur_col);

			if (num_target) {
				++ begin_target_idx;
				if (begin_target_idx > end_target_idx) {
					iHandle.close();
					return;
				}
			}
		}// while

		iHandle.close();
	}

  //	debug_print_alnfile (region, cols);
  //exit(1);
}//read_aln_file


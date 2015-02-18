//========================================================================
// Project     : VPhaser 2.0
// Name        : aln_manip.h
// Author      : Xiao Yang
// Created on  : Jul 3, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================


#ifndef ALN_MANIP_H_
#define ALN_MANIP_H_

#include <iostream>
#include "xutil.h"
#include "format.h"

/* This structure initially stores elements read from alingment file.
 * Then some entries are computed and used to store new values as noted
 */
struct entry_t{
//	bool is_firstMate;
//	char flanking_base;
//	int rcycle;
//	int qt;
//	std::string dt;

	bool is_D_entry;
	bool is_cons_norm;
	bool is_cons_poly;
	bool is_rv;
	int	rID;
	int eb_index;
	std::string cons_type;

//	double pe;

	entry_t () {
		//pe = 0.0;
		is_D_entry = false;
		is_cons_norm = false;
		is_cons_poly = false;
	}
};


// ref_pos, [entry0, 1, 2...]
// let consensus base to be the alphabetically smallest one when
// breaking the equality; only exception is for insertion column,
// insertion will be considered as variants against non-insertions,
// and will be reported as variant even if it is dominant
typedef struct {
	int ref_pos;
	//std::string consensus;
	std::vector<entry_t> entries;
} col_t;

//typedef std::vector<std::pair<int, std::string> > read_t;


//void read_aln_file (ipair_t& region, std::vector<col_t>& cols,
//		read_t& reads, const std::string& alnfile);
void read_aln_file (ipair_t& region, std::vector<col_t>& cols,
		const ivec_t& target_pos, const std::string& alnfile);

void read_jeb (GlobalParam& gParam, std::vector<jeb_t>& jeb);

void read_disjoint_eb (eb_t& disjoint_eb, const std::string& ebfile);
void read_pairs (std::vector<ipair_t>& lhs, std::istringstream& buf,
		std::ifstream& iHandle);

void col_profile (int& max_idx, ivec_t& profile, const std::vector<entry_t>& col);

#endif /* ALN_MANIP_H_ */

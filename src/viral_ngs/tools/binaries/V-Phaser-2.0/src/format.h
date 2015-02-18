//========================================================================
// Project     : VariantCaller
// Name        : format.h
// Author      : Xiao Yang
// Created on  : Jul 3, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================


#ifndef FORMAT_H_
#define FORMAT_H_

#include "api/BamReader.h"
#include "xutil.h"
#include "Parameter.h"
#include "xny/EnumParser.hpp"
#include "bam_manip.h"

using namespace BamTools;

/* joint error bucket */
struct jeb_t {
	uint32_t snp_err;
	uint32_t snp_sum;
	uint32_t lp_err;
	uint32_t lp_sum;
	jeb_t () {
		snp_err = 0;
		snp_sum = 0;
		lp_err = 0;
		lp_sum = 0;
	}
};


/* independent error bucket */
struct eb_t {
	std::vector<ipair_t> which_pair;
	std::vector<ipair_t> cycle;
	std::vector<ipair_t> qt;
	std::vector<ipair_t> dt;
	eb_t () {
		which_pair = std::vector<ipair_t> (2, ipair_t (0, 0));
		dt = std::vector<ipair_t> (36, ipair_t (0, 0));
	}
};

inline void debug_print_jeb (const std::vector<jeb_t>& jeb) {
	for (int i = 0; i < (int) jeb.size(); ++ i) {
		if (jeb[i].lp_sum != 0 || jeb[i].snp_sum != 0) {
			std::cout << i << ": \t";
			std::cout << jeb[i].lp_err << ",\t" << jeb[i].lp_sum << ",\t"
				<< jeb[i].snp_err << ",\t" << jeb[i].snp_sum << ",\t";
			if (jeb[i].lp_sum != 0) {
				std::cout << (jeb[i].lp_err * 1.0)/(jeb[i].lp_sum * 1.0) << ",\t";
			} else std::cout << " - ,\t";
			if (jeb[i].snp_sum != 0) {
				std::cout << (jeb[i].snp_err * 1.0)/(jeb[i].snp_sum * 1.0) << "\n\n";
			} else std::cout << "\n\n";
		}
	}
} // debug_print_jeb



void prep_aln_file (GlobalParam& gParam, const Parameter& myPara);


// ---------------- disjoint error probability --------------
void update_disjoint_eb_entry (eb_t& disjoint_eb, int isp0, int cycle,
		int dt, 	int qt,	bool is_cons);
void update_pairs (std::vector<ipair_t>& lhs, int index, int incre);

double calculate_disjoint_pe (int is_p0, int cycle,
		int dt, 	int qt, const eb_t& disjoint_eb);

// ---------------- joint err probability ------------------

void update_jeb (std::vector<jeb_t>& jeb, const std::vector<AlnColEntry>& col,
  bool is_ins, bool is_pre_ins, const Parameter& myPara, const GlobalParam& gParam);
void update_D_entry_dt (std::string& dt, const AlnColEntry& entry) ;
void resize_jeb (std::vector<jeb_t>& jeb, int idx);

void update_eb_entry (std::vector<ipair_t>& errBuckets, int index,
		int incre_mismatch, int incre_total);
int calculate_eb_index ( int isp0, int c, int maxRL, int di, int maxDI,
		int qt, const Parameter& myPara);

// ---------------- profile an alignment column --------------
void profiling (bool& is_diverse, std::vector<AlnColEntry>& col, bool is_ins);
std::string consensus_type (const AlnColEntry& entry);
void add_to_profile (std::map<std::string, int>& profile,
		const std::string& type);
void add_to_profile_strand (std::map<std::string, ipair_t>& profile,
		const std::string& type, bool is_rv);
void analyze_profile (std::string& cons, ivec_t& cnts,
		const std::map<std::string, int>& profile);

// ---------------- write to file ----------------------------
bool create_output_file (std::ofstream& oHandle, const std::string& filename);
void generate_mate_map (std::ofstream& oHandle, imap_t& mate_map,
		const std::map<std::string, ipair_t>& mateInfo);
void write_col (std::ofstream& oHandle, const std::vector<AlnColEntry>& col,
		bool is_ins, const imap_t& mate_map, int ref_pos,
		const Parameter& myPara, const GlobalParam& gParam);
void write_jeb (std::ofstream& oHandle, const std::vector<jeb_t>& jeb,
		uint32_t bonferroni_factor);
void write_disjoint_eb (std::ofstream& oHandle, const eb_t& disjoint_eb);
void write_pairs (std::ofstream& oHandle, const std::vector<ipair_t>& lhs);

void generate_covplot_Rscript (const ivec_t& cov, const std::string& path,
		const std::string& refName);

void get_mateInfo (std::map<std::string, ipair_t>& mateInfo,
	std::vector<BamAlignment>& alns, Ref::const_iterator& it_ref,
	GlobalParam& gParam, const Parameter& myPara);
void mate_to_fetch (
	std::vector<std::pair<std::string, MapRecord> >& to_fetch,
	const std::map<std::string, ipair_t>& mateInfo,
	const std::vector<BamAlignment>& alns,
	const GlobalParam& gParam, int delta);
void get_mate_locus (std::map<std::string, ivec_t>& ref_locus,
		const std::vector<std::pair<std::string, MapRecord> >& to_fetch,
		const Ref& reference);
bool get_refName (std::string& refName, int fileId, int refId,
		const Ref& reference);
void get_mates (std::vector<BamAlignment>& alns,
	const std::map<std::string, ivec_t>& ref_locus,
	const std::vector<std::pair<std::string, MapRecord> >& to_fetch,
	const GlobalParam& gParam, const Parameter& myPara);

void update_mateInfo (std::map<std::string, ipair_t>& mateInfo,
		std::vector<BamAlignment>& alns,
		const std::vector<BamAlignment>& mate_alns);

void get_aln_info (ivec_t& aln_ends, std::vector<ipair_t>& ranges,
		std::vector<BamAlignment>& alns,
		const std::map<std::string, MapEntry>& alnMap);

// unused function
void merge_readpairs (std::vector<BamAlignment>& alns, const ivec_t& aln_ends,
		const std::map<std::string, ipair_t>& mateInfo);

struct cmp_aln {
	bool operator () (const BamAlignment& lhs, const BamAlignment& rhs) {
		return lhs.Position < rhs.Position;
	}
};

struct cmp_alnentry{
	bool operator () (const AlnColEntry& lhs, const AlnColEntry& rhs) {
			return lhs.idx_aln_array < rhs.idx_aln_array;
	}
};

#endif /* FORMAT_H_ */

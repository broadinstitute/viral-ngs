//========================================================================
// Project     : VariantCaller
// Name        : bam_manip.h
// Author      : Xiao Yang
// Created on  : Apr 9, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================


#ifndef BAM_MANIP_H_
#define BAM_MANIP_H_

#include "api/BamReader.h"
#include "xutil.h"
#include "Parameter.h"
#include "jaz/string_add.hpp"
#include "jaz/fastx_iterator.hpp"
#include "xny/EnumParser.hpp"
#include "xny/seq_manip.hpp"

using namespace BamTools;

struct AlnColEntry{
	//	int idx_cigar;	   // index of cigar string
	//	int	cigar_pos;	   // the position of cigar
	int	idx_aln_array; // index of array that stores alignments
	int rcycle;	  	   // read cycle
	int qt; 			   // quantile of quality score
	bool is_firstMate;
	bool is_rv;
	char flanking_base;	// for dt = IxM or MDx or DxM--> fIxM or MDxf or fDxM
	bool is_cons_norm;		// if normal type is consensus
	bool is_cons_poly;  // if length polymorphic type is cons
	std::string dt;
	std::string cons_type; // consensus type: for IxM initially register
						   // as x inserted bases, then register as I+xbases
	AlnColEntry () {
		flanking_base = '\0';
		dt = cons_type = "";
		is_cons_norm = is_cons_poly =  false;
	};
};

inline void debug_print_alnentry (const std::vector<AlnColEntry>& col) {
	for (int d = 0; d < (int) col.size(); ++ d) {
		std::cout << col[d].idx_aln_array << ": " << col[d].dt << " , ";
		if (col[d].is_rv) std::cout << "- , ";
		else std::cout << "+ , ";
		std::cout << col[d].flanking_base << " , " << col[d].cons_type
				<< " , "	 << col[d].qt << ", " << col[d].is_cons_norm
				<< ", " << col[d].is_cons_poly << "\n";
	}
	std::cout << "\n\n";
}



/* reverse complementary aln */
struct rvc_aln_t {
	std::string raw_read;
	std::string raw_qual;
};

// identify the min and max quality scores using p% of input reads
void sampling (int& minQ, int& maxQ, int& maxRL, int& avgfragSz,
		int& stdfragSz, const strvec_t& filelist, int P);

void qqMap (imap_t& qq, int minQ, int maxQ, int quantile);

void parse_bam_header (GlobalParam& gParam, int& pSample);

void set_rmap_array (GlobalParam& gParam, const Parameter& myPara);

bool fetch_ref_seq (Ref& reference, const char* filepath) ;

// unused function
//void update_errBucket (std::vector<ipair_t>& errBuckets,
//		int& idx_ref, int maxRL, int maxQt, const BamAlignment& al,
//		const std::string& refseq, EnumParser<DiNt>& parser,
//		const imap_t& qq);


void udpate_alnMap (std::map<std::string, MapEntry>& alnMap,
		int fileId, const BamAlignment& al);
int calculate_aln_end (const BamAlignment& al) ;


void gather_alignments (std::vector<BamAlignment>& alns,
		Ref::const_iterator& it_ref, int start, int end,
		const GlobalParam& gParam, const Parameter& myPara);

void compute_rvc_alns (std::vector<rvc_aln_t>& rvc_alns,
		const std::vector<BamAlignment>& alns);

void get_aln_column (std::vector<AlnColEntry>& cur_col, int ref_pos,
	int ignbases, const std::vector<BamAlignment>& alns,
	const ivec_t& aln_ends,  const imap_t& qq);
void get_column_x (AlnColEntry& col, int ref_pos, const BamAlignment& al,
		int end, const imap_t& qq);

void get_ins_column (std::vector<AlnColEntry>& ins_col,
		const std::vector<AlnColEntry>& cur_col,
		const std::vector<AlnColEntry>& prev_col);
bool add_ins_entry (std::vector<AlnColEntry>& col, const AlnColEntry& entry);

void clean_aln_column (std::vector<AlnColEntry>& col);

void get_extended_cigar_string (cvec_t& ext_cigar,
		const std::vector<CigarOp>& cigar_data);
#endif /* BAM_MANIP_H_ */

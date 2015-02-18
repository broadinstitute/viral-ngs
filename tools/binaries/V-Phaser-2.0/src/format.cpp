//========================================================================
// Project     : VariantCaller
// Name        : format.cpp
// Author      : Xiao Yang
// Created on  : Jul 3, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================

#include "format.h"

 /* @brief	Given fileId and refId, identify the reference name
 * 			Return false if refName couldn't be found
 */
bool get_refName (std::string& refName, int fileId, int refId,
		const Ref& reference) {
	Ref::const_iterator it_r = reference.begin();
	for (; it_r != reference.end(); ++ it_r) {
		for (int i = 0; i < (int) it_r->second.bfID_refID.size(); ++ i) {
			if (it_r->second.bfID_refID[i] == ipair_t (fileId, refId)) {
				refName = it_r->first;
				return true;
			}
		}
	}
	return false;
} // id_refName

/* @brief	Identify which mates need to be fetched from bam file. This
 *			info is stored in [to_fetch].
 * @method	Identify if both mates have been already identified according
 * 			to [mateInfo]. If not, try to identify the mate according to
 * 			[gParam.alnMap] if the mate is mapped to the same ref with a
 * 			larger coordinate.
 * 			Note: currently we do NOT consider when pair mapped to different ref
 */
void mate_to_fetch (
	std::vector<std::pair<std::string, MapRecord> >& to_fetch,
	const std::map<std::string, ipair_t>& mateInfo,
	const std::vector<BamAlignment>& alns, const GlobalParam& gParam,
	int delta){

	std::map<std::string, ipair_t>::const_iterator it_mi = mateInfo.begin();
	for (; it_mi != mateInfo.end(); ++ it_mi) {

		// the mate is not already loaded in mateInfo
		if (it_mi->second.second == -1) {

			std::string rName = it_mi->first;
			int idx_0 = it_mi->second.first;
			int fileId_mate_0 = get_index<std::string> (
							alns[idx_0].Filename, gParam.bam_filelist);
			if (fileId_mate_0 == -1) {
				abording ("id_paired_alignments: can't find file "
											+ alns[idx_0].Filename);
			}
			int refId_mate_0 = alns[idx_0].RefID;
			bool is_firstMate = alns[idx_0].IsFirstMate();

			// search the read in the global map
			std::map<std::string, MapEntry>::const_iterator
				it_m = gParam.alnMap.find(rName);

			if (it_m == gParam.alnMap.end()) {
				abording ("id_paired_alignments: failed to find read " +	rName);
			}
			MapRecord  m0 = it_m->second.first,
 					   m1 = it_m->second.second;

			// the mate is mapped in the global map
			// For current pair m0, we only consider pair m1 that mapped after m0
			// if both mapped to the same reference
			if (m1.start != -1) {
				// mate_0 is the first entry in global read map [gParam.alnMap]
				if (fileId_mate_0 == m0.fileID && refId_mate_0 == m0.refID
						&& is_firstMate == m0.isfirstMate) {

					int dist = m1.start - m0.start;
					// distance constraint
					if (dist <= gParam.avgFragSz + gParam.stdFragSz*delta) {
						to_fetch.push_back(std::pair<std::string, MapRecord> (rName, m1));
					}

				} // mate_0 is the second entry in global read map; we don't consider its
				  // pair that mapped to a smaller coordinate on the same ref
				else if (fileId_mate_0 == m1.fileID &&
						refId_mate_0 == m1.refID &&	is_firstMate == m1.isfirstMate) {

					/* currently we do NOT consider when pair mapped to different ref
						if (fileId_mate_0 != m0.fileID || refId_mate_0 != m0.refID) {
							//the mate pair is mapped to diff fileId or refId
							to_fetch.push_back(std::pair<std::string, MapRecord>
											   (rName, m0));
						}
					*/
				} else {
					abording ("id_paired_alignments: the read " + rName +
							" is not registered in the global read map");
				}
			} // if (m1.start != -1)

		} // if (it_mi->second.second == -1)

	} // for (it_mi =
}//mate_to_fetch

/* @brief
 *
 */
void get_mate_locus (std::map<std::string, ivec_t>& ref_locus,
		const std::vector<std::pair<std::string, MapRecord> >& to_fetch,
		const Ref& reference){

	for (int i = 0; i < (int) to_fetch.size(); ++ i) {
		std::string refName;
		if (! get_refName (refName, to_fetch[i].second.fileID,
				to_fetch[i].second.refID, reference)){
			abording ("id_paired_alignments: can't id refName");
		}
		std::map<std::string, ivec_t>::iterator
		it_rl = ref_locus.find(refName);

		if (it_rl == ref_locus.end()) {
			ref_locus[refName] = ivec_t(1, i);
		} else {
			bool is_overlap = false;
			for (int j = 0; j < (int) it_rl->second.size(); ++ j) {
				// compare if current read aln overlaps w/ existing ones
				int idx = it_rl->second[j]; // existing idx
				if ((to_fetch[i].second.start <= to_fetch[idx].second.end) &&
					(to_fetch[i].second.end >= to_fetch[idx].second.start)) {

					is_overlap = true;
					break;
				}
			} // for
			if (! is_overlap) { // add a new aln range
				it_rl->second.push_back(i);
			}
		} // else
	} // for (int i

} // get_mate_locus

/* @brief	Fetch all alignments aligned to each locus on ref according
 * 			to [ref_locus], then keep only those alignments recorded in
 * 			[to_fetch], the results will be stored in [alns]
 *
 */
void get_mates (std::vector<BamAlignment>& alns,
	const std::map<std::string, ivec_t>& ref_locus,
	const std::vector<std::pair<std::string, MapRecord> >& to_fetch,
	const GlobalParam& gParam, const Parameter& myPara) {

	bool is_debug_print = false;

	alns.clear();

	/* gather all reads that needs to be fetched, stored as
	 * rName --> index on to_fetch */
	std::map<std::string, int> rNames;
	for (int i = 0; i < (int)  to_fetch.size(); ++ i) {
		rNames[to_fetch[i].first] = i;
	}

	if (is_debug_print){
		std::cout << "\t\t" << rNames.size() << " mates to be fetched\n";
	}

	std::map<std::string, ivec_t>::const_iterator
						it_rl = ref_locus.begin();
	for (; it_rl != ref_locus.end(); ++ it_rl) {
		std::string ref_name = it_rl->first;

		/* gather alignments for each genomic region */
		for (int i = 0; i < (int) it_rl->second.size(); ++ i) {
			int idx = it_rl->second[i];
			int window_start = to_fetch[idx].second.start,
				window_end = to_fetch[idx].second.end;

			/* stores temporary alignments wrt a particular ref region */
			std::vector<BamAlignment> tmp_alns;

			Ref::const_iterator it_ref = gParam.reference.find(ref_name);

			if (it_ref != gParam.reference.end()){

				gather_alignments (tmp_alns, it_ref, window_start,
						window_end, gParam, myPara);

				/*
				{	// debug print out all alignments
					std::cout << "wanted reads: " << it_rl->second[i].size() << "\n";
					std::cout << tmp_alns.size() << " alns between "
							<< window_start << "\t"	<< window_end << "\n\n";
					for (int i = 0; i < (int) tmp_alns.size(); ++ i) {
						std::cout << tmp_alns[i].Name << "\t"
								<< tmp_alns[i].Position << "\n";
					}
					std::cout << "\n";
					exit(1);
				} */

				for (int j = 0; j < (int) tmp_alns.size(); ++ j) {

					std::map<std::string, int>::iterator it_rn =
							rNames.find(tmp_alns[j].Name);

					if (it_rn != rNames.end()) {
						idx = it_rn->second;
						if (to_fetch[idx].second.isfirstMate ==
								tmp_alns[j].IsFirstMate()) {
							alns.push_back(tmp_alns[j]);
							rNames.erase(it_rn);
						}
					}
				} // for (int j
			} else {
				abording ("\tget_mates: can't identify ref: " + ref_name);
			}
		} // for (int i = 0
	} // for (; it_rl

	/*
	{ // debug print out
		std::cout << "before\n";
		for (int i = 0; i < alns.size(); ++ i) {
			std::cout << alns[i].Position << "\t";
			if (i != 0 && i % 20 == 0) std::cout << "\n";
		}
		std::cout << "--------------------\n\n";
	} */

	std::sort (alns.begin(), alns.end(), cmp_aln());
	if (is_debug_print){
		std::cout << "\t\t" << alns.size() << " mates fetched\n\n";
	}

} // get_mates

/* @brief	Append [mate_alns] to array [alns], and update [mateInfo]
 */
void update_mateInfo (std::map<std::string, ipair_t>& mateInfo,
		std::vector<BamAlignment>& alns,
		const std::vector<BamAlignment>& mate_alns){

	int start = alns.size();

	alns.insert(alns.end(), mate_alns.begin(), mate_alns.end());

	for (int i = start; i < (int) alns.size(); ++ i) {
		std::map<std::string, ipair_t>::iterator
			it = mateInfo.find(alns[i].Name);
		if (it == mateInfo.end()) {
			abording ("update_mateInfo: cannot find read " + alns[i].Name);
		} else {
			it->second.second = i;
		}
	} // for (int i = start
} // update_mateInfo

/* @brief	Given vector of alignments [alns], sorted according to the
 * 			alignment position wrt the reference sequence [it_ref],
 * 			identify mate pair information for any read involved in [alns].
 * 			When the mate alignment is not yet in [alns], they will be
 * 			retrieved from bam files. The final mate pair information is
 * 			recorded in [mateInfo].
 * @param
 * 	[mateInfo]: 	scan [alns], generate a map <readName, (idx_0, idx_1)>,
 * 	where idx_i	is the index of [alns] where readName maps. Initially,
 * 	idx_1 = -1,	it will be updated when the mate pair can be identified
 * 	read with idx_1 map to a larger coordinate on ref compared to idx_0
 */
void get_mateInfo (std::map<std::string, ipair_t>& mateInfo,
	std::vector<BamAlignment>& alns, Ref::const_iterator& it_ref,
	GlobalParam& gParam, const Parameter& myPara) {

	/* register information for mate pairs that already loaded */
	for (int i = 0; i < (int) alns.size(); ++ i) {
		std::map<std::string, ipair_t>::iterator
				it_mi = mateInfo.find(alns[i].Name);
		if (it_mi == mateInfo.end()) { // read is not found
			mateInfo[alns[i].Name] = ipair_t (i, -1);
		} else { // found the read, update the second index
			it_mi->second.second = i;
		}
	}

	/*
	{ // debug print mateinfo
		int cnt_equal_start = 0;
		std::map<std::string, ipair_t>::iterator
					it_m = mateInfo.begin();
		for (; it_m != mateInfo.end(); ++ it_m) {
			//std::cout << it_m->first << "\t";
			int idx_0 = m0, idx_1 = m1;

			if (idx_1 != -1) {
				//std::cout << alns[idx_0].Position
				//	<< "\t" << alns[idx_1].Position << "\n";
				if (alns[idx_0].Position == alns[idx_1].Position) {
					cnt_equal_start ++;
				}
			}
		}
		std::cout << "matepairs w/ same start:" << cnt_equal_start << "\n";
		std::cout << "\n\n\n";
		//exit(1);
	}*/

	/* 	scan [mateInfo], identify the mate pair that may need to be
	 *  loaded, given [gParam.alnMap]. The mate either has a larger
	 *  coordinate on the same ref or is located on a different ref
 	 */

	// the vector to store the mate that will be retrieved from bam files
	std::vector<std::pair<std::string, MapRecord> > to_fetch;
	mate_to_fetch (to_fetch, mateInfo, alns, gParam, myPara.delta);

	/*
	{ // debug print to_fetch vector
		std::cout << "alignments to fetch:\n";
		for (int i = 0; i < (int) to_fetch.size(); ++ i) {
			std::cout << to_fetch[i].first << "\t" << to_fetch[i].second.start
					<< "\t" << to_fetch[i].second.end << "\n";
		}
		exit(1);
	}*/

	/* Given a set of read information, pinpoint on ref genome, which
	 * regions to be used for retrieving the alignments.
	 * @method: a) select a read r from [to_fetch], discard all reads that
	 * overlaps with r; b) repeat a) until all reads  in [to_fetch] have
	 * been accounted for. The results is stored in
	 * [ref_locus]: storing a map from refName to indices of [to_fetch],
	 * where each index denotes a selected read in @method.
	 */
	std::map<std::string, ivec_t> ref_locus;
	get_mate_locus (ref_locus, to_fetch, gParam.reference);

	/* fetch the actual alignments, sorted in increasing order of starting
	 * position of mapping  */
	std::vector<BamAlignment> mate_alns;
	get_mates (mate_alns, ref_locus, to_fetch, gParam, myPara);

	// Append [mate_alns] to array [alns], and update [mateInfo]
	update_mateInfo (mateInfo, alns, mate_alns);

} // get_mateInfo

/* @brief	For each alignment in [alns] identify ending position
 * 			on ref and store in vector [aln_ends]; [ranges] stores
 * 			the start and end positions of ref regions covered by reads,
 * 			i.e. the positions to retrieve alignment columns.
 */
void get_aln_info (ivec_t& aln_ends, std::vector<ipair_t>& ranges,
	std::vector<BamAlignment>& alns, const std::map<std::string, MapEntry>& alnMap) {

	int sz = alns.size();
	std::map<std::string, MapEntry>::const_iterator it_m;
	for (int i = 0; i < sz; ++ i) {
		int start, end;
		it_m = alnMap.find(alns[i].Name);
		if (it_m != alnMap.end()) {
			if (it_m->second.first.isfirstMate == alns[i].IsFirstMate()) {
				start = it_m->second.first.start;
				end = it_m->second.first.end;
			} else if (it_m->second.second.isfirstMate == alns[i].IsFirstMate()){
				start = it_m->second.second.start;
				end = it_m->second.second.end;
			} else {
				std::cout << it_m->second.first.isfirstMate << ", "
						<< it_m->second.second.isfirstMate << ", "
						<< alns[i].IsFirstMate() << "\n";
				abording ("get_aln_ends: cannot find read " + alns[i].Name);
			}
			aln_ends.push_back(end);
			// merge [start, end] with ranges
			/*
			{ /////
			std::cout << "(" << start << ", " << end << ")\n";
			}*/
			if (!ranges.size()) ranges.push_back(ipair_t(start, end));
			else {
				bool is_inserted = false;
				std::vector<ipair_t>::iterator it_r = ranges.begin();
				for (; it_r != ranges.end(); ++ it_r) {
					if (end < it_r->first) {
						ranges.insert(it_r, 1, ipair_t(start, end));
						is_inserted = true;
						break;
					} else if ((start >= it_r->first && start <= it_r->second)
						|| (end >= it_r->first && end <= it_r->second)) {
						it_r->first = std::min(start, it_r->first);
						it_r->second = std::max(end, it_r->second);
						is_inserted = true;
						break;
					}
				}
				if( !is_inserted) ranges.push_back(ipair_t(start, end));
				/*{ ///////
				std::cout << "ranges: \t";
				for (int i = 0; i < (int) ranges.size(); ++ i) {
					std::cout << ranges[i].first << "\t" << ranges[i].second << "\n";
				}
				}*/

			}
		 } else {
			 abording ("get_aln_ends: cannot find read " + alns[i].Name);
		 }
	}
	/*
	{ // debug print aln_ends
		for (int i = 0; i < sz; ++ i) {
			std::cout << aln_ends[i] << "\t";
			if (i != 0 && i % 20 == 0) std::cout << "\n";
		}
		std::cout << "\n\n";
		exit(1);
	} */

} // get_aln_info

/* @brief	Identify all reads pairs that the forward strand overlaps with
 * 			the reverse strand. Merge both and modify entries in [alns]
 * @method	For each position in the overlapping region, If two bases
 * 			differ, retain the base with a higher quality, and substract
 * 			the quality score. Prefer base over gaps.
 *
 * 	This function is not used currently
 */
void merge_readpairs (std::vector<BamAlignment>& alns, const ivec_t& aln_ends,
		const std::map<std::string, ipair_t>& mateInfo){

	std::map<std::string, ipair_t>::const_iterator it_mi = mateInfo.begin();
	for (; it_mi != mateInfo.end(); ++ it_mi) {
		int idx_0 = it_mi->second.first, idx_1 = it_mi->second.second;
		if (idx_1 != -1 && aln_ends[idx_0] >= alns[idx_1].Position){
			// identify overlaps wrt reference
			int start = std::max(alns[idx_0].Position, alns[idx_1].Position),
				end = std::min(aln_ends[idx_0], aln_ends[idx_1]);


			if (alns[idx_0].Name.compare("A07N4111122:1:11:18392:9166") == 0) {
				std::cout << "here\n";
			}

			{ // debug print info
				if (alns[idx_0].Name.compare(alns[idx_1].Name) != 0) {
					abording("merge_readpairs p0 SC failed");
				} else {
					std::cout << "\n" <<  alns[idx_0].Name << "\t";
					std::cout << "idx_0, 1 = " << idx_0 << "," << idx_1 << "\n";
					std::cout << "len = " << alns[idx_0].AlignedBases.length()
							<< "," << alns[idx_1].AlignedBases.length() << "\n";
				}

				std::cout << " starts, ends: " << alns[idx_0].Position << ","
						<< aln_ends[idx_0]<< "," << alns[idx_1].Position
						<< "," <<  aln_ends[idx_1] << "\n";
				std::cout << start << "\t" << end << "\n";
			}

			cvec_t ext_cigar_0, ext_cigar_1;
			get_extended_cigar_string(ext_cigar_0, alns[idx_0].CigarData);
			get_extended_cigar_string(ext_cigar_1, alns[idx_1].CigarData);


			std::string aln0, aln1;
			int r0_refpos = alns[idx_0].Position;
			int r0_idx_extcigar = 0;
			for (; r0_idx_extcigar < (int) ext_cigar_0.size(); ++ r0_idx_extcigar) {
				if (r0_refpos == start) break;
				if (ext_cigar_0[r0_idx_extcigar] == 'I') continue;
				++ r0_refpos;
			}
			aln0 = alns[idx_0].AlignedBases.substr(0, r0_idx_extcigar);
			aln1 = std::string (r0_idx_extcigar, '_');

			int r1_idx_extcigar = 0;
			for (; r0_idx_extcigar < (int) ext_cigar_0.size() &&
				   r1_idx_extcigar < (int) ext_cigar_1.size();
					++ r0_idx_extcigar, ++ r1_idx_extcigar){
				char c0 = ext_cigar_0.at(r0_idx_extcigar),
					 c1 = ext_cigar_1.at(r1_idx_extcigar);
				if (c0 == 'I' && c1 != 'I') {
					aln0 += alns[idx_0].AlignedBases.at(r0_idx_extcigar);
					aln1 += '_';
					-- r1_idx_extcigar;
				} else if (c0 != 'I' && c1 == 'I') {
					aln0 += '_';
					aln1 += alns[idx_1].AlignedBases.at(r1_idx_extcigar);
					-- r0_idx_extcigar;
				} else {
					aln0 += alns[idx_0].AlignedBases.at(r0_idx_extcigar);
					aln1 += alns[idx_1].AlignedBases.at(r1_idx_extcigar);
				}
			}


			std::cout << aln0 << "\n" << aln1 << "\n";
			std::cout << "\n";


			{
				if (alns[idx_0].Name.compare("A07N4111122:1:11:18392:9166") == 0) {
					for (int d = 0; d < ext_cigar_0.size(); ++ d) {
						std::cout << ext_cigar_0[d] ;
					}
					std::cout << "\n";
					for (int d = 0; d < ext_cigar_1.size(); ++ d) {
						std::cout << ext_cigar_1[d] ;
					}
					std::cout << "\n";
					exit(1);
				}
			}

		}
	}
	/*
	{ // debug print out mate pair information
		// --->       --->    --->     --->
		// <---     <---        <---         <---
		int cnt_equal = 0, cnt_larger = 0, cnt_overlap = 0, cnt_no = 0;
		std::map<std::string, ipair_t>::iterator it_mi = mateInfo.begin();
		for (; it_mi != mateInfo.end(); ++ it_mi) {
			int idx_0 = it_mi->second.first, idx_1 = it_mi->second.second;

			std::cout << it_mi->first << ":\t\t";
			if (idx_1 != -1) {
				std::cout << alns[idx_0].Position
							<< "\t\t" << alns[idx_1].Position << "\t";
				if (alns[idx_0].Position == alns[idx_1].Position){
					cnt_equal ++;
					std::cout << "+++++++++++++++++++++++++++++++++\n";
				} else if (alns[idx_0].Position > alns[idx_1].Position){
					cnt_larger ++;
					std::cout << "---------------------------------\n";
				} else if ((alns[idx_0].Position < alns[idx_1].Position) &&
						aln_ends[idx_0] >= alns[idx_1].Position) {
					cnt_overlap ++;
					std::cout << aln_ends[idx_0] - alns[idx_1].Position + 1 << "\t";
					std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
				} else {
					cnt_no ++;
					std::cout << "\n";
				}
			} else {
				std::cout << alns[idx_0].Position << "\n";
			}
		}
		std::cout << "mate_pair positions: equal, larger, overlap, not-overlap: " <<
				cnt_equal << "," << cnt_larger << "," << cnt_overlap << ","
				<< cnt_no << "\n";
		std::cout << "\n";
		exit(1);
	}*/
} // merge_readpairs

/* @brief	return false if the output file already exists
 */
bool create_output_file (std::ofstream& oHandle, const std::string& filename) {
	oHandle.clear();
	/* check if file exists already */
	// [undo this]
	/*
	std::fstream fin;
	fin.open(filename.c_str());
	if (fin.is_open()) {
		fin.close();
		std::cout << "\n\t\tfile: " << filename << " exists !\n\n";
		return false;
	}*/

	std::cout << "\n\t\tcreate file: " << filename << "\n";
	oHandle.open(filename.c_str());
	if (!oHandle.is_open()) {
		abording ("create_output_file: can't open file " + filename);
	}
	return true;
} // create_output_file

/* @brief	For each entry of mateInfo: rName->(p0, p1), assuming p1 > p0,
 * 			write to file "p0 \t rName" and
 * 			generate a map [mate_map]: p1 -> p0, which is used when output
 * 			aln columns to file, p1 will be recored as p0 to count for
 * 			read pair info
 */
void generate_mate_map (std::ofstream& oHandle, imap_t& mate_map,
		const std::map<std::string, ipair_t>& mateInfo) {

//	oHandle << mateInfo.size() << "\n"; // disable writing read name
	std::map<std::string, ipair_t>::const_iterator it_mi (mateInfo.begin());
	for (; it_mi != mateInfo.end(); ++ it_mi) {
		int lower = it_mi->second.first;
		if (it_mi->second.second != -1) {
			int larger = it_mi->second.second;
			if (lower > larger){
				lower = larger;
				larger = it_mi->second.first;
			}
			mate_map[larger] = lower;
		}
//		oHandle << lower << "\t" <<it_mi->first << "\n";  // disable writing read name
	}
}// generate_mate_map

/* @brief	Writing col information to file
 */
void write_col (std::ofstream& oHandle, std::vector<AlnColEntry>& col,
		bool is_ins, const imap_t& mate_map, int ref_pos,
		const Parameter& myPara, const GlobalParam& gParam) {

	EnumParser<DiNt> parser;

	int sz = col.size();
	for (int i = 0; i < sz; ++ i) {
		// replace the larger ID of the paired read to the smaller one
		// using mate_map information
		imap_t::const_iterator it = mate_map.find(col[i].idx_aln_array);
		if (it != mate_map.end()) col[i].idx_aln_array = it->second;
	} // for (int i = 0

	std::sort(col.begin(), col.end(), cmp_alnentry());

	if (is_ins) oHandle << (2*ref_pos - 1) << "\t" << sz << "\n";
 	else oHandle << 2*ref_pos << "\t" << sz << "\n";

	for (int i = 0; i < sz; ++ i) {

		std::string dt = col[i].dt;

		// DxM or MDx, update dt to exclude 'Dx' and include flanking base
		update_D_entry_dt (dt, col[i]);

		int idx_jeb = calculate_eb_index (col[i].is_firstMate,
				col[i].rcycle, gParam.maxRL, parser.parseEnum(dt),
				DiNtSize, col[i].qt, myPara);

//		oHandle << col[i].idx_aln_array << "\t" << col[i].is_firstMate
//				<< "\t" << col[i].rcycle << "\t" << col[i].qt << "\t"
//				<< col[i].dt << "\t" << col[i].cons_type;
//		if (!is_ins && col[i].flanking_base != '\0') {
//			oHandle << "\t" << col[i].flanking_base << "\n";
//		}
//		else oHandle << col[i].flanking_base << "\n";

		oHandle << col[i].idx_aln_array << "\t" << idx_jeb
			<< "\t" << col[i].cons_type;
		if (col[i].dt.at(0) == 'D') oHandle << "D";
		oHandle << "\t" << col[i].is_rv << "\n";
	}
} //write_col

//
void write_jeb (std::ofstream& oHandle, const std::vector<jeb_t>& jeb,
		uint32_t bonferroni_factor) {
	int sz = jeb.size();
	oHandle << sz << "\t" << bonferroni_factor << "\n";
	for (int i = 0; i < sz; ++ i) {
		oHandle << jeb[i].lp_err << "\t" << jeb[i].lp_sum << "\t"
				<< jeb[i].snp_err << "\t" << jeb[i].snp_sum << "\n";
	}
} // write_errbucket

//
void write_disjoint_eb (std::ofstream& oHandle, const eb_t& disjoint_eb) {
	write_pairs (oHandle, disjoint_eb.which_pair);
	write_pairs (oHandle, disjoint_eb.dt);
	write_pairs (oHandle, disjoint_eb.qt);
	write_pairs (oHandle, disjoint_eb.cycle);
} // write_disjoint_eb

void write_pairs (std::ofstream& oHandle, const std::vector<ipair_t>& lhs) {
	int sz = lhs.size();
	oHandle << sz << "\n";
	for (int i = 0; i < sz; ++ i) {
		oHandle << lhs[i].first << "\t"	<< lhs[i].second << "\n";
	}
} // write_pairs
/* @brief	Parsing bam file, and create a list of files that stores aln
 * 			columns; this is where initila error bucket get calculated
 */
void prep_aln_file (GlobalParam& gParam, const Parameter& myPara) {

	bool is_debug_print = false;

	Ref::const_iterator it_ref = gParam.reference.begin();

	std::stringstream ss;
	std::ofstream oHandle;

	std::vector<jeb_t> jeb; // joint error bucket
	eb_t disjoint_eb;

	int debug = 0;
	/* handle each reference sequence */
	for (; it_ref != gParam.reference.end(); ++ it_ref) {

		/*if (debug == 0) {
			debug ++;
			continue;
		}*/

		int wstart = 0;

		int reflen = it_ref->second.length;

		std::cout << " Ref: " << it_ref->first << " , len = "
				<< reflen << "\n\n";

		ivec_t cov; // the vector to store coverage for all positions
		/* handle each reference sequence window by window */

		strvec_t alnfiles; // stores all files created for current ref [it_ref]

		while (wstart + 1 < reflen) {

			/*---------- create output file for current window -------- */
			int wend = std::min (wstart + myPara.windowSz - 1, reflen - 1);

			if (is_debug_print){
				std::cout << "\tprocess region " << wstart << ", " << wend << "\n";
			}

			ss.str("");
			ss << myPara.oDirNm << "/" << it_ref->first << "." << wstart
												<< "." << wend << ".region";
			std::string filename = ss.str();
			alnfiles.push_back(filename);
			if (!create_output_file (oHandle, filename)) {
				// file already exists, no creation needed
				wstart += myPara.windowSz;
				continue;
			}

			/*------ writing to file ref start, ref end ------*/
			if (wstart == 0) oHandle << 2*wstart << "\t" << 2*wend << "\n";
			else oHandle << 2*wstart - 1 << "\t" << 2*wend << "\n";

			/* ------- gather read alignment for current window ------ */
			std::vector<BamAlignment> alns;
			if (wstart == 0) {
				gather_alignments (alns, it_ref, wstart, wend, gParam, myPara);
			} else { // get the previous column so that ins can be captured
				gather_alignments (alns, it_ref, wstart - 1, wend, gParam, myPara);
			}

			if (is_debug_print){
				std::cout << "\t\t" << alns.size() << " alignments gathered\n";
			}

			/*{
				std::cout << "alns sz: " << alns.size() << "\n";
				for (int i = 0; i < alns.size(); ++ i) {
					std::cout << alns[i].Name << "\t" << alns[i].Position << "\n";
				}
			}*/

			/* read pair alignment exists */
			std::map<std::string, ipair_t> mateInfo;
			if (!gParam.alnMap.empty()) {
				get_mateInfo (mateInfo, alns, it_ref, gParam, myPara);
			}

			imap_t mate_map;
			generate_mate_map (oHandle, mate_map, mateInfo);

			/* aln_ends[i]: stores the end position of aln[i] wrt the ref
			 * [ranges]: stores the start and end positions of ref regions
			 *  covered by reads, i.e. the positions to retrieve aln cols. */
			ivec_t aln_ends;
			std::vector<ipair_t> ranges;
			get_aln_info (aln_ends, ranges, alns, gParam.alnMap);

			/* for each aln that is rvc mapped, compute original read and
			 * cigar data [not sure needed]*/
			//std::vector<rvc_aln_t> rvc_alns;
			//compute_rvc_alns (rvc_alns, alns);

			/* adjust ranges and make the min range starting at window_start */
			if (ranges.size()) {
				if (ranges[0].first < wstart) ranges[0].first = wstart;
				if (ranges[0].second < wstart) abording ("prep_aln_file: SC failed");
			}

			// print out ranges
			if (is_debug_print){
				std::cout << "\t\tCollect aln cols in " << ranges.size() << " ranges: \t";
				for (int i = 0; i < (int) ranges.size(); ++ i) {
					std::cout << "[" <<  ranges[i].first << ","	<< ranges[i].second << "]\t";
				}
				std::cout << "\n\n";
			}
			// merge read pairs if they overlap
			//merge_readpairs (alns, aln_ends, mateInfo);

			// obtain aln info for each column wrt [pos] on ref
			int cnt_debug = 0;
			for (int i = 0; i < (int) ranges.size(); ++ i) {

				std::vector<AlnColEntry> prev_col;
				// get prev_col for purpose of calculating ins b4 first of cur_col
				if (ranges[0].first > 0) {
					get_aln_column (prev_col, ranges[0].first - 1,
						myPara.ignoreBases, alns, aln_ends, gParam.qq);
				}

				for (int pos = ranges[i].first; pos <= ranges[i].second; ++ pos) {


					// given a pos on ref, retrieve aln column
					std::vector<AlnColEntry> cur_col, ins_col;

					// calculate current aln column: fIxM, MDxf, fDxM
					// No 'DD' or 'non-nt' entry recorded

					get_aln_column (cur_col, pos, myPara.ignoreBases,
									alns, aln_ends, gParam.qq) ;

					// coverage plots
					if (pos >= wstart && pos <= wend) {
						cov.push_back(cur_col.size());
					}

					// calculate inserted column b4 cur_col, if applicable
					// No DxM or MDx entry recorded
					get_ins_column (ins_col, cur_col, prev_col);

					prev_col = cur_col;

					// clean up cur_col by removing ins bases and update
					// dt for every IxM entry
					clean_aln_column (cur_col);

					// column profiling -- set .cons_type and .is_cons fields
					bool is_diverse;
					profiling (is_diverse, ins_col, true); // for ins_col

					/*
					{ // debug print out inserted column
						if (ins_col.size()) {
							std::cout << "\n insert column: " << pos << ", is_diverse = "
									<< is_diverse << "\n";
							debug_print_alnentry (ins_col);
						}
					}*/

					// update joint error probability for complete alignment columns
					// only (not for partial ones as identified by paired reads
					if (pos >= wstart && pos <= wend) {
						update_jeb (jeb, ins_col, true, false, myPara, gParam);
					}

					if (is_diverse) {
						write_col (oHandle, ins_col, true, mate_map, pos,
								myPara, gParam);
					}

					profiling (is_diverse, cur_col, false); // for cur_col

					if (pos >= wstart && pos <= wend) {
						update_jeb (jeb, cur_col, false,
								ins_col.size(), myPara, gParam);
					}

					if (is_diverse) {
						if (pos >= wstart && pos <= wend) {
							++ gParam.bonferroni;
						}
						write_col (oHandle, cur_col, false, mate_map, pos,
								myPara, gParam);
					}

				} //for (int pos = ranges[i].first
			}// for (int i = 0

			//abording ("requested exit");
			wstart += myPara.windowSz;

			oHandle.close();

			//exit(1);

		} //while (window_start + 1 < reflen)

		gParam.alnfiles.push_back(alnfiles);

		generate_covplot_Rscript (cov, myPara.oDirNm, it_ref->first);
	} // for (; it_ref

	{
		//debug_print_jeb (jeb);
	}

	/*** write to file joint error bucket ***/
	ss.clear();
	ss.str("");
	ss << myPara.oDirNm << "/" <<  gParam.reference.begin()->first << ".eb";
	gParam.ebfile = ss.str();
	if (create_output_file (oHandle, gParam.ebfile)) {
		write_jeb (oHandle, jeb, gParam.bonferroni);
		oHandle.close();
	}
	/*
	std::string disjoint_ebfile = gParam.ebfile + ".disjoint";
	if (create_output_file (oHandle, disjoint_ebfile)) {
		write_disjoint_eb (oHandle, disjoint_eb);
		oHandle.close();
	}*/

} // prep_aln_file

/* @brief	Create a coverage plot R script.
 */
void generate_covplot_Rscript (const ivec_t& cov, const std::string& path,
		const std::string& refName){

	double avgcov = 0;
	std::stringstream ss;
	std::ofstream oHandle;
	ss.str("");
	ss << path << "/" << refName << ".covplot.R";
	if (!create_output_file (oHandle, std::string (ss.str()))) {
		// file already exists, no creation needed
		warning ("can't generate coverage plot R script:" + ss.str());
		return;
	}
	int sz = cov.size();
	if (! sz) {
		oHandle << "No cov data has been collected\n";
		return;
	}

	oHandle << "x=c(";
	for (int i = 1; i < sz; ++ i) oHandle << i << ",";
	oHandle << sz << ")\n";

	oHandle << "y=c(";
	for (int i = 0; i < sz - 1; ++ i) {
		oHandle << cov[i] << ",";
		avgcov += cov[i];
	}
	oHandle << cov[sz - 1] << ")\n";
	avgcov += cov[sz - 1];

	// average coverage
	avgcov /= (1.0 * sz);
	std::cout << "\t\t(#pos, Avg Cov) = (" << sz << ", " << avgcov << ")\n";

	oHandle << "plot(x,y,col=\"red\", xaxt=\"n\", yaxt=\"n\", "
			"xlab=\"Genomic Position\",ylab=\"Coverage\","
			"type=\"l\",lwd=2,ann=T,cex.lab=0.8, ylim=c(0,max(y)+10))\n";
	oHandle << "axis(1, at=500*0:max(x))\n";
	oHandle << "axis(2, at=50*0:max(y))\n";
	// draw a horizontal line denoting avg cov
	oHandle << "abline(h=" << avgcov << ",col=4,lty=2)\n";
	oHandle.close();
} // generate_covplot_Rscript

/* @brief	This is during preparation stage.
 * 			First calculate consensus for each aln column.
 * 			1) for ins column, treat length polymorphic case only, different
 * 			Ix entries are treated differently while all non I entries are
 * 			treated the same. In other words, we count for length polymorphic
 * 			or not
 * 			2) for non-ins column, first treat length polymorphic case (a)
 * 			then treat non-length polymorphic case (b)
 * 			 	for (a), we treat different Dx entries differently, and
 * 			 	non Dx entries the same.
 * 			 	for (b) disregard Dx entries, and consider M entries only
 * 			Note: consensus is taken as alphabetically smaller one in tie
 *
 * 			Set [.cons_type], [.is_cons_norm], and [.is_cons_poly] fields
 * 			for each col entry
 */
void profiling (bool& is_diverse, std::vector<AlnColEntry>& col, bool is_ins) {

	// 1. generate profile
	std::map<std::string, int> profile, pprofile;

	for (int i = 0; i < (int) col.size(); ++ i) {
		col[i].cons_type = consensus_type (col[i]);

		if (col[i].cons_type.at(0) == 'D' || col[i].cons_type.at(0) == 'I') {
			add_to_profile (pprofile, col[i].cons_type);
		} else {
			if (is_ins)	add_to_profile (pprofile, "d");
			else {
				add_to_profile (pprofile, "i");
				add_to_profile (profile, col[i].cons_type);
			}
		}
	}

	// 2. analyze profile,
	//    set is_diverse and .is_cons_norm and .is_cons_poly fields
	is_diverse = false;

	std::string cons;
	ivec_t cnts;
	if (pprofile.size() > 1) { // analyze length polymorphic case
		analyze_profile (cons, cnts, pprofile);
		// set .is_cons_poly field
		for (int i = 0; i < (int) col.size(); ++ i) {
			if (cons.at(0) == 'I' || cons.at(0) == 'D') {
				if (col[i].cons_type.compare(cons) == 0) {
					col[i].is_cons_poly = true;
				}
			} else if (col[i].cons_type.at(0) != 'I'
						&& col[i].cons_type.at(0) != 'D') {
				col[i].is_cons_poly = true;
			}
		}

		if (cnts[cnts.size() - 2] >= 2) is_diverse = true;

	} else {
		for (int i = 0; i < (int) col.size(); ++ i) {
			col[i].is_cons_poly = true;
		}
	}

	if (is_ins) return;

	if (profile.size() > 1) { // non-length polymorphic entries
		analyze_profile (cons, cnts, profile);
		for (int i = 0; i < (int) col.size(); ++ i) {
			if (col[i].cons_type.compare(cons) == 0) {
				col[i].is_cons_norm = true;
			}
		}
		if (cnts[cnts.size() - 2] >= 2) is_diverse = true;
	} else {
		for (int i = 0; i < (int) col.size(); ++ i) {
			col[i].is_cons_norm = true;
		}
	}

} // profiling

/* @brief	Given an alignment entry, determine it's profile type for
 * 			consensus calculation.
 *
 * @method  for insertions, 				 type = I + inserted_str
 * 			for dt = M || MM || DxM: 	 type = M
 * 			for dt = MDx: 		  		 type = Dx
 * 			Note
 * 			 a) if entry.is_rv, inserted_str or M needs to be rvc
 * 			 b) ignore M if it is not nt, e.g. M = D/N
 * 			 c) x could be of multi-digit e.g. x = 11
 */
std::string consensus_type (const AlnColEntry& entry) {
	std::string type;
	if (entry.cons_type.size()) { // insertion
		type += 'I';
		if (entry.is_rv) type += xny::get_rvc_str(entry.cons_type);
		else type += entry.cons_type;
	} else { // NOT ins
		if (!entry.dt.size()) abording ("consensus_type: dt empty !");

		if (entry.dt.size() == 1) { // beginning or end of a read
			type += entry.is_rv ? xny::get_rvc_base(entry.dt.at(0))
							    : entry.dt.at(0);
		} else {
			type += entry.is_rv ? xny::get_rvc_base(entry.dt.at (1))
			 	 	 	 	 	: entry.dt.at(1);
			if (entry.dt.at(0) == 'D') { //DxM: type = M
				type = entry.is_rv ? xny::get_rvc_base(*entry.dt.rbegin())
										: *entry.dt.rbegin();
			} else { // MDx or MM
				if (entry.dt.at(1) == 'D'){ // MDx: type = Dx
					type = entry.dt.substr(1, entry.dt.size() - 1);
				}
			}
		}
	}
	return type;
} // cosensus_type

/* @brief
 */
void add_to_profile (std::map<std::string, int>& profile,
		const std::string& type){
	std::map<std::string, int>::iterator it = profile.find(type);
	if (it != profile.end()) ++ it->second;
	else profile[type] = 1;
} // add_to_profile

/* @brief	Counting strand information as well.
 * 			profile: type --> (#fwd, #rv)
 */
void add_to_profile_strand (std::map<std::string, ipair_t>& profile,
		const std::string& type, bool is_rv){
	std::map<std::string, ipair_t>::iterator it = profile.find(type);
	if (it != profile.end()) {
		if (is_rv) ++ it->second.second;
		else ++ it->second.first;
	} else {
		if (is_rv) profile[type] = ipair_t (0, 1);
		else profile[type] = ipair_t (1, 0);
	}
} // add_to_profile_strand


/* @brief	Given [profile], generate consensus [cons]
 * 			Record all profile counts in a sorted vector [cnts]
 */
void analyze_profile (std::string& cons, ivec_t& cnts,
		const std::map<std::string, int>& profile) {

	cons.clear();
	cnts.clear();

	int cons_cnt = -1;
	std::map<std::string, int>::const_iterator it_p;
	for (it_p = profile.begin(); it_p != profile.end(); ++ it_p) {
		cnts.push_back (it_p->second);
		if (it_p->second > cons_cnt) {
			cons_cnt = it_p->second;
			cons = it_p->first;
		} else if (it_p->second == cons_cnt
				&& cons.compare(it_p->first) > 0) {
			cons = it_p->first;
		}
	}
	std::sort (cnts.begin(), cnts.end());
} // analyze_profile

// ----------------------------------------------------------------------------
//			following functions are used to update error buckets
// ----------------------------------------------------------------------------

/* @brief	Update joint error bucket
 *
 */
void update_jeb (std::vector<jeb_t>& jeb, const std::vector<AlnColEntry>& col,
  bool is_ins, bool is_pre_ins, const Parameter& myPara, const GlobalParam& gParam){

	EnumParser<DiNt> parser;

	if (is_ins) { // insertion column: measure poly once
		for (int i = 0; i < (int) col.size(); ++ i) {
			int idx_jeb = calculate_eb_index (col[i].is_firstMate,
				col[i].rcycle, gParam.maxRL, parser.parseEnum(col[i].dt),
				DiNtSize, col[i].qt, myPara);

			resize_jeb (jeb, idx_jeb);

			{/*
				if(idx_jeb == 50295 || idx_jeb == 49821) {
					std::cout << idx_jeb << ":\t";
					std::cout << col[i].is_firstMate << ", " <<
							col[i].rcycle << ", " <<  col[i].dt
							<< ", " << col[i].qt << ", I\n";
				}*/
			}

			if (! col[i].is_cons_poly) { // NOT consensus
				++ jeb[idx_jeb].lp_err;
			}
			++ jeb[idx_jeb].lp_sum;
		}
	} else{ // normal column
		for (int i = 0; i < (int) col.size(); ++ i) {

			std::string dt = col[i].dt;

			// DxM or MDx, update dt to exclude 'Dx' and include flanking base
			update_D_entry_dt (dt, col[i]);

			int idx_jeb = calculate_eb_index (col[i].is_firstMate,
				col[i].rcycle, gParam.maxRL, parser.parseEnum(dt),
				DiNtSize, col[i].qt, myPara);

			/*
			{
				if(idx_jeb == 50295 || idx_jeb == 49821) {
					std::cout << idx_jeb << ":\t";
					std::cout << col[i].is_firstMate << ", " <<
					 col[i].rcycle << ", " <<  col[i].dt
					<< ", " << col[i].qt << "\n";
				}
			} */

			resize_jeb (jeb, idx_jeb);
			/*if (col[i].dt.at(0) == 'D') { // DxM: measure M as substitution once
				if (!col[i].is_cons_norm) ++ jeb[idx_jeb].snp_err;
				++ jeb[idx_jeb].snp_sum;
			} else if (col[i].dt.size() > 1 && col[i].dt.at(1) == 'D') {
				// MDx: measure as polymorphic once
				if (!col[i].is_cons_poly) ++ jeb[idx_jeb].lp_err;
				++ jeb[idx_jeb].lp_sum;
			} else { // no 'D' involved
				//if (col[i].dt.size() == 1) {
					// the first M in forward and the last M in rv: measure M as subst once
				//	if (!col[i].is_cons_norm) ++ jeb[idx_jeb].snp_err;
				//	++ jeb[idx_jeb].snp_sum;
				//} else { // MM or  the first M in forward and the last M in rv
					// exists ins column b4 *this: measure M as subst once
					if (is_pre_ins) {
						if (!col[i].is_cons_norm) ++ jeb[idx_jeb].snp_err;
						++ jeb[idx_jeb].snp_sum;
					} else {
						 // NOT exists ins b4 *this: measure once as poly, once as subst
						if (!col[i].is_cons_norm) ++ jeb[idx_jeb].snp_err;
						++ jeb[idx_jeb].snp_sum;
						if (!col[i].is_cons_poly) ++ jeb[idx_jeb].lp_err;
						++ jeb[idx_jeb].lp_sum;
					}
				//}
			}*/
			if (col[i].dt.size() > 1 && col[i].dt.at(1) == 'D') {
				// MDx: measure as polymorphic once
				if (!col[i].is_cons_poly) ++ jeb[idx_jeb].lp_err;
				++ jeb[idx_jeb].lp_sum;
			} else { // DxM, MM or the first M in forward and the last M in rv
				     // regardless of whether pre_col is ins or not
				if (!col[i].is_cons_norm) ++ jeb[idx_jeb].snp_err;
				++ jeb[idx_jeb].snp_sum;
				if (!col[i].is_cons_poly) ++ jeb[idx_jeb].lp_err;
				++ jeb[idx_jeb].lp_sum;

			}
		} // for (int i = 0
	} // normal column
} // update_jeb

/* @brief	for MDx, DxM entry, update dt to exclude 'Dx' and include
 * 			the flanking base
 */
void update_D_entry_dt (std::string& dt, const AlnColEntry& entry) {
	if (entry.dt.at(0) == 'D') { // DxM
		dt.clear();
		dt += entry.flanking_base;
		dt += *entry.dt.rbegin();
	} else if (entry.dt.size() > 1 && entry.dt.at(1) == 'D') { // MDx
		dt.clear();
		dt += entry.dt.at(0);
		dt += entry.flanking_base;
	}
} // update_D_entry_dt

//
void resize_jeb (std::vector<jeb_t>& jeb, int idx){
	int diff = idx + 1 - (int) jeb.size();
	if (diff > 0) {
		std::vector<jeb_t> tmp (diff);
		jeb.insert(jeb.end(), tmp.begin(), tmp.end());
	}
}// resize_jeb

/* @brief	Transform (isp0, c, di, q) to the errBucket index
 */
int calculate_eb_index (int isp0, int c, int maxRL, int di, int maxDI,
		int qt, const Parameter& myPara) {
	if (c > maxRL) {
		abording ("found cycle > max read length, increase parameter"
				" pSample to get the max read length");
	}
	if (!myPara.var_cycle) {
		c = 0;
		maxRL = 1;
	}
	if (!myPara.var_dt) {
		di = 0;
		maxDI = 1;
	}
	if (!myPara.var_matepair) isp0 = 0;

	return isp0*maxRL*maxDI*myPara.var_qt + c*maxDI*myPara.var_qt
			+ di*myPara.var_qt + qt;
} // cal_index

/* @brief	Update the independent error bucket entry
 */
void update_disjoint_eb_entry (eb_t& disjoint_eb, int is_p0, int cycle,
		int dt, 	int qt,	bool is_cons) {

	int incre = is_cons ? 0 : 1;

	// is_p0
	++ disjoint_eb.which_pair[is_p0].second;
	disjoint_eb.which_pair[is_p0].first += incre;

	// dt
	if (dt >= (int) disjoint_eb.dt.size()) {
		std::cout << dt << "\n";
		abording ("update_disjoint_eb_entry: SC failed, dt too large");
	}

	update_pairs (disjoint_eb.dt, dt, incre);

	update_pairs (disjoint_eb.qt, qt, incre);

	update_pairs (disjoint_eb.cycle, cycle, incre);

} // update_disjoint_eb_entry

/* @brief	Called by update_disjoint_eb_entry
 */
void update_pairs (std::vector<ipair_t>& lhs, int index, int incre) {
	int sz = index + 1 - (int) lhs.size();
	if (sz > 0) {
		std::vector<ipair_t> tmp (sz, ipair_t (0, 0));
		lhs.insert(lhs.end(), tmp.begin(), tmp.end());
	}
	lhs[index].first += incre;
	++ lhs[index].second;
} //  update_pairs

/* @brief	This is for treating err variants as independent variables.
 */
double calculate_disjoint_pe (int is_p0, int cycle,
		int dt, 	int qt, const eb_t& disjoint_eb) {

	double rslt = 1.0;

	if (is_p0 != 0 && is_p0 != 1) { // SC
		std::cout << "is_p0 = " << is_p0 << "\n";
		abording ("calculate_disjoint_pe: is_p0 is out of scope");
	}

	rslt *= (disjoint_eb.which_pair[is_p0].first + 1.0) /
			(disjoint_eb.which_pair[is_p0].second + 1.0);

	if (cycle >= (int) disjoint_eb.cycle.size()) {
		std::cout << "cycle = " << cycle << "\n";
		abording ("calculate_disjoint_pe: cycle is out of scope");
	}

	rslt *= (disjoint_eb.cycle[cycle].first + 1.0) /
			(disjoint_eb.cycle[cycle].second + 1.0);

	if (dt >= (int) disjoint_eb.dt.size()) {
		std::cout << "dt = " << dt << "\n";
		abording ("calculate_disjoint_pe: dt is out of scope");
	}

	rslt *= (disjoint_eb.dt[dt].first + 1.0) /
			(disjoint_eb.dt[dt].second + 1.0);

	if (qt >= (int) disjoint_eb.qt.size()) {
		std::cout << "qt = " << qt << "\n";
		abording ("calculate_disjoint_pe: qt is out of scope");
	}

	rslt *= (disjoint_eb.qt[qt].first + 1.0) /
			(disjoint_eb.qt[qt].second + 1.0);

	return rslt;
} //  calculate_disjoint_pe

/* @brief	if the index is larger than the current size, increase then
 * 			initialize the errBuckets. Mismatch value is incremented
 * 			by [val_mismatch] and total value is incremented by
 * 			[val_total], where these values can be negative.
 *
 */
void update_eb_entry (std::vector<ipair_t>& errBuckets, int index,
		int incre_mismatch, int incre_total) {

	int sz = index + 1 - (int) errBuckets.size();
	if (sz > 0) {
		std::vector<ipair_t> tmp (sz, ipair_t (0,0));
		errBuckets.insert(errBuckets.end(), tmp.begin(), tmp.end());
	}
	errBuckets[index].first += incre_mismatch;
	errBuckets[index].second += incre_total;
} // update_eb_entry




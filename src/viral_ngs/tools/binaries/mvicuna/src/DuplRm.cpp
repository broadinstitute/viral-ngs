//========================================================================
// Project     : M-Vicuna
// Name        : DuplRm.cpp
// Author      : Xiao Yang
// Created on  : Jun 25, 2013
// Version     : 1.0
// Copyright   : The Broad Institute
//  				 SOFTWARE COPYRIGHT NOTICE AGREEMENT
// 				 This software and its documentation are copyright (2013)
//				 by the Broad Institute. All rights are reserved.
//
// 				 This software is supplied without any warranty or 
//				 guaranteed support whatsoever. The Broad Institute cannot 
//				 be responsible for its use,	misuse, or functionality.
// Description :
//========================================================================

#include "DuplRm.h"

/** Function duplicate_removal ()
 *
 * Input: a list of paired fastq files
 * Output: a list of paired fastq files with duplicate removed. If the
 * num of output files is equal to the input, then dupl removal is applied
 * to each input pair; otherwise, 2 paired fastq output files should be
 * specified, and the output are combined in these two files.
 *
 * Randomly mutate 'N's in each read into 'A'. This will
 * result in reads containing a small number of Ns being considered for
 * clustering whereas reads containing a lot of Ns will not be considered
 * anyway.
 */
void duplicate_removal (const strvec_t& ifqs, const drm_t& drm, int w,
		int w2, xny::low_complexity& lc, int batch, bool silent) {


	// sanity check
	if ((drm.op.size() != ifqs.size()) && (drm.op.size () != 2)) {
		abording ("duplicate_removal ofqs.size() != ifqs.size() and"
				"ofqs.size() != 2");
	}
	// first obtain the length of a read, double it to be fragment length
	std::ifstream fh_tmp;
	xny::openfile<std::ifstream>(fh_tmp, ifqs[0]);
	bio::fastq_input_iterator<> fq(fh_tmp);
	int frag_len = 2 * (std::get<1>(*fq)).length();
	xny::closefile(fh_tmp);

	// calculate the upper bound of mismatches can be tolerated,
	// number of seeds & seed length; max seed length will be bounded by 31
//	int ub_mismatch = std::min (drm.max_mismatch, frag_len * (100 - drm.perc_sim) /100),
	int ub_mismatch = frag_len * (100 - drm.perc_sim) /100,
		//num_seed = ub_mismatch + 1,
		seed_len = std::min(frag_len / (ub_mismatch + 1), 31);
		//num_seed = std::max(num_seed, frag_len/seed_len);

	//std::cout << max_mismatch << ", " << frag_len << ", " << perc_sim << "\n";
	//std::cout << "ub_mismatch = " << ub_mismatch << "\n";

	// process every pair of files
	int num_file_pairs = ifqs.size()/2;

	xny::sketch_list slistgen (w, false);
	xny::super_sketch ssgen (w2);

	std::ofstream ofhfq, ofhfq2;
	if (drm.op.size() == 2) {
		xny::openfile<std::ofstream> (ofhfq, drm.op[0]);
		xny::openfile<std::ofstream> (ofhfq2, drm.op[1]);
	}

	for (int i = 0; i < num_file_pairs; ++ i) {

		int fID = 2*i;

		if (! silent) {
			std::cout << "\tprocess files: " << ifqs[fID] << " and "
					<< ifqs[fID + 1] << "\n\n";
		}

		// --------------- generate seeds for each fragment -------------
		// ----------- a compressed form to represent fragments ---------
		if (! silent) std::cout << "\tgenerate seeds...\n";
		ii64vec_t list_seeds; // stores the list of seeds per fragment,
							  // the last element stores the fragment ID
		get_seeds (list_seeds, ifqs[fID], ifqs[fID + 1], seed_len,
				batch, silent);

		// ---- initialize the global union-find structure ----
		ivec_t uf_clst (list_seeds.size());
		for (unsigned int j = 0; j < list_seeds.size(); ++ j) uf_clst[j] = j;

		// ------------------ clustering via ss -------------------------

		if (! silent) std::cout << "\tclustering via super sketches ...\n";

		clustering_via_ss (uf_clst, list_seeds, ifqs[fID], ifqs[fID + 1],
				batch, slistgen, ssgen, ub_mismatch, silent);

		// --------------- clustering via seeds -----------------
		if (!silent) std::cout << "\tclustering via seeds ...\n";
		clustering_via_seeds (uf_clst, list_seeds, frag_len/seed_len,
				ub_mismatch, silent);


		// -------- generate final union find clusters --------------
		iivec_t clusters;
		uf_generate_cls (clusters, uf_clst);

		// ---- generate the duplicated fragment IDs ----------------
		iset_t duplIDs; // duplicate fragIDs
		int debug_counter = 0;
		for (int j = 0; j < (int) clusters.size(); ++ j) {
			duplIDs.insert(clusters[j].begin() + 1, clusters[j].end());

			/*{ // debug code print out clusters
				if (clusters[j].size() > 100) {
					debug_print_fragments (clusters[j], ifqs[fID], ifqs[fID+1]);
					++ debug_counter;
					if (debug_counter > 10) exit(1);
				}
			}*/
		}
		clusters.clear();

		if (!silent) std::cout << "\n\t\tnum duplicate frags: " << duplIDs.size()
				<< "(" << 100 * duplIDs.size()/ list_seeds.size() << "% total)" << "\n\n";

		if (!silent) std::cout << "\toutput non-redundant read-pairs...\n";

		// ----- output non-redundant read-pairs ---------
		if (drm.op.size() > 2) {

			xny::openfile<std::ofstream> (ofhfq, drm.op[fID]);
			xny::openfile<std::ofstream> (ofhfq2, drm.op[fID + 1]);
		}

		clean_dupl_frag (ifqs[fID], ifqs[fID+1], ofhfq, ofhfq2,
				 duplIDs, lc, batch);

		if (drm.op.size() > 2) {
			xny::closefile(ofhfq);
			xny::closefile(ofhfq2);
		}

	} // for (int i = 0

} // duplicate_removal

/** Function debug_print_fragments ()
 * Given fragment IDs then print concatenated fragments from input fastq files
 */
void debug_print_fragments (const ivec_t& fragIDs, const std::string& fq,
		const std::string& fq2) {

	std::cout << "\nnum fragments: " << fragIDs.size() << "\n";
	iset_t ids (fragIDs.begin(), fragIDs.end());
	std::ifstream ifhfq, ifhfq2;
	xny::openfile<std::ifstream>(ifhfq, fq);
	xny::openfile<std::ifstream>(ifhfq2, fq2);
	bio::fastq_input_iterator<> iter_fq (ifhfq), end, iter_fq2(ifhfq2);
	int fragID = 0;
	for (; iter_fq != end, iter_fq2 != end; ++ iter_fq, ++ iter_fq2) {
		if (ids.count(fragID)) {
			std::string frag = std::get<1>(*iter_fq) + std::get<1> (*iter_fq2);
			std::cout << frag << "\n";
		}
		++ fragID;
	}
	xny::closefile(ifhfq);
	xny::closefile(ifhfq2);
} //debug_print_fragments

/**	Function clustering_via_ss
 *
 *	Input 1) fragments in binary representation [list_seeds]
 *		  2) super_sketches for each fragment
 */
void clustering_via_ss (ivec_t& uf_clst, const ii64vec_t& list_seeds,
	const std::string& f, const std::string& f2, int batch,
	xny::sketch_list& slistgen, xny::super_sketch& ssgen, int max_mismatch,
	bool silent) {

	int sz = list_seeds.size();
	if (sz == 0) return;

	// ------- multiple iterations of sketching ------------------

	int pre_sz = 1, iter = 0;

	while (true) {

		if (! silent) std::cout << "\n\t\tsketching iteration: " << iter << "\n";
		++ iter;

		jaz::murmur264 hashfunc (rand() % RAND_MAX);

		// -----------------------  sketching ---------------------------
		if (! silent) std::cout << "\t\t\tgenerate super sketches ...\n";
		std::vector<sketch_t> super_sketches;
		get_super_sketches (super_sketches, f, f2, slistgen, ssgen,
				hashfunc, batch, silent);
		std::sort(super_sketches.begin(), super_sketches.end(),
					xny::cmp_sketch());

		if (sz != list_seeds.size()) abording ("clusetring_via_ss SC failed.");


		// In [init_clusters], each 1d elem stores the indices of
		// [list_seeds] that share the same super sketch, where in [list_seeds],
		// the index i should be equal to fragID
		iivec_t init_clusters (1, ivec_t{super_sketches[0].second});
		for (int i = 1; i < sz; ++ i) {
			if (super_sketches[i].first == super_sketches[i-1].first) {
				init_clusters.rbegin()->push_back(super_sketches[i].second);
			} else init_clusters.push_back({super_sketches[i].second});
		}

		int init_sz = init_clusters.size();

		if (!silent){
			std::cout << "\t\t\t" << init_sz << " clusters to validate\n";
		}

		validate_clusters (uf_clst, init_clusters, list_seeds, max_mismatch,
				INT_MAX);

		// ---- generate final union find clusters ------
		iivec_t clusters;
		uf_generate_cls (clusters, uf_clst);

		// ---- generate the duplicated fragment IDs ----------------
		// iteration ending criteria % duplicate ID increase < 5%
		iset_t duplIDs;
		for (int i = 0; i < (int) clusters.size(); ++ i) {
			duplIDs.insert(clusters[i].begin() + 1, clusters[i].end());
		}
		int perc_incr = 100 * (duplIDs.size() - pre_sz)/pre_sz;
		if (!silent){
			std::cout << "\n\t\t\tduplicates: " << duplIDs.size() << " (" <<
					perc_incr << " % increase)\n" ;
		}
		if (perc_incr < 5) break;
		else pre_sz = duplIDs.size();

		break;
	}// 	while (true)

} //

/**	Function get_super_sketches
 *
 */
void get_super_sketches (std::vector<sketch_t>& super_sketches,
	const std::string& f, const std::string& f2, xny::sketch_list& slistgen,
	xny::super_sketch& ssgen, jaz::murmur264& hashfunc, int batch, bool silent) {

	std::ifstream fh, fh2;
	xny::openfile<std::ifstream>(fh, f);
	xny::openfile<std::ifstream>(fh2, f2);
	bio::fastq_input_iterator<> fq(fh), end, fq2(fh2);

	int total_read_pairs = 0;
	strvec_t pairs;

	while (fq != end && fq2 != end) {

 		add_fq_reads_only (pairs, batch/2, fq, end);
		add_fq_reads_only (pairs, batch/2, fq2, end);

		generate_super_sketches (super_sketches, pairs, slistgen,
				ssgen, hashfunc);

		total_read_pairs += pairs.size()/2;

		pairs.clear();

	} // while

	if (!silent) {
		std::cout << "\t\t\ttotal pairs, super_sketches: " << total_read_pairs
				<< ", " 	<< super_sketches.size() <<  "\n";
	}

	xny::closefile(fh);
	xny::closefile(fh2);
} // get_super_sketches

/**	Function generate_super_sketches
 *
 * Input: a list of read-pair sequences
 * Output: super sketches for each fragment
 */
void generate_super_sketches (std::vector<sketch_t>& super_sketches,
		const strvec_t& pairs, xny::sketch_list& slistgen,
		xny::super_sketch ssgen, jaz::murmur264& hashfunc) {
	bool debug = false;

	int num_frag = pairs.size()/2;
	int fragID = super_sketches.size();

	std::vector<sketch_t> batch_sketches (num_frag);

	#pragma omp parallel for
	for (int i = 0; i < num_frag; ++ i) {
		std::string frag = pairs[i] + pairs[i + num_frag];
		if (frag.size()) {
			if (std::isupper(frag.at(0))) {
				std::replace (frag.begin(), frag.end(), 'N', 'A');
			} else std::replace (frag.begin(), frag.end(), 'n', 'a');

			std::vector<sketch_t> slist = slistgen (frag, hashfunc);
			std::sort (slist.begin(), slist.end(), xny::cmp_sketch());
			ssgen (batch_sketches[i], slist, hashfunc);
			batch_sketches[i].second = fragID + i;
		}
	}

	super_sketches.insert (super_sketches.end(),
			batch_sketches.begin(), batch_sketches.end());
} //generate_super_sketches

/**	Function clustering ()
 *
 *	Given fragments, each represented by a list of seeds, represented by
 *	[list_seeds]. Cluster the fragments by pairwise comparison of ones
 *	sharing the same seed. Two fragments are included in the same cluster
 *	if they meet max_mismatch criterion. Output [fragID2clsID], where
 *	the clstID is the smallest fragID in the cluster.
 *
 *	Method: a global union find structure [uf_clust] is initialized to be
 *	the number of fragments. Clustering iterates through each seed:
 *	boundaries in [list_seeds] are identified such that each chunk of
 *	fragments share the same seed. Then clustering can be applied in parallel:
 *	within the clustering function, a local union find structure is generated
 *	and [uf_clust] is only checked but not updated for fragment comparison.
 *	Once local clusters were generated, [uf_clust] is then updated to reflect
 *	the clustering. The purpose of using this approach is to use OMP
 */
void clustering_via_seeds (ivec_t& uf_clst, ii64vec_t& list_seeds,
		int num_seed, int max_mismatch, bool silent) {
	bool debug = false;

	if (list_seeds.size() == 0 || list_seeds[0].size() == 0) {
		abording ("DuplRm.cpp -- clustering () SC failed");
	}

	int sz = list_seeds.size();

	// -------- cluster according to seed i ---------
	int num_seed_to_check = std::min(num_seed, max_mismatch + 1);
	num_seed_to_check = std::min (5, num_seed_to_check); // cap at 5 iterations
	for (int seed_i = 0; seed_i < num_seed_to_check; ++ seed_i) {

		if (!silent) {
			std::cout << "\t\tcluster by seed " << seed_i << "\n";
		}

		std::sort (list_seeds.begin(), list_seeds.end(), cmp_seed(seed_i));

		// linear scan the sorted [list_seeds] wrt the ith seed, and
		// generate 2d vector, where each dimension stores the indices of
		// [list_seeds] that share the same seed
		iivec_t init_clusters (1, ivec_t{0});
		for (int i = 1; i < sz; ++ i) {
			if (list_seeds[i][seed_i] == list_seeds[i - 1][seed_i]) {
				init_clusters.rbegin()->push_back(i);

				/*
				if (debug) { // print out substring w/ large count
					if (init_clusters.rbegin()->size() > 100000) {
						std::string fwd_str = xny::ID2Str<int64_t>
							(list_seeds[i][seed_i], 31);
						std::cout << fwd_str << "\t";
						std::cout << xny::get_rvc_str(fwd_str) << "\n";
						exit(1);
					}
				} */
			} else init_clusters.push_back({i});
		}

		int init_sz = init_clusters.size();

		if (!silent){
			std::cout << "\t\t" << init_sz << " clusters to validate\n";
		}
		//------- generate clusters: parallel clustering for each chunk
		// of boundary then merge to the global cluster -------------
		validate_clusters (uf_clst, init_clusters, list_seeds, max_mismatch,
				20000);

	} // for (int seed_i = 0; seed_i < num_seeds; ++ seed_i) {

} // clustering_via_seeds

/**	Function validate_clusters ()
 *
 */
void validate_clusters (ivec_t& uf_clst, const iivec_t& init_clusters,
	 const ii64vec_t& list_seeds, int max_mismatch, int max_cls_sz) {


	iivec_t global_clusters;
	int init_sz = init_clusters.size();

	ivec_t cls_sz (init_sz, 0);

	#pragma omp parallel
	{
		iivec_t private_clusters;
		#pragma omp for // private (private_clusters) -- do not put private here otherwise it is private to the block only
		for (int c = 0; c < init_sz; ++ c) {
			if (init_clusters[c].size() > max_cls_sz) continue;

			cls_sz[c] = init_clusters[c].size();

			make_cluster (private_clusters, list_seeds,
					init_clusters[c], max_mismatch, uf_clst);
		}
		#pragma omp critical
		{
			global_clusters.insert (global_clusters.end(),
					private_clusters.begin(), private_clusters.end());
		}
	} // #pragma omp parallel

	for (int i = 1; i < init_sz; ++ i) cls_sz[0] = std::max(cls_sz[0], cls_sz[i]);
	std::cout << "\t\t\tmax cls found: " << cls_sz[0] << "\n";

	// now update [uf_clst] according to [global_clusters]
	for (int i = 0; i < (int) global_clusters.size(); ++ i) {
		update_uf (uf_clst, global_clusters[i]);
	}
} // validate_clusters

/**	Function update_uf()
 * 	Given fragment IDs of the same cluster, update the union find
 * 	structure [uf_clust] to reflect this information
 */
void update_uf (ivec_t& uf_clst, const ivec_t& clusters) {
	int sz = clusters.size ();
	for (int i = 0; i < sz - 1; ++ i) {
		for (int j = i + 1; j < sz; ++ j) {
			int fragID_i = clusters[i],
				fragID_j = clusters[j];
			int root_i = uf_find (fragID_i, uf_clst),
				root_j = uf_find (fragID_j, uf_clst);
			uf_clst[root_j] = root_i;
		}
	}
} // update_uf

/** Function make_cluster()
 *
 *  Given a list of fragments denoted by seeds, make pairwise comparison
 *  and clustering conforming max_mismatch criteria
 *
 *  Output: clusters in 2d vector format, where each row of the vector
 *  stores the clustered fragment IDs.
 */
void make_cluster (iivec_t& clusters, const ii64vec_t& list_seeds,
		const ivec_t& init_cluster, int max_mismatch, const ivec_t& uf_clst) {

	if (list_seeds.size() == 0) {
		abording ("DuplRm.cpp -- make_cluster(): SC failed");
	}

	//--------- union find: (1) initialize the cluster ---------
	int sz = init_cluster.size();
	bvec_t visited (sz, false);
	ivec_t clst (sz);
	for (int i = 0; i < sz; ++ i) clst[i] = i;

	//---------  pairwise comparison ---------
	for (int i = 0; i < sz - 1; ++ i) {
		if (visited[i]) continue; // to speed up

		int idx_i = init_cluster[i];
		for (int j = i + 1; j < sz; ++ j) {

			if (visited[j]) continue; // to speed up, avoid of comparison
									  // if this is already clustered
			int idx_j = init_cluster[j];

			// check global uf structure according to fragID
			int root_uf_i = uf_clsfind ((int) list_seeds[idx_i].back(), uf_clst),
				root_uf_j = uf_clsfind ((int) list_seeds[idx_j].back(), uf_clst);

			if (root_uf_i != root_uf_j) {
				int root_i = uf_find (i, clst),
					root_j = uf_find (j, clst);

				if (root_i != root_j) {
					if (is_similar (list_seeds[idx_i], list_seeds[idx_j],
							max_mismatch)) {
						clst[root_j] = root_i;
						visited[j] = true;
					}
				}
			} // if
		} // for (int j = i + 1
	} // for (int i = 0

	//----- generate final cluster { clusterID --> fragment IDs } ------
	std::map<int, ivec_t> clstID_fragIDs;
	std::map<int, ivec_t>::iterator it;
	for (int i = 0; i < sz; ++ i) {
		int idx_i = init_cluster[i];
		int fragID = list_seeds[idx_i].back();
		int root_i = uf_clsfind (i, clst);
		it = clstID_fragIDs.find (root_i);
		if (it != clstID_fragIDs.end()) it->second.push_back(fragID);
		else clstID_fragIDs[root_i] = ivec_t (1, fragID);
	} // for (int i = 0

	// go through the map and produce clusters in sorted vector format
	for (it = clstID_fragIDs.begin(); it != clstID_fragIDs.end(); ++ it) {
		if (it->second.size() > 1) {
			std::sort (it->second.begin(), it->second.end());
			clusters.push_back(it->second);
		}
	} // for (it

} // make_cluster

/** Function is_similar ()
 *
 * Compare if two fragments (denoted by seeds) are similar
 */
bool is_similar (const i64vec_t& s0, const i64vec_t& s1, int max_mismatch){
	int num_seeds = s0.size() - 1;
	int mismatches = 0;
	for (int s = 0; s < num_seeds; ++ s) {
		if (s0[s] != s1[s]) {
			int diff = xny::hdlet<int64_t>(s0[s], s1[s], max_mismatch);
			mismatches += diff;
			if (diff == -1 || mismatches > max_mismatch) return false;
		}
	} // for (int s
	return true;
} // is_similar

/** Function get_seeds()
 *
 * Produce seeds: rID-> (s0, s1, s2...) for each fragment given two
 * fq files f and f2
 */
void get_seeds (ii64vec_t& list_seeds, const std::string& f,
		const std::string& f2, int seed_len, int batch, bool silent) {

	std::ifstream fh, fh2;
	xny::openfile<std::ifstream>(fh, f);
	xny::openfile<std::ifstream>(fh2, f2);
	bio::fastq_input_iterator<> fq(fh), end, fq2(fh2);

	int total_read_pairs = 0;
	strvec_t pairs;

	while (fq != end && fq2 != end) {

 		add_fq_reads_only (pairs, batch/2, fq, end);
		add_fq_reads_only (pairs, batch/2, fq2, end);

		generate_seeds (list_seeds, pairs, seed_len);

		total_read_pairs += pairs.size()/2;

		pairs.clear();

	} // while

	if (!silent) {
		std::cout << "\t\ttotal frags: " << total_read_pairs <<  "\n";
	}

	xny::closefile(fh);
	xny::closefile(fh2);
} // get_seeds

/** Function generate_seeds ()
 *
 * Given [pairs] which stores 2*n number of reads, where pairs[i] and
 * pairs[i + n] for 0 <= i <= n - 1 form a read pair
 *
 * Output: [list_seeds] the non-overlapping seeds for each read pair
 */
void generate_seeds (ii64vec_t& list_seeds, const strvec_t& pairs,
		int seed_len) {
	bool debug = false;

	int num_pairs = pairs.size()/2;
	int fragID = list_seeds.size();

	ii64vec_t batch_seeds (num_pairs);

	int debug_cnt = 0;

	#pragma omp parallel for
	for (int i = 0; i < num_pairs; ++ i) {
		std::string frag = pairs[i] + pairs[i + num_pairs];
		if (frag.size()) {
			if (std::isupper(frag.at(0))) {
				std::replace (frag.begin(), frag.end(), 'N', 'A');
			} else std::replace (frag.begin(), frag.end(), 'n', 'a');
		}
		get_seeds_per_fragment (batch_seeds[i], frag, seed_len);


		if (debug) { // print out frag and its seeds

			std::cout << frag << "\n";
			int tmp_last_seed_len = frag.length() - seed_len * (frag.length()/seed_len);
			std::cout << "fraglen, seed_len, last_seed_len = " << frag.length() << ", "
					<< seed_len << ", " << tmp_last_seed_len << "\n";
			for (auto& x : batch_seeds[i]) {
				std::string d_fwd = xny::ID2Str<int64_t>(x, seed_len);
				std::cout << d_fwd << "\t" << xny::get_rvc_str(d_fwd) << "\n";
			}
			if (tmp_last_seed_len != 0) {
				std::cout << "last str: ";
				std::string d_fwd = xny::ID2Str<int64_t>(batch_seeds[i].back(),
						tmp_last_seed_len);
				std::cout << d_fwd << "\t" << xny::get_rvc_str(d_fwd) <<  "\n";
			}
			debug_cnt ++;
			if (debug_cnt > 2) exit(1);
		}

		// append to the seed list the fragment ID
		batch_seeds[i].push_back(fragID + i);

		if (debug) { // debug print
			std::cout << frag << "\n";
			for (int j = 0; j < batch_seeds[i].size(); ++ j) {
				std::cout << batch_seeds[i][j] << "  ";
			}
			std::cout << "\n\n";
			++ debug_cnt;
			if (debug_cnt > 50) exit(1);
		}
	}

	list_seeds.insert (list_seeds.end(),
			batch_seeds.begin(), batch_seeds.end());
} // generate_seeds

/**	Function get_seeds_per_fragment ()
 *
 * For each fragment, generate non-overlapping seeds according to
 * num_seed and seed_len. Consider both forward and reverse complementary
 * for each seed and select whichever is larger
 *
 * Note: the last seed may be shorter than seed_len
 */
void get_seeds_per_fragment (i64vec_t& seeds, const std::string& frag,
		int seed_len) {

	int fraglen = frag.length();
	int num_seed = fraglen/seed_len;

	/* first go through seeds with length [seed_len] */
	seeds.resize(num_seed, -1);

	for (int i = 0; i < num_seed; ++ i) {
		if (i*seed_len + seed_len > frag.length()) {
			abording ("get_seeds_per_fragment -- SC: seed out of range");
		}
		int64_t id;
		if (xny::str2ID<int64_t> (id, frag.substr(i*seed_len, seed_len))) {
			seeds [i] = std::max(id, xny::get_rvc_bits (id, seed_len));
		}
	}

	/* possible last remaining seed with len < seed_len */
	int last_seed_len = fraglen - num_seed * seed_len;
	if (last_seed_len != 0) {
		int64_t id;
		if (xny::str2ID<int64_t> (id,
				frag.substr(fraglen - last_seed_len, last_seed_len))) {
			seeds.push_back(std::max(id,
					xny::get_rvc_bits (id, last_seed_len)));
		} else seeds.push_back(-1);
	} // if (last_seed_len != 0)

} // get_seeds_per_fragment

/**	Function clean_dupl_frag ()
 *
 */
void clean_dupl_frag (const std::string& ifq, const std::string& ifq2,
	std::ofstream& ofhfq, std::ofstream& ofhfq2, const iset_t& duplIDs,
	xny::low_complexity& lc, int batch){

	std::ifstream ifhfq, ifhfq2;
	xny::openfile<std::ifstream>(ifhfq, ifq);
	xny::openfile<std::ifstream>(ifhfq2, ifq2);
	bio::fastq_input_iterator<> fq(ifhfq), end, fq2(ifhfq2);

	int total_read_pairs = 0;
	std::vector<fqtuple_t> pairs;

	int num_lc = 0;
	int fragID = 0;
	while (fq != end && fq2 != end) {

 		add_fq_reads (pairs, batch/2, fq, end);
		add_fq_reads (pairs, batch/2, fq2, end);
		iset_t low_complex_frag;
		check_low_complexity (low_complex_frag, fragID, pairs, lc);
		num_lc += low_complex_frag.size();
		int fragnum = pairs.size()/2;
		for (int i = 0; i < fragnum; ++ i) {
			if ( (!duplIDs.count(fragID)) &&
				 (!low_complex_frag.count(fragID))) { // output
				ofhfq << "@" << std::get<0>(pairs[i]) << "\n";
				ofhfq << std::get<1>(pairs[i]) << "\n";
				ofhfq << "+\n";
				ofhfq << std::get<2>(pairs[i]) << "\n";

				ofhfq2 << "@" << std::get<0>(pairs[i + fragnum]) << "\n";
				ofhfq2 << std::get<1>(pairs[i + fragnum]) << "\n";
				ofhfq2 << "+\n";
				ofhfq2 << std::get<2>(pairs[i + fragnum]) << "\n";
			}
			++ fragID;
		}

		total_read_pairs += pairs.size()/2;

		pairs.clear();

	} // while

	std::cout << "\t\tlow complexity fragments: " << num_lc << "\n\n";
	xny::closefile(ifhfq);
	xny::closefile(ifhfq2);

} // clean_dupl_frag


void check_low_complexity (iset_t& low_complex_frag, int start_fragID,
		const std::vector<fqtuple_t>& pairs, xny::low_complexity& lc){
	int num_frag = pairs.size()/2;
	bvec_t is_lc (num_frag);
	#pragma omp parallel for
	for (int i = 0; i < num_frag; ++ i) {
		std::string frag = std::get<1> (pairs[i]) + std::get<1>(pairs[i + num_frag]);
		is_lc[i] = lc (frag);
	}

	for (int i = 0; i < num_frag; ++ i) {
		if (is_lc[i]) low_complex_frag.insert(i + start_fragID);
	}
} //check_low_complexity

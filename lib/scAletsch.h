#ifndef __GTFMERGE_H__
#define __GTFMERGE_H__

#include "genome1.h"

using namespace std;

//int do_union(const string &file, genome1 &gm);
int do_union(const vector<string> &v, genome1 *gm, const string &suffix);
int load_genome(const string &v, genome1 *gm);
int load_genome1(const string &v, genome1 *gm);
int load_genome1(const string &v, genome1 *gm); // store intron pair
int solve_one_to_multiple_pair(map<string, set<string> > &sgp1, map<string, set<string> > &mgp1, map< vector<string>, int> &sgp_n, map< vector<string>, int> &mgp_n);
int add_in_and_out_edge(transcript &t1, transcript &t2);
string empty_side(transcript &t);

class scAletsch
{
public:
	genome1 gm;

	genome1 gm1;
	genome1 gm2;
	genome1 gm3;
	
	vector< vector<int> > paths;
	vector< vector<double> > pscores;

	int cnum; // num of individual gtf
	vector< set<int> > subgraphs;
	vector< set<PI32> > sgj;
	vector<bool> sg_covered;
	
	vector<vector<double>> features; // features of the final meta-level transcripts

	vector<transcript> ft; // store all full-length transcripts

public:
	int build_union(const string &file);
	int build_union1(const string &file); // add cell_id info
	int gene_pair(const string &file);
	int gene_isoform_distribution(vector<transcript> &vref);
	int gene_distribution(vector<transcript> &vref);
	int gene_pair1(const string &file);
	int gene_pair_comb(const string &ref_file, const string &file);
    int intron_retention_filter(const string &file, const string &fo, const string &fo1);
	int gene_pair_mb(const string &file);
	int gi_distribution(const string &file);
	int split_target_transcript_biotype(const string &file, const string &target_biotype);
	int split_coverage(const string &ifile, const string &dcov, const string &lfile, const string &rfile); 
	int link_merge(const string &prefix);
	int rm_cover();
	int meta_vote(const string &prefix);
	int group_vote(const string &glist, const string &prefix);
	int split_single(const string &file, const string &fo);
	int filter_junction(const string &fi, const string &fo);
	int filter_boundary(const string &fi, const string &fo);
	int split_empty(const string &fi);
	int split_support_cell(const string &prefix);
	int split_match(const string &fi);
	int write_union(const string &fo);

	int write_individual(const string &prefix);
	int write_individual_feature(const string &prefix);
	int calculate_cell_specific_features(genome1 &ori_gm, genome1 &pre_gm, transcript &t, int cid, vector<double> & specific_features);
};

#endif

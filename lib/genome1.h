#ifndef __GENOME1_H__
#define __GENOME1_H__

#include "genome.h"

using namespace std;

typedef pair<string, int> PSI;
typedef map<string, int> MSI;

class genome1
{
public:
	genome1();
	genome1(const string &file);
	genome1(const string &file, const int cid);

public:
	vector<transcript> transcripts;
	MSI intron_hashing;
	map<int, int> cell_trst_num; //input trst num for each cell
	map<int, int> cell_gene_num; //input gene num for each cell

	map<string, vector<int> > gt; // gene_id and its transcripts list

	vector<transcript> stranscripts; // single-exon
	vector<vector<string>> mgenes; // <<gene_id, seqname>>, genes containing multi-exon transcripts
	vector< set<int32_t> > mgblist; // boundry list of genes containing multi-exon transcripts
	vector< vector< pair<int32_t, int32_t> > > mgilist; // intron list of genes containing multi-exon transcripts
	vector< vector<double> > mgilist_cov; // coverage of intron list 'mgilist'
	vector< set<string> > mgiclist; // intron-chain list of genes containing multi-exon transcripts
	map<string, int> mgi; // <gene_id, index> index in mgenes
	vector<transcript> irtranscripts; // intron-retention transcripts
	vector<transcript> irftranscripts; // intron-retention filtered transcripts

	map< string, set<int> > cjlist; // cell list of junction
	map< string, set<int> > celist; // cell list of exon
	set< pair<int32_t, int32_t> > ilist;
	map< string, set<int> > tjlist; // transcript list of junction
	map< string, set<int> > telist; // transcript list of exon 	
	vector< vector<int> > subgraph;

	map<string, vector<int32_t>> jblist;

	map< string, set<int> > jtocnum; // cell list of each junction
public:
	int add_transcript(const transcript &t);
	int add_transcript_b(const transcript &t);
	int add_transcript_simple(const transcript &t);
	int build(const string &file);
	int build(const vector<transcript> &v);
	int build_union(const genome1 &gm);
	int build_intersection(const genome1 &gm, genome1 &out);
	int clear();
	int write(const string &file);
	int add_suffix(const string &p);
	int remove_redundancy();
	int print(int index);
	int print_hashing();

	int build1(const string &file);
	string compute_intron_hashing(const transcript &t);
	int build2(const string &file);
	int build3(const string &file);
	int write_individual(const vector<string> &ogtf);
	
	int build_cid(const string &file, int cid);
	int add_cjlist(const transcript &t);
	int add_celist(const transcript &t);
	int add_ilist(const transcript &t);
	int add_tjlist(const transcript &t, int tidx);

private:
	int build_multiexon_transcripts(const string &file);
};

string tostring(int p);
//string compute_intron_hashing(const transcript &t);

#endif

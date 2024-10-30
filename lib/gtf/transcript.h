/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __GTF_TRANSCRIPT_H__
#define __GTF_TRANSCRIPT_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include "item.h"

using namespace std;

typedef pair<int32_t, int32_t> PI32;

class transcript
{
public:
	transcript(const item &ie);
	transcript();
	~transcript();

public:
	bool operator< (const transcript &t) const;

public:
	string seqname;
	string source;
	string feature;
	string gene_id;
	string transcript_id;
	string gene_type;
	string transcript_type;
	int32_t start;
	int32_t end;
	double score;
	char strand;
	int frame;
	double coverage;
	double covratio;
	double RPKM;
	double FPKM;
	double TPM;

	vector<PI32> exons;

	vector< vector<int> > nfp;
	set<int> si;
	set<int> so;
	map<int, int> siw;
	map<int, int> sow;
	set<int> clist; // cell id list
	map<int, int> cjcount; // cjcout: pair of cell id and the number of overlapping junctions in this cell
	map<int, int> cw; // cell weight: sum of support junction; e.g. <1, 2>, <2, 3> -> <<1,1>,<2,2>,<3,1>>
	map<int, int> cwe; // cell weight: sum of support exon

	string path;
	vector<int> fragPath;
	vector<double> scoreList;
	vector< vector<double> > empty_vote; // < <left_non_empty, left_empty>, <right_non_empty, right_empty> >
	bool left_empty;
	bool right_empty;
	vector<double> empty_style; // coverage of {both nonempty, both empty, left empty only, right empty only}
	double hcov; // highest cov in all t
	vector<transcript> tlist;
	int32_t sb; // boundary interval of start
	int32_t eb; // boundary interval of end

public:
	int add_exon(int s, int t);
	int add_exon(const item &e);
	int assign_RPKM(double factor);
	int sort();
	int clear();
	int shrink();
	int assign(const item &e);
	int length() const;
	PI32 get_bounds() const;
	PI32 get_first_intron() const;
	vector<PI32> get_intron_chain() const;
	bool intron_chain_match(const transcript &t) const;
	string label() const;
	int write(ostream &fout) const;
};

#endif

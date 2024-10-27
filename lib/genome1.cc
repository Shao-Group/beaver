#include "genome1.h"
#include "config.h"
#include <cassert>
#include <fstream>
#include <cmath>

genome1::genome1()
{}

genome1::genome1(const string &file)
{
	build(file);
}

genome1::genome1(const string &file, const int cid)
{
	build_cid(file, cid);
}

int genome1::build1(const string &file)
{
	clear();
	genome gm(file);
	int tcount = 0;
	for(int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];
		tcount += g.transcripts.size();
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			transcript t = g.transcripts[k];
			if(t.gene_id == "") continue;
			if(t.exons.size() <= 1) continue;
			if(merge_coverage_as_counts == true) t.coverage = 1.0;
			else if(merge_coverage_log == true) t.coverage = log(1 + t.coverage);
			add_transcript(t);
		}
	}
	//printf("total # gene: %d\ntotal # transcript: %d\n", gm.genes.size(), tcount);
	for(int k = 0; k < transcripts.size(); k++)
	{
		transcript t = transcripts[k];
		string gid = t.gene_id;
		gt[gid].push_back(k);
	}

	/*
	// print name of all genes
	set<string> gene_list;

	for(int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];

                for(int k = 0; k < g.transcripts.size(); k++)
                {
			transcript t = g.transcripts[k];
			if(t.gene_id == "") continue;
			gene_list.insert(t.gene_id);
		}
	}

	printf("total # gene: %d\n", gene_list.size());
	set<string>::iterator it;
	for(it = gene_list.begin(); it != gene_list.end(); it++)
	{
		printf("%s\n", (*it).c_str());
	}
	*/
	

        return 0;
}

int genome1::build2(const string &file)
{
        clear();
        genome gm(file);
	int tcount = 0;
	bool add;
	set<int32_t> blist;
	string gid;
	string seqname;

	// add all transcripts
        for(int i = 0; i < gm.genes.size(); i++)
        {
                const gene &g = gm.genes[i];
                for(int k = 0; k < g.transcripts.size(); k++)
                {
                        transcript t = g.transcripts[k];
                        if(t.exons.size() <= 1) continue;
                        if(merge_coverage_as_counts == true) t.coverage = 1.0;
                        else if(merge_coverage_log == true) t.coverage = log(1 + t.coverage);
                        add_transcript(t);
                }
        }

        for(int i = 0; i < gm.genes.size(); i++)
        {
                const gene &g = gm.genes[i];
		tcount += g.transcripts.size();
		blist.clear();
		add = false;
                for(int k = 0; k < g.transcripts.size(); k++)
                {
                        transcript t = g.transcripts[k];
                        if(t.exons.size() <= 1) {
				stranscripts.push_back(t);
				continue;
			}
			else {
				if(add == false) {
					add = true;
					gid = t.gene_id;
					seqname = t.seqname;
				}
				blist.insert(t.exons[0].second);
				for(int j = 1; j < t.exons.size(); j++)
				{
					blist.insert(t.exons[j].first);
					blist.insert(t.exons[j].second);
				}
			}
                }
		if(add == true)
		{
			vector<string> mg;
			mg.push_back(gid);
			mg.push_back(seqname);
			mgenes.push_back(mg);
			mgblist.push_back(blist);
		}
        }
	printf("total # genes = %ld, # single-exon genes = %ld, # multi-exon genes = %ld\ntotal # transcripts = %d, # single-exon transcripts = %ld, # multi-exon transcripts = %ld\n", gm.genes.size(), gm.genes.size() - mgenes.size(), mgenes.size(), tcount, stranscripts.size(), tcount - stranscripts.size());
	return 0;

}

int genome1::build3(const string &file)
{
        clear();
        genome gm(file);
        int tcount = 0;
	int icount = 0;
        bool add;
        vector< pair<int32_t, int32_t> > ilist;
	vector<double> ilist_cov;
	set<int32_t> blist;
	set<string> iclist;

        string gid;
        string seqname;

	// add all transcripts
        for(int i = 0; i < gm.genes.size(); i++)
        {
                const gene &g = gm.genes[i];
                for(int k = 0; k < g.transcripts.size(); k++)
                {
					transcript t = g.transcripts[k];
					//if(t.exons.size() <= 1) continue;
					if(t.exons.size() < 1) continue;
					if(merge_coverage_as_counts == true) t.coverage = 1.0;
					else if(merge_coverage_log == true) t.coverage = log(1 + t.coverage);
					add_transcript(t);
                }
        }

	int total_b = 0;
	int total_i = 0;
	int total_ic = 0;

	// add boundary list: blist, add intron list: ilist, add intron-chain list: iclist
        for(int i = 0; i < gm.genes.size(); i++)
        {
			const gene &g = gm.genes[i];
			tcount += g.transcripts.size();
			ilist.clear();
			ilist_cov.clear();
			blist.clear();
			iclist.clear();
                add = false;
                for(int k = 0; k < g.transcripts.size(); k++)
                {
                        transcript t = g.transcripts[k];
                        if(t.exons.size() > 1) {
                                if(add == false) {
                                        add = true;
                                        gid = t.gene_id;
                                        seqname = t.seqname;
                                }
                                for(int j = 0; j < t.exons.size() - 1; j++)
                                {
					ilist.push_back({t.exons[j].second, t.exons[j+1].first});
					ilist_cov.push_back(t.coverage);
					blist.insert(t.exons[j].second);
					blist.insert(t.exons[j+1].first);
					icount++;
                                }
				iclist.insert(compute_intron_hashing(t));
                        }
                }
                if(add == true)
                {
                        vector<string> mg;
                        mg.push_back(gid);
                        mg.push_back(seqname);
                        mgenes.push_back(mg);
                        mgilist.push_back(ilist);
			mgilist_cov.push_back(ilist_cov);
			mgblist.push_back(blist);
			mgiclist.push_back(iclist);
			mgi[gid] = mgenes.size() - 1;

			total_b += blist.size();
			total_i += ilist.size();
			total_ic += iclist.size();
                }
        }
		
	printf("total # genes = %ld, # multi-exon genes = %ld, # total boundaries = %d, # total junctions = %d, # total intron-chains = %d \n", gm.genes.size(), mgenes.size(), total_b, total_i, total_ic);
        return 0;
}

int genome1::build_cid(const string &file, int cid)
{
        clear();
        genome gm(file);
		int total_trst_num = 0;

        for(int i = 0; i < gm.genes.size(); i++)
        {
			const gene &g = gm.genes[i];
			for(int k = 0; k < g.transcripts.size(); k++)
			{
				transcript t = g.transcripts[k];
							
				if(t.exons.size() <= 1) continue;
				//if(t.exons.size() < 1) continue;
				if(t.coverage < min_input_transcript_coverage) continue;
				total_trst_num++;

				if(merge_coverage_as_counts == true) t.coverage = 1.0;
				else if(merge_coverage_log == true) t.coverage = log(1 + t.coverage);

				t.clist.insert(cid);
				t.tlist.push_back(t);

				/*	
				t.empty_style = {0,0,0,0};
				if(t.nfp.size() == 0) t.empty_style[0] += t.coverage;
				if(t.nfp.size() == 2) t.empty_style[1] += t.coverage;
				if(t.nfp.size() == 1) 
				{
					if(t.nfp[0][0] <= t.start + 1) t.empty_style[2] += t.coverage;
					if(t.nfp[0][1] == t.end) t.empty_style[3] += t.coverage;
				}

				t.empty_vote = {{0,0},{0,0}};
				if(t.nfp.size() == 2)
				{
					t.empty_vote[0][1] += t.coverage;
					t.empty_vote[1][1] += t.coverage;
				}
				else if(t.nfp.size() == 1)
				{
					if(t.nfp[0][0] <= t.start + 1) t.empty_vote[0][1] += t.coverage;
					if(t.nfp[0][1] == t.end) t.empty_vote[1][1] += t.coverage;
				}
				
				else if(t.nfp.size() == 0)
				{
					t.empty_vote[0][0] += t.coverage;
					t.empty_vote[1][0] += t.coverage;
				}
				*/

				

				add_transcript(t);

				add_cjlist(t);
				//add_celist(t);

				//add_tjlist(t, k);
				//add_ilist(t);
			}
        }

		cell_trst_num[cid] = total_trst_num;
		cell_gene_num[cid] = gm.genes.size();
		printf("Loading cell %d: %d transcripts in %d genes.\n", cid, cell_trst_num[cid], cell_gene_num[cid]);

		return 0;
}

int genome1::build(const string &file)
{
	clear();
	genome gm(file);
	for(int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			transcript t = g.transcripts[k];
			if(t.exons.size() <= 1) continue;
			if(merge_coverage_as_counts == true) t.coverage = 1.0;
			else if(merge_coverage_log == true) t.coverage = log(1 + t.coverage);
			add_transcript(t);
		}
	}
	return 0;
}

int genome1::build(const vector<transcript> &v)
{
	clear();
	for(int i = 0; i < v.size(); i++)
	{
		add_transcript(v[i]);
	}
	return 0;
}

int genome1::clear()
{
	transcripts.clear();
	intron_hashing.clear();
	gt.clear();
	stranscripts.clear();
	mgenes.clear();
	mgblist.clear();
	mgilist.clear();
	mgilist_cov.clear();
	mgiclist.clear();
	mgi.clear();
	irtranscripts.clear();
	irftranscripts.clear();

	cjlist.clear();
	celist.clear();
	ilist.clear();

	subgraph.clear();

	return 0;
}

int genome1::add_ilist(const transcript &t)
{
	if(t.exons.size() <= 1) return 0;
	for(int i = 0; i < t.exons.size() - 1; i++)
        {
		ilist.insert({t.exons[i].second, t.exons[i+1].first});
	}

	return 0;
}

int genome1::add_tjlist(const transcript &t, int tidx)
{
	if(t.exons.size() <= 1) return 0;

	for(int i = 0; i < t.exons.size() -1; i++)
	{
		string h = t.seqname;
                h.append(tostring(t.exons[i].second));
                h.append(tostring(t.exons[i+1].first));

                tjlist[h].insert(tidx);
	}
	return 0;
}

int genome1::add_cjlist(const transcript &t)
{
	if(t.exons.size() <= 1) return 0;
	
	for(int i = 0; i < t.exons.size() - 1; i++)
	{
		string h = t.seqname;
		h.append(tostring(t.exons[i].second));
		h.append(tostring(t.exons[i+1].first));
		
		for(auto c : t.clist)
		{
			cjlist[h].insert(c);
		}
	}

	return 0;
}

int genome1::add_celist(const transcript &t)
{
        if(t.exons.size() <= 1) return 0;

        for(int i = 0; i < t.exons.size(); i++)
        {
                string h = t.seqname;
                h.append(tostring(t.exons[i].first));
                h.append(tostring(t.exons[i].second));

                for(auto c : t.clist)
                {
                    celist[h].insert(c);
                }
        }

        return 0;
}

int genome1::add_transcript(const transcript &t)
{
        string s = compute_intron_hashing(t);
        if(intron_hashing.find(s) == intron_hashing.end())
        {
			intron_hashing.insert(PSI(s, transcripts.size()));
			transcripts.push_back(t);
        }
        else
        {
			int k = intron_hashing[s];
			assert(k >= 0 && k < transcripts.size());
                
			if(transcripts[k].hcov >= t.hcov)
			//if(transcripts[k].length() >= t.length())
			{
				transcripts[k].coverage += t.coverage;
				for(auto c : t.clist) transcripts[k].clist.insert(c);
				for(auto c : t.tlist) transcripts[k].tlist.push_back(c);
				//for(int i = 0; i < 4; i++) transcripts[k].empty_style[i] += t.empty_style[i];

				/*
				if(t.empty_vote.size() != 0)
				{	
					transcripts[k].empty_vote[0][0] += t.empty_vote[0][0];
					transcripts[k].empty_vote[0][1] += t.empty_vote[0][1];
					transcripts[k].empty_vote[1][0] += t.empty_vote[1][0];
					transcripts[k].empty_vote[1][1] += t.empty_vote[1][1];
				}
				*/
			}
			else
			{
				transcripts[k].hcov = t.hcov;
				double c = transcripts[k].coverage + t.coverage;
				set<int> pre_clist = transcripts[k].clist;
				vector<transcript> pre_tlist = transcripts[k].tlist;
				//vector<double> pre_es =  transcripts[k].empty_style;
				//vector<vector<double> > pre_ev = transcripts[k].empty_vote;
                        
				transcripts[k] = t;
				transcripts[k].coverage = c;
				for(auto c : pre_clist) transcripts[k].clist.insert(c);
				for(auto c : pre_tlist) transcripts[k].tlist.push_back(c);
				//for(int i = 0; i < 4; i++) transcripts[k].empty_style[i] += pre_es[i];

				/*
				if(t.empty_vote.size() != 0)
				{
					transcripts[k].empty_vote[0][0] += pre_ev[0][0];
					transcripts[k].empty_vote[0][1] += pre_ev[0][1];
					transcripts[k].empty_vote[1][0] += pre_ev[1][0];
					transcripts[k].empty_vote[1][1] += pre_ev[1][1];
				}
				*/
		}

        }
        return 0;
}

int genome1::add_transcript_simple(const transcript &t)
{
	transcripts.push_back(t);
	return 0;
}

int genome1::add_transcript_b(const transcript &t)
{
	string ss = compute_intron_hashing(t);
	string s = to_string(t.sb) + ss + to_string(t.eb);

	if(intron_hashing.find(s) == intron_hashing.end())
	{
		intron_hashing.insert(PSI(s, transcripts.size()));
		transcripts.push_back(t);
	}
	else
	{
		int k = intron_hashing[s];
		assert(k >= 0 && k < transcripts.size());

		if(transcripts[k].hcov >= t.hcov)
		{
			transcripts[k].coverage += t.coverage;
			for(auto c : t.clist) transcripts[k].clist.insert(c);
			for(auto c : t.tlist) transcripts[k].tlist.push_back(c);
		}
		else
		{
			transcripts[k].hcov = t.hcov;
			double c = transcripts[k].coverage + t.coverage;
			set<int> pre_clist = transcripts[k].clist;
			vector<transcript> pre_tlist = transcripts[k].tlist;
			//vector<double> pre_es =  transcripts[k].empty_style;
			//vector<vector<double> > pre_ev = transcripts[k].empty_vote;

			transcripts[k] = t;
			transcripts[k].coverage = c;
			for(auto c : pre_clist) transcripts[k].clist.insert(c);
			for(auto c : pre_tlist) transcripts[k].tlist.push_back(c);
			//for(int i = 0; i < 4; i++) transcripts[k].empty_style[i] += pre_es[i];

			/*
			if(t.empty_vote.size() != 0)
			{
				transcripts[k].empty_vote[0][0] += pre_ev[0][0];
				transcripts[k].empty_vote[0][1] += pre_ev[0][1];
				transcripts[k].empty_vote[1][0] += pre_ev[1][0];
				transcripts[k].empty_vote[1][1] += pre_ev[1][1];
			}
			*/
		}

    }
	return 0;
}

int genome1::build_intersection(const genome1 &gm, genome1 &out)
{
	out.clear();
	for(MSI::iterator it = intron_hashing.begin(); it != intron_hashing.end(); it++)
	{
		string s = it->first;
		int k1 = it->second;
		transcript t = transcripts[k1];
		MSI::const_iterator x = gm.intron_hashing.find(s);
		if(x == gm.intron_hashing.end()) continue;
		int k2 = x->second;
		t.coverage += gm.transcripts[k2].coverage;
		out.add_transcript(t);
	}
	return 0;
}

int genome1::build_union(const genome1 &gm)
{
	// add transcripts
	for(int k = 0; k < gm.transcripts.size(); k++)
	{
		const transcript &t = gm.transcripts[k];
		add_transcript(t);
		add_cjlist(t);
	}

	// update cell raw information
	for(auto p : gm.cell_trst_num)
	{
		int cid = p.first;
		int num = p.second;

		if(cell_trst_num.find(cid) != cell_trst_num.end())
		{
			printf("Error: Duplicate cell files unioned.\n");
		}

		cell_trst_num.insert(p);
	}
	for(auto p : gm.cell_gene_num)
	{
		int cid = p.first;
		int num = p.second;

		if(cell_gene_num.find(cid) != cell_gene_num.end())
		{
			printf("Error: Duplicate cell files unioned.\n");
		}

		cell_gene_num.insert(p);
	}

	return 0;
}

int genome1::add_suffix(const string &p)
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		t.transcript_id.append("-").append(p);
		t.gene_id.append("-").append(p);
	}
	return 0;
}

int genome1::print(int index)
{
	printf("genome %d: %lu transcripts, %lu distinct first intron\n", index, transcripts.size(), intron_hashing.size());
	return 0;
}

int genome1::print_hashing()
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		string s = compute_intron_hashing(transcripts[i]);
		printf("hash = %s\n", s.c_str());
	}
	return 0;
}

int genome1::write(const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail()) return 0;
	
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		if(t.coverage < min_transcript_coverage) continue;
		t.write(fout);
	}
	fout.close();
	return 0;
}

int genome1::write_individual(const vector<string> &ogtf)
{
	for(int i = 0; i < ogtf.size(); i++)
	{
		ofstream fout(ogtf[i].c_str());
		if(fout.fail()) return 0;
		//printf("writing idividual gtf # %d...\n", i+1);

		for(int j = 0; j < transcripts.size(); j++)
		{
			transcript &t = transcripts[j];
			if(t.coverage < min_transcript_coverage) continue;

			if(t.clist.find(i+1) == t.clist.end()) continue;
			
			t.write(fout);
		}
		
		fout.close();
	}
	return 0;
}

string tostring(int p)
{
	char buf[10240];
	sprintf(buf, "%d", p);
	return string(buf);
}

string genome1::compute_intron_hashing(const transcript &t)
{
	string h = t.seqname;
	
	// do not consider strand
	/*
	if(t.strand == '.') h.append("0");
	if(t.strand == '+') h.append("1");
	if(t.strand == '-') h.append("2");
	*/

	// single exon, seqname + start + (end - start)
	if(t.exons.size() <= 1) 
	{
		h.append(tostring(t.start));
		h.append(tostring(t.end - t.start));
		return h;
	}

	// multi exon
	int32_t p = t.exons[0].second;
	h.append(tostring(p));

	for(int k = 1; k < t.exons.size(); k++)
	{
		int32_t q1 = t.exons[k].first;
		int32_t q2 = t.exons[k].second;
		h.append(tostring(q1 - p));
		if(k == t.exons.size() - 1) break;
		h.append(tostring(q2 - q1));
		p = q2;
	}
	return h;
}

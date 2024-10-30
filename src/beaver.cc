#include "genome1.h"
#include "beaver.h"
#include "config.h"
#include <string>
#include <thread>
#include <queue>
#include <algorithm>
#include <cmath>

int do_union(const vector<string> &v, genome1 *gm, const string &suffix)
{
	gm->clear();
	int cnt = 0;	
	for(int k = 0; k < v.size(); k++)
	{
		string s = v[k];
		genome1 g(s);
		g.add_suffix(suffix + "-" + tostring(k));
		gm->build_union(g);
		cout << "union genome " << ++cnt << " " << s.c_str() << endl << flush;
	}
	return 0;
}

int do_union1(const vector< pair<string, int> > &v, genome1 *gm, const string &suffix)
{
        gm->clear();
        int cnt = 0;
        for(int k = 0; k < v.size(); k++)
        {
			string s = v[k].first;
			int cid = v[k].second;
			genome1 g(s, cid);
			
			g.add_suffix(suffix + "-" + tostring(k));
			gm->build_union(g);
			cout << "union genome " << ++cnt << " " << s.c_str() << endl << flush;
        }
        return 0;
}

int load_genome(const string &s, genome1 *gm)
{
	gm->clear();
	gm->build1(s);
	return 0;
}

int load_genome1(const string &s, genome1 *gm)
{
        gm->clear();
        gm->build2(s);
        return 0;
}

int load_genome2(const string &s, genome1 *gm)
{
        gm->clear();
        gm->build3(s);
        return 0;
}

int load_genome_all(const string &s, genome1 *gm, int cid)
{
	gm->clear();
	gm->build_cid(s, cid);
	return 0;
}


int solve_one_to_multiple_pair(map<string, set<string> > &sgp1, map<string, set<string> > &mgp1, map< vector<string>, int> &sgp_n, map< vector<string>, int> &mgp_n)
{
	printf("solve one to multiple...\n");
	set<string> glist1;
	set<string> glist2;
	for(map<string, set<string> >::iterator it = mgp1.begin(); it != mgp1.end(); it++)
	{
		if(it->second.size() > 1)
		{
			string g1 = it->first;
			glist1.insert(g1);
			for(set<string>::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++)
			{
				string g2 = *it1;
				glist2.insert(g2);
				if(it->second.size() > 9)
					cout << it->first << ":" << g2 << "# matching = " << mgp_n[{g1,g2}] << endl;
			}
		}
	}
	printf("# g1 = %ld, # g2 = %ld\n", glist1.size(), glist2.size());
	return 0;
}

int Beaver::build_union1(const string &file)
{
	vector<genome1> gv(num_threads);

	vector< vector< pair<string, int> > > sv1(num_threads);
	
	ifstream fin(file.c_str());
	if(fin.fail()) return 0;
	string line;
	int cnt = 0;
	while(getline(fin, line))
	{
		if(line == "") continue;
		int k = cnt % num_threads;
		sv1[k].push_back({line, cnt+1});
		cnt++;
	}
	fin.close();
	cnum = cnt;

	vector<thread> threads;
	for(int k = 0; k < sv1.size(); k++)
	{
		printf("thread %d processes %lu genomes\n", k, sv1[k].size());
		threads.emplace_back(do_union1, sv1[k], &(gv[k]), string("x").append(tostring(k)));
	}
	for(int k = 0; k < threads.size(); k++)
	{
		threads[k].join();
	}

	gm.clear();
	for(int k = 0; k < gv.size(); k++)
	{
		printf("final union genome %d\n", k);
		gm.build_union(gv[k]);
	}
	
	return 0;
}

int Beaver::build_union(const string &file)
{
	vector<genome1> gv(num_threads);
	vector< vector<string> > sv(num_threads);

	ifstream fin(file.c_str());
	if(fin.fail()) return 0;
	string line;
	int cnt = 0;
	while(getline(fin, line))
	{
		if(line == "") continue;
		int k = cnt % num_threads;
		sv[k].push_back(line);
		cnt++;
	}
	fin.close();

	vector<thread> threads;
	for(int k = 0; k < sv.size(); k++)
	{
		printf("thread %d processes %lu genomes\n", k, sv[k].size());
		threads.emplace_back(do_union, sv[k], &(gv[k]), string("x").append(tostring(k)));
	}
	for(int k = 0; k < threads.size(); k++)
	{
		threads[k].join();
	}

	gm.clear();
	for(int k = 0; k < gv.size(); k++)
	{
		printf("final union genome %d\n", k);
		gm.build_union(gv[k]);
	}

	return 0;
}

int Beaver::gene_pair(const string &file)
{
	vector<string> sv;
	gm1.clear();
	gm2.clear();

	ifstream fin(file.c_str());
        if(fin.fail()) return 0;
        string line;
        int cnt = 0;
        while(getline(fin, line))
        {
                if(line == "") continue;
                sv.push_back(line);
                cnt++;
        }
        fin.close();

	load_genome(sv[0], &gm1);
	printf("gm1: # gene = %ld, # transcript = %ld\n", gm1.gt.size(), gm1.transcripts.size());
	load_genome(sv[1], &gm2);
	printf("gm2: # gene = %ld, # transcript = %ld\n", gm2.gt.size(), gm2.transcripts.size());

	map<vector<string>, vector< vector<string> >> gp;// map<gene_pair, vector<transcript_pair>>
	vector<transcript> tlist;
	tlist.clear();
	for(int i = 0; i < gm1.transcripts.size(); i++)
	{
		transcript t = gm1.transcripts[i];
		string s = gm1.compute_intron_hashing(t);
		if(gm2.intron_hashing.find(s) != gm2.intron_hashing.end())
		{
			vector<string> curg;
			curg.push_back(t.gene_id);
			int idx = gm2.intron_hashing[s];
			curg.push_back(gm2.transcripts[idx].gene_id);
			tlist.push_back(t);

			vector<string> curt;
			curt.push_back(t.transcript_id);
			curt.push_back(gm2.transcripts[idx].transcript_id);
			gp[curg].push_back(curt);
		}
	}	

	printf("total # gene pair = %ld\n", gp.size());
	int count = 0;
	int tc = 0;
	for(map<vector<string>, vector< vector<string> >>::iterator it = gp.begin(); it != gp.end(); it++)
	{
		count++;
		tc += it->second.size();
		//printf("#%d g1 = %s, g2 = %s, # identical = %d, # t1 = %d, # t2 = %d\n", count, it->first[0].c_str(), it->first[1].c_str(), it->second.size(), gm1.gt[it->first[0]].size(), gm2.gt[it->first[1]].size());
		//printf("%s %s %d %d %d ", it->first[0].c_str(), it->first[1].c_str(), it->second.size(), gm1.gt[it->first[0]].size(), gm2.gt[it->first[1]].size());
		printf("%s, %s\n", it->first[0].c_str(), it->first[1].c_str());
		/*
		vector<vector<string> > curtp = it->second;
		for(int k = 0; k < curtp.size(); k++)
		{
			if(k != curtp.size()-1) printf("%s:%s,", curtp[k][0].c_str(), curtp[k][1].c_str());
			else printf("%s:%s", curtp[k][0].c_str(), curtp[k][1].c_str());
		}
		printf("\n");
		*/
	}
	printf("total # transcript pair = %d\n", tc);
	

	//gene_distribution(tlist);
	//gene_isoform_distribution(tlist);

	return 0;
}

int Beaver::gene_distribution(vector<transcript> &vref)
{
        map<string, int> glist; // <gid, # transcripts>
        map<string, int> gtype; // <gene_biotype, # genes in this biotype>
        for(int k = 0; k < vref.size(); k++)
        {
                glist[vref[k].gene_id]++;
                if(glist[vref[k].gene_id] == 1) gtype[vref[k].gene_type]++;
        }
        printf("total number of Genes: %ld\n", glist.size());
        for(map<string, int>::iterator it1 = gtype.begin(); it1 != gtype.end(); it1++) {
                cout << it1->first << " " << it1->second << "\n";
        }
        return 0;
}

int Beaver::gene_isoform_distribution(vector<transcript> &vref)
{
	printf("total number of transcripts: %ld\n", vref.size());
        map< string, map<string, int> > gcount;// <gene_biotype, <transcript_biotype, # transcripts>>
        map<string, int> tcount; // <transcript_biotype, # transcripts>
        for(int k = 0; k < vref.size(); k++)
        {
                string gbio = vref[k].gene_type;
                string tbio = vref[k].transcript_type;
                gcount[gbio][tbio]++;
                tcount[tbio]++;
        }
        for(map<string, map<string, int> >::iterator it1 = gcount.begin(); it1 != gcount.end(); it1++) {
                for(map<string, int>::iterator it2 = tcount.begin(); it2 != tcount.end(); it2++) {
                        string tbio = it2->first;
                        if(!it1->second[tbio]) it1->second[tbio] = 0;
                }
        }
        for(map<string, int>::iterator it1 = tcount.begin(); it1 != tcount.end(); it1++) {
                printf("%s ", it1->first.c_str());
        }
        printf("\n");
        for(map<string, map<string, int> >::iterator it1 = gcount.begin(); it1 != gcount.end(); it1++) {
                printf("%s ", it1->first.c_str());
                for(map<string, int>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                        printf("%d ", it2->second);
                }
                printf("\n");
        }
        return 0;
}

int Beaver::gene_pair1(const string &file)
{
        vector<string> sv;
        gm1.clear();
        gm2.clear();

        ifstream fin(file.c_str());
        if(fin.fail()) return 0;
        string line;
        int cnt = 0;
        while(getline(fin, line))
        {
                if(line == "") continue;
                sv.push_back(line);
                cnt++;
        }
        fin.close();

	printf("loading gm1...\n");
        load_genome1(sv[0], &gm1);
	printf("loading gm2...\n");
        load_genome1(sv[1], &gm2);

	set< vector<string> > sgp; // single-exon gene pairs
	set< vector<string> > mgp; // multi-exon gene pairs
	map< vector<string>, int> sgp_n; // single-exon gene pairs, <<gene pair>, match boundary number>
	map< vector<string>, int> mgp_n; // multi-exon gene pairs, <<gene pair>, match boundary number>

	int len, olen, x1, y1, x2, y2 = 0;
	for(int i = 0; i < gm1.stranscripts.size(); i++)
	{
		x1 = gm1.stranscripts[i].start;
		y1 = gm1.stranscripts[i].end;
		for(int j = 0; j < gm2.stranscripts.size(); j++)
		{
			if(gm1.stranscripts[i].seqname != gm2.stranscripts[j].seqname) continue;
			olen = 0;
			len = y1 - x1;
			x2 = gm2.stranscripts[j].start;
                	y2 = gm2.stranscripts[j].end;
			if(y2 - x2 < len) len = y2 - x2;
			if(x2 >= x1 && x2 <= y1) {
				if(y1 > y2) olen = y2 - x2;
				else olen = y1 - x2;
			}
			else if(x1 >= x2 && x1 <= y2) {
				if(y1 > y2) olen = y2 - x1;
				else olen = y1 - x1;
			}
			if(olen > 0.8 * len) 
			{
				vector<string> cur_gp = {gm1.stranscripts[i].gene_id, gm2.stranscripts[j].gene_id};
				sgp.insert(cur_gp);
				if(sgp_n.find(cur_gp) != sgp_n.end()) sgp_n[cur_gp]++;
				else sgp_n[cur_gp] = 1;
			}
		}
	}

	set<int32_t> l1, l2;
	for(int i = 0; i < gm1.mgenes.size(); i++) {
		l1 = gm1.mgblist[i];
		for(int j = 0; j < gm2.mgenes.size(); j++) {
			if(gm1.mgenes[i][1] != gm2.mgenes[j][1]) continue;
			l2 = gm2.mgblist[j];
			for(set<int32_t>::iterator it = l1.begin(); it != l1.end(); it++) {
				if(l2.find(*it) != l2.end()) {
					vector<string> cur_gp = {gm1.mgenes[i][0], gm2.mgenes[j][0]};
					mgp.insert(cur_gp);
					if(mgp_n.find(cur_gp) != mgp_n.end()) mgp_n[cur_gp]++;
					else mgp_n[cur_gp] = 1;
					//break;
				}
			}
		}
	}

	printf("one-one: # total = %ld, # single-exon gene pairs = %ld, # multi-exon gene pairs = %ld\n", sgp.size()+mgp.size(), sgp.size(), mgp.size());

	map<string, set<string> > sgp1;
	map<string, set<string> > mgp1;
	for(set< vector<string> >::iterator it = sgp.begin(); it != sgp.end(); it++) {
		vector<string> gp_id = *it;
		sgp1[gp_id[0]].insert(gp_id[1]);
	}
	for(set< vector<string> >::iterator it = mgp.begin(); it != mgp.end(); it++) {
                vector<string> gp_id = *it;
                mgp1[gp_id[0]].insert(gp_id[1]);
        }
	printf("one-multiple: # total = %ld, # single-exon gene pairs = %ld, # multi-exon gene pairs = %ld\n", sgp1.size()+mgp1.size(), sgp1.size(), mgp1.size());

	vector<int> sgp1_distribution(101,0);
	vector<int> mgp1_distribution(101,0);
	for(map<string, set<string> >::iterator it = sgp1.begin(); it != sgp1.end(); it++) {
		sgp1_distribution[it->second.size()]++;
	}
	
	for(map<string, set<string> >::iterator it = mgp1.begin(); it != mgp1.end(); it++) {
                mgp1_distribution[it->second.size()]++;
		/*
		if(it->second.size() >= 9) 
		{
			printf("%s\n", it->first.c_str());
			set<string> cur = it->second;
			for(set<string>::iterator it1 = cur.begin(); it1 != cur.end(); it1++){
				string cc = *it1;
				printf("%s\n", cc.c_str());
			}
		}
		*/
        }
	
	printf("distribution of sgp: \n");
	int count = 0;
	for(int i = 1; i < 101; i++) {
		if(sgp1_distribution[i] != 0) printf("%d = %d\n", i, sgp1_distribution[i]);
		count += sgp1_distribution[i];
	}
	printf("total = %d\n", count);
	printf("distribution of mgp: \n");
	count = 0;
        for(int i = 1; i < 101; i++) {
		if(mgp1_distribution[i] != 0) printf("%d = %d\n", i, mgp1_distribution[i]);
		count += mgp1_distribution[i];
	}
	printf("total = %d\n", count);

	//solve_one_to_multiple_pair(sgp1, mgp1, sgp_n, mgp_n);

	
	// print all gene pairs
	printf("list of single-exon gene pairs:\n");
	for(set< vector<string> >::iterator it = sgp.begin(); it != sgp.end(); it++) {
		vector<string> cgp = *it;
		printf("%s, %s\n", cgp[0].c_str(), cgp[1].c_str());
	}
	printf("list of multi-exon gene pairs:\n");
        for(set< vector<string> >::iterator it = mgp.begin(); it != mgp.end(); it++) {
                vector<string> cgp = *it;
                printf("%s, %s\n", cgp[0].c_str(), cgp[1].c_str());
        }
	

	return 0;
}

int Beaver::gene_pair_comb(const string &ref_file, const string &file)
{
        vector<string> sv;
        gm1.clear();
        gm2.clear();

        ifstream fin(file.c_str());
        if(fin.fail()) return 0;
        string line;
        int cnt = 0;
        while(getline(fin, line))
        {
                if(line == "") continue;
                sv.push_back(line);
                cnt++;
        }
        fin.close();

        load_genome(sv[0], &gm1);
	//printf("gm1: # gene = %d, # transcript = %d\n", gm1.gt.size(), gm1.transcripts.size());
	//map<string, vector<int> >::iterator it;
	//for(it = gm1.gt.begin(); it != gm1.gt.end(); it++)
	//{
		//printf("%s\n", it->first.c_str());
	//}
        //load_genome(sv[1], &gm2);
	//printf("gm2: # gene = %d, # transcript = %d\n", gm2.gt.size(), gm2.transcripts.size());
	//for(it = gm2.gt.begin(); it != gm2.gt.end(); it++)
        //{
                //printf("%s\n", it->first.c_str());
        //}
	return 0;
}

int Beaver::intron_retention_filter(const string &file, const string &fo, const string &fo1)
{
        load_genome2(file.c_str(), &gm1);
	bool if_ir = false;
	int mcount = 0;
	printf("# total transcripts = %ld, ", gm1.transcripts.size());
	for(int i = 0; i < gm1.transcripts.size(); i++)
        {
		if_ir = false;
		transcript t = gm1.transcripts[i];
		if(t.exons.size() <= 1) 
		{
			gm1.irftranscripts.push_back(t);
			continue;
		}
		mcount++;

		double t_cov = t.coverage;

		int mgindex = gm1.mgi[t.gene_id];
		vector< pair<int32_t, int32_t> > ilist = gm1.mgilist[mgindex];
		vector<double> ilist_cov = gm1.mgilist_cov[mgindex];

		double cov_threshold = 2.0 ;

		for(int j = 0; j < t.exons.size(); j++)
		{
			int32_t lpos = t.exons[j].first;
			int32_t rpos = t.exons[j].second;
			
			
			// check if first/last exon included in an intron
			if(j == 0)
                        {
                                for(int k = 0; k < ilist.size() - 1; k++)
                                {
					//if(ilist[k].second - ilist[k].first < 80) continue;
					if(t_cov > cov_threshold * ilist_cov[k]) continue;
                                        if(lpos > ilist[k].first && lpos < ilist[k].second && rpos == ilist[k+1].first)
                                        {
						//if((rpos - lpos) / (ilist[k].second - ilist[k].first) > 4.0) continue;
                                                gm1.irtranscripts.push_back(t);
                                                if_ir = true;
                                                break;
                                        }
                                }
                        }
                        else if(j == t.exons.size() - 1)
                        {
                                for(int k = 1; k < ilist.size(); k++)
                                {
					//if(ilist[k].second - ilist[k].first < 80) continue;
                                        if(t_cov > cov_threshold * ilist_cov[k]) continue;
					if(rpos > ilist[k].first && rpos < ilist[k].second && lpos == ilist[k-1].second)
                                        {
						//if((rpos - lpos) / (ilist[k].second - ilist[k].first) > 4.0) continue;
                                                gm1.irtranscripts.push_back(t);
                                                if_ir = true;
                                                break;
                                        }
                                }
                        }
			if(if_ir == true) break;

			// exon overlap whole intron
			for(int k = 0; k < ilist.size(); k++)
			{
				//if(ilist[k].second - ilist[k].first < 80) continue;
				//if(ilist[k].second - ilist[k].first > 1500) continue;
				if(t_cov > cov_threshold * ilist_cov[k]) continue;
				if(ilist[k].first > lpos && ilist[k].second < rpos)
				{
					// skip if exon length >> intron length
					if((rpos - lpos) / (ilist[k].second - ilist[k].first) > 4.0) continue;
					gm1.irtranscripts.push_back(t);
					if_ir = true;
					break;
				}
			}
			if(if_ir == true) break;
		}
		if(if_ir == false) gm1.irftranscripts.push_back(t);
	}
	printf("# multi-exon transcripts = %d, # intron-retention transcripts = %ld, # transcripts after intron-retention filtering = %ld\n", mcount, gm1.irtranscripts.size(), gm1.irftranscripts.size());

	ofstream fout(fo.c_str());
	for(int i = 0; i < gm1.irtranscripts.size(); i++)
	{
		transcript t = gm1.irtranscripts[i];
		t.write(fout);
	}

	ofstream fout1(fo1.c_str());
	for(int i = 0; i < gm1.irftranscripts.size(); i++)
	{
		transcript t = gm1.irftranscripts[i];
		t.write(fout1);
	}

	return 0;
}

int Beaver::gene_pair_mb(const string &file)
{
        vector<string> sv;
        gm1.clear();
        gm2.clear();

        ifstream fin(file.c_str());
        if(fin.fail()) return 0;
        string line;
        int cnt = 0;
        while(getline(fin, line))
        {
                if(line == "") continue;
                sv.push_back(line);
                cnt++;
        }
        fin.close();

        printf("loading gm1...\n");
        load_genome2(sv[0], &gm1);
        printf("loading gm2...\n");
        load_genome2(sv[1], &gm2);

	set< vector<string> > mgp; // multi-exon gene pairs

	// find multi-exon gene pairs using boundary matching
	set<int32_t> l1, l2;
	for(int i = 0; i < gm1.mgenes.size(); i++)
	{
		l1 = gm1.mgblist[i];
                for(int j = 0; j < gm2.mgenes.size(); j++) {
                        if(gm1.mgenes[i][1] != gm2.mgenes[j][1]) continue;
                        l2 = gm2.mgblist[j];
                        for(set<int32_t>::iterator it = l1.begin(); it != l1.end(); it++) {
                                if(l2.find(*it) != l2.end()) {
                                        vector<string> cur_gp = {gm1.mgenes[i][0], gm2.mgenes[j][0]};
                                        mgp.insert(cur_gp);
					break;
                                }
                        }
                }
	}


	printf("# multi-exon gene pairs = %ld\n", mgp.size());

	// count 
	set< pair<int, string> > shared_ic;
	set< pair<int, pair<int32_t, int32_t> > > shared_j;
	set< pair<int, int32_t> > shared_b;

	for(set< vector<string> >::iterator it = mgp.begin(); it != mgp.end(); it++) {
                vector<string> cgp = *it;
		int idx1 = gm1.mgi[cgp[0]];
		int idx2 = gm2.mgi[cgp[1]];

		// count shared intron-chains
		set<string> ic1 = gm1.mgiclist[idx1];
		set<string> ic2 = gm2.mgiclist[idx2];
		int ic_count = 0;
		for(set<string>::iterator it_ic = ic1.begin(); it_ic != ic1.end(); it_ic++){
			if(ic2.find(*it_ic) != ic2.end()) {
				ic_count++;
				shared_ic.insert({idx1, *it_ic});
			}
		}
		
		// count shared junctions
		set< pair<int32_t, int32_t> > j1(gm1.mgilist[idx1].begin(), gm1.mgilist[idx1].end());
		set< pair<int32_t, int32_t> > j2(gm2.mgilist[idx2].begin(), gm2.mgilist[idx2].end());
		//set< pair<int32_t, int32_t> > j2(i2.begin(), i2.end());
		int j_count = 0;

		for(set< pair<int32_t, int32_t> >::iterator it_j = j1.begin(); it_j != j1.end(); it_j++){
			if(j2.find(*it_j) != j2.end()) {
				j_count++;
				shared_j.insert({idx1, *it_j});
			}
		}
		
		// count shared boundaries
		set<int32_t> b1 = gm1.mgblist[idx1];
		set<int32_t> b2 = gm2.mgblist[idx2];
		int b_count = 0;

		for(set<int32_t>::iterator it_b = b1.begin(); it_b != b1.end(); it_b++){
			if(b2.find(*it_b) != b2.end()) {
				b_count++;
				shared_b.insert({idx1, *it_b});
			}
		}

                printf("%s, %s, %d, %ld, %ld, %d, %ld, %ld, %d, %ld, %ld\n", cgp[0].c_str(), cgp[1].c_str(), ic_count, ic1.size(), ic2.size(), j_count, j1.size(), j2.size(), b_count, b1.size(), b2.size());
        }

	printf("# total shared boundaries = %ld, # total shared junctions = %ld, # total shared intron-chains = %ld\n", shared_b.size(), shared_j.size(), shared_ic.size());
	return 0;
}

int Beaver::gi_distribution(const string &file)
{	
	load_genome1(file.c_str(), &gm1);

	printf("single-exon transcript distribution:\n");
	gene_isoform_distribution(gm1.stranscripts);
	
	printf("multi-exon transcript distribution:\n");
	gene_isoform_distribution(gm1.transcripts);

	return 0;
	}

int Beaver::split_target_transcript_biotype(const string &file, const string &target_biotype)
{
        load_genome1(file.c_str(), &gm1);

	ofstream fout(target_biotype + "-" + file);
        for(int i = 0; i < gm1.transcripts.size(); i++)
        {
                transcript t = gm1.transcripts[i];
		if(t.transcript_type == target_biotype.c_str()) t.write(fout);
        }
        return 0;

}

int Beaver::split_coverage(const string &ifile, const string &dcov, const string &lfile, const string &rfile)
{
	load_genome_all(ifile.c_str(), &gm, 0);
	
	double cov_threshold = stod(dcov);

	ofstream fout1(lfile);
	ofstream fout2(rfile);
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		if(gm.transcripts[i].coverage < cov_threshold) gm.transcripts[i].write(fout1);
		else gm.transcripts[i].write(fout2);
	}
	return 0;
}

int Beaver::split_single(const string &file, const string &fo)
{
	load_genome_all(file.c_str(), &gm, 0);
	
	ofstream fout(fo.c_str());
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		if(gm.transcripts[i].exons.size() == 1) gm.transcripts[i].write(fout);
	}

	fout.close();

	return 0;
}

int remove_ir(genome1 &gm, genome1 &gm1)
{
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript &t = gm.transcripts[i];
		bool if_ir = false;
		set<int> t_clist = t.clist;

		for(int j = 0; j < t.exons.size(); j++)
		{
			for(auto intron : gm.ilist)
			{
				if(intron.first > t.exons[j].first && intron.second < t.exons[j].second)
				{
					string h = t.seqname;
					h.append(tostring(intron.first));
					h.append(tostring(intron.second));
					
					//if(gm.cjlist[h].size() < t_clist.size()) continue;

					if_ir = true;
					break;
				}
			}
			
			if(if_ir == true) break;
		}

		if(if_ir == false) gm1.add_transcript(t);
	}

	return 0;
}

int store_ft(genome1 &gm, vector<transcript> &ft)
{
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript &t = gm.transcripts[i];
		if(t.nfp.size() == 0) 
		{
			//printf("transcript #%d is full-length...\n",i);
			ft.push_back(t);
		}
	}

	return 0;
}

int addback_ft(genome1 &gm, vector<transcript> &ft)
{
	for(int i = 0; i < ft.size(); i++)
	{
		gm.add_transcript(ft[i]);
	}
	return 0;
}

int add_transcript_cw(genome1 &gm)
{
	map < string, set<int> > &cjlist = gm.cjlist;

	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript &t = gm.transcripts[i];
		if(t.exons.size() <= 1) continue;

		map<int, int> temp_cw;
		
		for(int i = 0; i < t.exons.size() - 1; i++)
		{
			string h = t.seqname;
			h.append(tostring(t.exons[i].second));
			h.append(tostring(t.exons[i+1].first));
			
			set<int> cl = cjlist[h];
			for(int c : cl)
			{
				if(temp_cw.find(c) == temp_cw.end()) temp_cw[c] = 1;
				else temp_cw[c]++;
			}
		}

		t.cw = temp_cw;
	}

	return 0;
}

int add_transcript_cwe(genome1 &gm)
{
        map < string, set<int> > &celist = gm.celist;

        for(int i = 0; i < gm.transcripts.size(); i++)
        {
                transcript &t = gm.transcripts[i];
                if(t.exons.size() <= 1) continue;

                map<int, int> temp_cwe;

                for(int i = 0; i < t.exons.size(); i++)
                {
                        string h = t.seqname;
                        h.append(tostring(t.exons[i].first));
                        h.append(tostring(t.exons[i].second));

                        set<int> cl = celist[h];
                        for(int c : cl)
                        {
                                if(temp_cwe.find(c) == temp_cwe.end()) temp_cwe[c] = 1;
                                else temp_cwe[c]++;
                        }
                }

                t.cwe = temp_cwe;
        }

        return 0;
}

int32_t find_interval(int32_t &p, vector<int32_t> &L)
{
        int32_t before = -1;

        int low = 0;
        int high = L.size() - 1;

        while(low <= high)
        {
                int mid = low + (high - low) / 2;

                if (L[mid] == p) {
                        before = p;
                        break;
                } else if (L[mid] < p) {
                        before = L[mid];
                        low = mid + 1;
                } else {
                        high = mid - 1;
                }
        }

        return before;
}

int get_edge_weight(transcript &t1, transcript &t2, int same_exon_num)
{
	int w = 0;
	
	map<int, int> cw1 = t1.cw;
	map<int, int> cw2 = t2.cw;

	for(auto i : cw1)
	{
		if(cw2.find(i.first) != cw2.end())
		{
			w += cw2[i.first];
		}
	}

	map<int, int> cwe1 = t1.cwe;
    map<int, int> cwe2 = t2.cwe;

	for(auto i : cwe1)
    {
		if(cwe2.find(i.first) != cwe2.end())
		{
				w += cwe2[i.first];
		}
	}

	//w *= same_exon_num + 1;

	return w;
}

string empty_side(transcript &t)
{
	
	string ep = "undecide"; // empty position, left, right, both
        
	
	if(t.nfp.size() == 1)
	{
			if(t.nfp[0][0] <= t.start + 1) ep = "left";
			else if(t.nfp[0][1] == t.end) ep = "right";
			else ep = "middle";
	}
	else if(t.nfp.size() == 2) ep = "both";
	
	/*	
	if(t.left_empty == true && t.right_empty == true) ep = "both";
	else if(t.left_empty == true && t.right_empty == false) ep = "left";
	else if(t.left_empty == false && t.right_empty == true) ep = "right";
	*/

	return ep;
}

int share_empty_info(genome1 &gm)
{
	map<string, set<int32_t> > left_nonempty;
	map<string, set<int32_t> > right_nonempty;

	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript &t = gm.transcripts[i];
		string es = empty_side(t);

		if(es == "left")
		{
			right_nonempty[t.seqname].insert(t.exons[t.exons.size()-1].first);
			//left_empty[t.seqname].insert(t.exons[0].second);
		}
		else if(es == "right")
		{
			left_nonempty[t.seqname].insert(t.exons[0].second);
			//right_empty[t.seqname].insert(t.exons[t.exons.size()-1].first);
		}
		/*
		else if(es == "both")
		{
			left_empty[t.seqname].insert(t.exons[0].second);
			right_empty[t.seqname].insert(t.exons[t.exons.size()-1].first);
		}
		*/
	}

	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript &t = gm.transcripts[i];
		string es = empty_side(t);

		if(es == "left")
		{
			if(left_nonempty[t.seqname].find(t.exons[0].second) != left_nonempty[t.seqname].end())
			{
                        	t.nfp.clear();
			}
		}
		if(es == "right")
		{
			if(right_nonempty[t.seqname].find(t.exons[t.exons.size()-1].first) != right_nonempty[t.seqname].end())
			{
                        	t.nfp.clear();
			}
		}
		if(es == "both")
		{
			if(right_nonempty[t.seqname].find(t.exons[t.exons.size()-1].first) != right_nonempty[t.seqname].end()) t.nfp.pop_back();
			if(left_nonempty[t.seqname].find(t.exons[0].second) != left_nonempty[t.seqname].end()) t.nfp.pop_back();
		}

	}

	return 0;
}

vector< vector<int> > share_empty_left(transcript t1, transcript t2)
{
	vector< vector<int> > temp = t1.nfp;

	// if t2 is incompte in left, t1 also incompte in left
	if(empty_side(t2) == "both" || empty_side(t2) == "left")
	{
		if(empty_side(t1) != "left" && empty_side(t1) != "both")
		{
			vector<int> target;
			target.push_back((int)t1.exons[0].first);
			target.push_back((int)t1.exons[0].second);

			temp.insert(temp.begin(), target);
		}
	}
	return temp;
}

vector< vector<int> > share_empty_right(transcript t1, transcript t2)
{
        vector< vector<int> > temp = t2.nfp;

        // if t1 is incompte in right, t2 also incompte in right
        if(empty_side(t1) == "both" || empty_side(t1) == "right")
        {
                if(empty_side(t2) != "right" && empty_side(t2) != "both")
                {
			vector<int> target;
                        target.push_back(t2.exons[t2.exons.size()-1].first);
                        target.push_back(t2.exons[t2.exons.size()-1].second);
                        temp.push_back(target);
                }
        }
        return temp;
}

bool if_pos_similar(int32_t a, int32_t b)
{
	int shift = max_shift_bp;
	if(a == b) return true;
	if(a > b && a - b <= shift) return true;
	if(b > a && b - a <= shift) return true;
	
	return false;
}

int print_transcript(transcript &t, int idx)
{
	string ep = empty_side(t);

	cout<<"transcript "<<idx<<": chr "<<t.seqname<<", tid "<<t.transcript_id<<", ("<<t.start<<","<<t.end<<"), "<<"empty-side "<<ep<<endl;
	printf("# cell = %ld, cell list: ", t.clist.size());
	for(auto it = t.clist.begin(); it != t.clist.end(); ++it)
	{
		cout<<*it<<", ";
	}
	printf("\n");

	printf("# exons = %ld, exon list: ", t.exons.size());
	for(int i = 0; i < t.exons.size(); i++)
        {
                printf("(%d,%d),",t.exons[i].first, t.exons[i].second);
        }
        printf("\n");

	printf("# si = %ld, si list: ", t.si.size());
	set<int>::iterator itr; 
	for(itr = t.si.begin(); itr != t.si.end(); itr++)
	{
		cout << "(" << *itr << ", " << t.siw[*itr] << "), ";
	}
	printf("\n");

	printf("# so = %ld, so list: ", t.so.size());
        for(itr = t.so.begin(); itr != t.so.end(); itr++)
        {
                cout << "(" << *itr << ", " << t.sow[*itr] << "), ";
        }
        printf("\n");

	return 0;
}

int print_transcript_pair(transcript &t1, transcript &t2)
{
	string ep1 = empty_side(t1);
        string ep2 = empty_side(t2);

        cout<<"t1 "<<t1.transcript_id<<" ("<<t1.start<<","<<t1.end<<")"<<" "<<ep1<<" "<<t1.nfp[0][0]<<","<<t1.nfp[0][1]<<endl;
        cout<<"t2 "<<t2.transcript_id<<" ("<<t2.start<<","<<t2.end<<")"<<" "<<ep2<<" "<<t2.nfp[0][0]<<","<<t2.nfp[0][1]<<endl;

        for(int i = 0; i < t1.exons.size(); i++)
        {
                printf("(%d,%d),",t1.exons[i].first, t1.exons[i].second);
        }
        printf("\n");
        for(int i = 0; i < t2.exons.size(); i++)
        {
                printf("(%d,%d),",t2.exons[i].first, t2.exons[i].second);
        }
        printf("\n");
	
	return 0;
}

int print_path(int idx, vector<int> path, genome1 &gm)
{
	printf("path %d: # nodes = %ld, (", idx, path.size());
	for(int i = 0; i < path.size(); i++)
	{
		printf("%d, ", path[i]);
	}
	printf(")\n");
	for(int i = 0; i < path.size(); i++)
	{
		print_transcript(gm.transcripts[path[i]], path[i]);
	}
	return 0;
}

int is_covered(transcript &t1, transcript &t2)
{
	bool if_covered = false;
	
	if(t1.exons[0].second < t2.exons[0].second || t1.exons.back().first > t2.exons.back().first) return 0;

	vector< pair<int32_t, int32_t> > exons1 = t1.exons;
	vector< pair<int32_t, int32_t> > exons2 = t2.exons;
	int idx1 = 0;
	int idx2 = 0;

	int start_overlap_idx = 0;

	while(idx2 < exons2.size())
	{
		if(exons2[idx2].second < exons1[idx1].first)
		{
				idx2++;
				continue;
		}

		start_overlap_idx = idx2;

		// find overlap start, t1 need to be exactly part of t2
		while(idx1 < exons1.size() && idx2 < exons2.size())
		{
			// if t1 first exon, only need to match in the end
			if(idx1 == 0)
			{
				if(if_pos_similar(exons1[idx1].second, exons2[idx2].second))
				//if(exons1[idx1].second == exons2[idx2].second)
				{
						idx1++;
						idx2++;
						continue;
				}
				else break;
			}
			// if t1 last exon
			if(idx1 == exons1.size() - 1)
			{
				if(if_pos_similar(exons1[idx1].first, exons2[idx2].first))
				//if(exons1[idx1].first == exons2[idx2].first)
				{
						//printf("covered!\n");
						if_covered = true;
				}
				break;
			}

			if(if_pos_similar(exons1[idx1].first, exons2[idx2].first) && if_pos_similar(exons1[idx1].second, exons2[idx2].second))
			//if(exons1[idx1].first == exons2[idx2].first && exons1[idx1].second == exons2[idx2].second)
			{
					idx1++;
					idx2++;
					continue;
			}
			else break;
         }
        break;
    }

	//int max_diff_bp = 10;
        if(if_covered == true) 
	{
		/*
		if(idx2 != exons2.size()-1)
		{
			if(t1.end >= exons2[idx2].second + max_diff_bp) return 0;
		}	
		else
		{
			if(t1.start < exons2[start_overlap_idx].first - max_diff_bp) return 0;
		}
		*/

		return 1;
	}
	
        return 0;
}

bool is_extendable(transcript &t1, transcript &t2)
{
	bool if_extendable = false;

	// t1 need to be extenable to t2
	if(t1.exons[0].second >= t2.exons[0].second || t1.exons.back().first >= t2.exons.back().first) return 0;

	vector< pair<int32_t, int32_t> > exons1 = t1.exons;
        vector< pair<int32_t, int32_t> > exons2 = t2.exons;
        int idx1 = 0;
        int idx2 = 0;

	int start_overlap_idx = 0;

        while(idx2 < exons2.size())
        {
		if(exons2[idx2].second < exons1[idx1].first)
                {
                        idx2++;
                        continue;
                }

		start_overlap_idx = idx2;

                // find overlap start, t1 need to be exactly part of t2
                while(idx1 < exons1.size() && idx2 < exons2.size())
                {
                        // if t1 first exon, only need to match in the end
                        if(idx1 == 0)
                        {
                                if(if_pos_similar(exons1[idx1].second, exons2[idx2].second))
                                //if(exons1[idx1].second == exons2[idx2].second)
                                {
                                        idx1++;
                                        idx2++;
                                        continue;
                                }
                                else break;
                        }
			// if t1 last exon
                        if(idx1 == exons1.size() - 1)
                        {
                                if(if_pos_similar(exons1[idx1].first, exons2[idx2].first))
                                //if(exons1[idx1].first == exons2[idx2].first)
                                {
                                        if_extendable = true;
                                }
                                break;
                        }

                        if(if_pos_similar(exons1[idx1].first, exons2[idx2].first) && if_pos_similar(exons1[idx1].second, exons2[idx2].second))
                        //if(exons1[idx1].first == exons2[idx2].first && exons1[idx1].second == exons2[idx2].second)
                        {
                                idx1++;
                                idx2++;
                                continue;
                        }
                        else break;
                }
                break;

	}

	return if_extendable;
}

int add_in_and_out_edge(transcript &t1, transcript &t2, int idx1, int idx2)
{	
	int s1 = t1.start;
	int s2 = t2.start;
	int e1 = t1.end;
	int e2 = t2.end;
	vector< pair<int32_t, int32_t> > exons1 = t1.exons;
	vector< pair<int32_t, int32_t> > exons2 = t2.exons;

	//bool if_add_edge = false;
	int same_exon_num = 0;

	// check overlap
	int i = 0;
	while(i < exons1.size())
	{
		// move next until find first overlap
		if(exons1[i].second < exons2[0].first) 
		{
			i++;
			continue;
		}

		// find overlap, then check if overlap parts are the same
		int j = 0;
		same_exon_num = 0;
		bool if_match = true;
		int num_unmatch_exon = 0;

		//print_transcript_pair(t1, t2);
		//printf("overlap start exon pair: exon #%d (%d, %d), exon #%d (%d, %d)\n", i, exons1[i].first, exons1[i].second,j, exons2[j].first, exons2[j].second);

		while(i < exons1.size() && j < exons2.size())
		{
			// i is the last exon, and j is the first exon
			if(i == exons1.size() - 1 && j == 0)
			{
				if(exons2[j].first >= exons1[i].first && exons2[j].first <= exons1[i].second && exons2[j].second >= exons1[i].second) break;
				else
				{
					num_unmatch_exon++;

					if(num_unmatch_exon >= max_unmatch_exon) 
					{
						if_match = false;
						break;
					}
					//break;
				}
			}

			// j is the first exon, only need to match in the end
			if(j == 0) 
			{
				if(if_pos_similar(exons1[i].second, exons2[j].second))
				{
					i++;
					j++;
					same_exon_num++;
					continue;
				}
				else 
				{
					num_unmatch_exon++;

					if(num_unmatch_exon >= max_unmatch_exon) 
					{
						if_match = false;
						break;
					}
					//break;
				}
			}

			// last exon pair
			if(j == exons2.size() - 1 || exons2[j].second > t1.end || i == exons1.size() - 1)
			//if(j == exons2.size() - 1 || exons2[j].second > t1.end)
			{
				if(if_pos_similar(exons1[i].first, exons2[j].first))
                                {
                                        i++;
                                        j++;
                                        same_exon_num++;
					continue;
                                }
                                else 
				{
					num_unmatch_exon++;

					if(num_unmatch_exon >= max_unmatch_exon) 
					{
						if_match = false;
						break;
					}
					//break;
				}
			}

			if(if_pos_similar(exons1[i].first,exons2[j].first) && if_pos_similar(exons1[i].second, exons2[j].second))
			{
				i++;
				j++;
				same_exon_num++;
			}
			else
			{
				num_unmatch_exon++;
				if(num_unmatch_exon > max_unmatch_exon) 
				{
					if_match = false;
					break;
				}
				//break;
			}
		}

		if(if_match == true) 
		{
			//print_transcript_pair(t1, t2);
			
			// only add in/out edge if can extend
			if(same_exon_num < t1.exons.size() && same_exon_num < t2.exons.size())
			{
				int w = get_edge_weight(t1, t2, same_exon_num);
				t1.so.insert(idx2);
				t1.sow[idx2] = w;
				t2.si.insert(idx1);
				t2.siw[idx1] = w;

				//if_add_edge = true;
				return 1;
			}
				
			if(is_covered(t1, t2))
			{
				// t1 is covered by t2 and in head, then add edge t1->t2
				if(t1.exons[0].second == t2.exons[0].second)
				{
					int w = get_edge_weight(t1, t2, same_exon_num);
					t1.so.insert(idx2);
					t1.sow[idx2] = w;
					t2.si.insert(idx1);
					t2.siw[idx1] = w;

					//t1.nfp.clear();
					//t1.nfp = share_empty_left(t1, t2);
					return 1;

				}
			}
			if(is_covered(t2, t1))
			{
				// if t2 is covered by t1 and in tail, then add edge t1 -> t2
				if(t1.exons[t1.exons.size()-1].first == t2.exons[t2.exons.size()-1].first)
				{
					int w = get_edge_weight(t1, t2, same_exon_num);
					t1.so.insert(idx2);
					t1.sow[idx2] = w;
					t2.si.insert(idx1);
					t2.siw[idx1] = w;

					//t2.nfp.clear();
					//t2.nfp = share_empty_right(t1, t2);
					return 1;
				}
			}
			
			
			

		}
		/*	
		if(is_covered(t2, t1) || is_covered(t1, t2))
		{	
			// t1 is covered by t2 and in head, then add edge t1->t2
			if(t1.exons[0].second == t2.exons[0].second)
			{
				int w = get_edge_weight(t1, t2, same_exon_num);
				t1.so.insert(idx2);
                                t1.sow[idx2] = w;
                                t2.si.insert(idx1);
                                t2.siw[idx1] = w;
			}
			
			// if t2 is covered by t1 and in tail, then add edge t2 -> t1
			else if(t1.exons[t1.exons.size()-1].first == t2.exons[t2.exons.size()-1].first)
                        {
				int w = get_edge_weight(t1, t2, same_exon_num);
                                t2.so.insert(idx1);
                                t2.sow[idx1] = w;
                                t1.si.insert(idx2);
                                t1.siw[idx2] = w;
                        }
			return 1;	
		}
		*/
		
		break;
	}

	//printf("\n");
	//if(if_add_edge == true) return 1;
	
	/*
	// add edge if only few exons are different	
	if(e1 >= e2) return 0;
	
	set<PI32> es2(exons2.begin(), exons2.end());
	double diff_num = 0;
	double overlap_num = 0;
	
	for(int i = 0; i < exons1.size() - 1; i++)
	{
		if(exons1[i].second < exons2[0].first) continue;

		if(es2.find(exons1[i]) == es2.end()) diff_num++;
		overlap_num++;

	}

	if(diff_num/overlap_num <= 0.2)
	{
		int w = get_edge_weight(t1, t2, same_exon_num);
                t1.so.insert(idx2);
                t1.sow[idx2] = w;
                t2.si.insert(idx1);
                t2.siw[idx1] = w;
	}
	*/
	

	return 0;
}


bool cmp(pair<int, int> &a, pair<int, int> &b)
{
	return a.second > b.second;
}

set<int> get_top(map<int, int> &m)
{
	vector< pair<int,int> > A;
	for(auto it : m)
	{
		A.push_back(it);
	}
	sort(A.begin(), A.end(), cmp);

	set<int> s;
	for(int i = 0; i < max_path; i++) s.insert(A[i].first);

	return s;
}

void dfs_left(genome1 &gm, vector< vector<int> > &paths, vector<int> &temp, int i, set<int> &visited)
{
	visited.insert(i);
	temp.push_back(i);
	bool can_extend = false;
	
	transcript t = gm.transcripts[i];
	set<int>::iterator itr;

	set<int> s = t.so;
	if(s.size() > max_path) s = get_top(t.sow);

	for (itr = s.begin(); itr != s.end(); itr++)
	{
		if(visited.count(*itr) == 0)
		{
			can_extend = true;
			dfs_left(gm, paths, temp, *itr, visited);
		}
	}

	if(!can_extend && temp.size() > 1)
	{
		// if last t is still incomplete, skip
		transcript t_last = gm.transcripts[temp[temp.size()-1]];
		
		if(empty_side(t_last) != "both" && empty_side(t_last) != "right")
		{
			paths.push_back(temp);
		}
	}

	visited.erase(i);
	temp.pop_back();

}

int print_subgraph(genome1 &gm, int subgraph_idx, string ofile)
{
	ofstream fout(ofile.c_str());
	vector<int> sg = gm.subgraph[subgraph_idx];

	string vname; // v + v_idx
	string tname; // tid

	string vinfo;
	string einfo;

	fout<<"digraph{"<<endl;

	for(int i = 0; i < sg.size(); i++)
	{
		transcript t = gm.transcripts[sg[i]];

		vinfo = "\t" + t.transcript_id;
		//vinfo = "\tv" + tostring(i);
		fout<<vinfo<<endl;
	}

	for(int i = 0; i < sg.size(); i++)
        {
                transcript t = gm.transcripts[sg[i]];

		for(auto j : t.so)
		{
			einfo = "\t" + t.transcript_id + " -> " + gm.transcripts[j].transcript_id + " [label=" + tostring(t.sow[j]) + "]";
                	//einfo = "\tv" + tostring(i) + " -> v" + tostring(j) + " [label=" + tostring(t.sow[j]) + "]";
                	fout<<einfo<<endl;
		}
        }

	fout<<"}"<<endl;

	return 0;
}

int print_graph(genome1 &gm, string prefix)
{
        string d = prefix + "_subgraph";
        string cmd = "mkdir -p " + d;
        const char *command = cmd.c_str();
        int sys = system(command);

        for(int i = 0; i < gm.subgraph.size(); i++)
        {
                string f = d + "/geneI" + tostring(i) + ".gv";
                print_subgraph(gm, i, f);
        }

        return 0;
}

int build_exonlist(genome1 &gm, const vector<int> &p, vector< pair<int32_t, int32_t> > &exons)
{
        int idx1 = 0;
        int idx2 = 1;

        exons.push_back(gm.transcripts[p[idx1]].exons[0]);

        while(idx2 < p.size())
        {
                vector< pair<int32_t, int32_t> > exons1 = gm.transcripts[p[idx1]].exons;
                vector< pair<int32_t, int32_t> > exons2 = gm.transcripts[p[idx2]].exons;

                if(gm.transcripts[p[idx1]].end < gm.transcripts[p[idx2]].start) 
                {
                        for(int i = 1; i < exons1.size(); i++) exons.push_back(exons1[i]);
                        exons[exons.size()-1].second = exons2[0].second;
                        idx1++;
                        idx2++;
                        continue;
                }

                double cov1 = gm.transcripts[p[idx1]].coverage;
                double cov2 = gm.transcripts[p[idx2]].coverage;

                for(int i = 1; i < exons1.size(); i++)
                {

                        if(exons1[i].second <= exons2[0].second + max_shift_bp)
                        {
                                exons.push_back(exons1[i]);
                                
                                
                                if(exons1[i].second == exons2[0].second) break;
                                if(i == exons1.size() - 1) exons[exons.size()-1].second = exons2[0].second;
                        }
                        else
                        {
                                //printf("error! exons1[i].second > exons2[0].second + max_shift_bp\n");
                                if(i == exons1.size() - 1 && exons1[i].first < exons2[0].second) 
                                {
                                        exons.push_back(exons1[i]);
                                        exons[exons.size()-1].second = exons2[0].second;
                                }
                                break;
                        }
                        
                }

                idx1++;
                idx2++;
        }

        vector< pair<int32_t, int32_t> > exons1 = gm.transcripts[p[idx1]].exons;
        for(int i = 1; i < exons1.size(); i++) exons.push_back(exons1[i]);
        
        return 0;
}

void dfs_assign_gid(int subgraph_idx, genome1 &gm, set<int> &visited, int i, set<int> &subgraph)
{
	transcript &t = gm.transcripts[i];

	visited.insert(i);
	subgraph.insert(i);

	t.gene_id = "geneI" + tostring(subgraph_idx);
	if(gm.subgraph.size() < subgraph_idx + 1) gm.subgraph.push_back({i});
	else gm.subgraph[subgraph_idx].push_back(i);

	t.transcript_id = "geneI" + tostring(subgraph_idx) + "I" + tostring(gm.subgraph[subgraph_idx].size());

	for(int n : t.si)
	{
		if(visited.count(n) == 0) dfs_assign_gid(subgraph_idx, gm, visited, n, subgraph);
	}
	for(int n : t.so)
	{
		if(visited.count(n) == 0) dfs_assign_gid(subgraph_idx, gm, visited, n, subgraph);
	}

}

int assign_gid(genome1 &gm, vector< set<int> > &subgraphs)
{
	set<int> visited;
    int subgraph_idx = 0;
	set<int> subgraph;

	int iso_cnt = 0;
	int multi_cnt = 0;

	for(int i = 0; i < gm.transcripts.size(); i++)
	{	

		if(gm.transcripts[i].si.size() == 0 && gm.transcripts[i].so.size() == 0)
		{
			iso_cnt++;
		}
		//if(gm.transcripts[i].si.size() == 0 && gm.transcripts[i].so.size() == 0) continue;

		if(visited.count(i) == 0)
		{
			dfs_assign_gid(subgraph_idx, gm, visited, i, subgraph);
			subgraphs.push_back(subgraph);

			subgraph.clear();
			subgraph_idx++;
		}
	}

	printf("# graphs = %ld, isolated node/subgraph = %d\n", subgraphs.size(), iso_cnt);

	return 0;
}

bool is_sg_covered(const set<PI32>& a, const set<PI32>& b) 
{
	return false;

	int max_diff_junc_num = 1;
    	int diff_junc_num = 0;

	if (a.size() > b.size()) {
        	return false;
    	}
    	
	for (PI32 x : a) {
        	if (b.find(x) == b.end()) {
        		diff_junc_num++;    	
        	}
    	}
	
    	if(diff_junc_num <= max_diff_junc_num) 
	{
		return true;
	}

	return false;
}

int rm_duplicate_subgraph(genome1 &gm, vector< set<int> > &subgraphs, vector<bool> &sg_covered, vector< set<PI32> > &sgj)
{
	// keep all subgraph
	for(int i = 0; i < subgraphs.size(); i++) sg_covered.push_back(false);
	return 0;

	for(int k = 0; k < subgraphs.size(); k++)
	{
		set<PI32> jlist;
		set<int> &sg = subgraphs[k];
		
		for(int v : sg)
		{
			if(v == -1 || v == -2) continue;
			vector<PI32> &exons = gm.transcripts[v].exons;
                	for(int i = 0; i < exons.size()-1; i++)
                	{
                        	jlist.insert({exons[i].second, exons[i+1].first});
			}
                }

		sgj.push_back(jlist);
	}

	for(int i = 0; i < subgraphs.size(); i++) sg_covered.push_back(false);

	
	for(int i = 0; i < subgraphs.size(); i++)
	{
		if(sg_covered[i]) continue;

		for(int j = i+1; j < subgraphs.size(); j++)
		{
			if(is_sg_covered(sgj[i], sgj[j]))
			{
				sg_covered[i] = true;
				printf("subgraph %d is covered by %d, %ld, %ld\n", i, j, sgj[i].size(), sgj[j].size());
				//for(PI32 jc : sgj[i]) cout<<"("<<jc.first<<","<<jc.second<<") ";
				//printf("\n");
				//for(PI32 jc : sgj[j]) cout<<"("<<jc.first<<","<<jc.second<<") "<<endl;
				//printf("\n");
				break;
			}
			if(is_sg_covered(sgj[j], sgj[i])) 
			{
				sg_covered[j] = true;
				//printf("subgraph %d is covered by %d, %d, %d\n", j, i, sgj[j].size(), sgj[i].size());
				//for(PI32 jc : sgj[j]) cout<<"("<<jc.first<<","<<jc.second<<") ";
                                //printf("\n");
                                //for(PI32 jc : sgj[i]) cout<<"("<<jc.first<<","<<jc.second<<") "<<endl;
                                printf("\n");
			}
		}
	}
	

	return 0;
}

void run_dp(genome1 &gm, vector<int> &sg, map< PI32, double> &W, map< int, set<PI32> > &LV, vector<int> &p, set< PI32 > &js, set<int> &vtot)
{
	/*		
	for(int a : sg) 
	{
		if(a == -1 || a == -2) cout<<a<<" ";
		else cout<<gm.transcripts[a].transcript_id<<" ";
	}
	cout<<"; vtot size = "<<vtot.size()<<endl;
	*/

	int n = sg.size();
	
	map<int, int> vIdx;
	for(int i = 0; i < n; i++) vIdx[sg[i]] = i;

	vector<double> dp(n, -100);
	vector<double> dpcov(n, -100);
	vector< set<PI32> > dpj(n);
	vector<int> prev(n, -100);

	dp[0] = 0;
	dpcov[0] = 0;

	for(int i = 1; i < n; i++)
	{
		int v = sg[i];
		set<PI32> vj = LV[v];

		set<int> iedges = vtot;
		if(v != -2) iedges = gm.transcripts[v].si;

		for(auto u : iedges)
		{
			int uidx = vIdx[u];
			set<PI32> uniqj = {};
			set<PI32> uj = dpj[uidx];

			//printf("checking %d -> %d\n", u, v);

			double addw = 0;

			for(auto j : vj)
			{
				if(uj.find(j) == uj.end()) 
				{
					uniqj.insert(j);
					addw += W[j];
				}
			}
			
			//cout<<dp[i]<<" vs "<<dp[uidx] + addw<<endl;

			if(dp[i] < dp[uidx] + addw)
			{
				dp[i] = dp[uidx] + addw;
				prev[i] = uidx;
				//cout<<"update prev["<<i<<"] = "<<prev[i]<<endl;
				dpj[i] = dpj[uidx];
				dpj[i].insert(uniqj.begin(), uniqj.end());
				
				if(u == -1) dpcov[i] = gm.transcripts[v].coverage;
				else if(v == -2) dpcov[i] = dpcov[uidx];
				else dpcov[i] = min(gm.transcripts[v].coverage, dpcov[uidx]);
			}
			else if(dp[i] == dp[uidx] + addw && dpcov[uidx] > dpcov[prev[i]])
			{
				dp[i] = dp[uidx] + addw;
                                prev[i] = uidx;
                                //cout<<"update prev["<<i<<"] = "<<prev[i]<<endl;
                                dpj[i] = dpj[uidx];
                                dpj[i].insert(uniqj.begin(), uniqj.end());

				if(u == -1) dpcov[i] = gm.transcripts[v].coverage;
                                else if(v == -2) dpcov[i] = dpcov[uidx];
                                else dpcov[i] = min(gm.transcripts[v].coverage, dpcov[uidx]);
			}
		}
	}

	js = dpj[n-1];
	int t = n-1;

	while(t != 0)
	{
		p.push_back(sg[prev[t]]);
		t = prev[t];
	}
	p.pop_back();
	reverse(p.begin(), p.end());

	/*
	printf("find path: ");
	for(int i : p) printf("%s, ", gm.transcripts[i].transcript_id.c_str());
	printf("\n");
	*/
}

void get_path_in_subgraph(genome1 &gm, vector<int> &sg, vector< vector<int> > &paths, set<int> &vtot, set<PI32> L)
{
	map< PI32, double> W; // weights of all junctions
	map< int, set<PI32> > LV; // junction list of each v
	LV[-1] = {};
	LV[-2] = {};

	for(int v : sg)
	{
		if(v == -1 || v == -2) continue;

		vector<PI32> &exons = gm.transcripts[v].exons;
		for(int i = 0; i < exons.size()-1; i++)
		{
			if(LV.find(v) == LV.end()) LV[v] = {{exons[i].second, exons[i+1].first}};
			else LV[v].insert({exons[i].second, exons[i+1].first});
		}
	}

	for(PI32 l : L) W[l] = 1;

	set< PI32 > visited;
	set< vector<int> > ps;

	/*
	set<int32_t> B;
	for(PI32 l : L) 
	{
		B.insert(l.first);
		B.insert(l.second);
	}
	set<int32_t> visitedB;

	while(visitedB.size() < B.size())
	*/

	while(visited.size() < L.size() -1)
	{
		//printf("visited size = %d, L size = %d\n", visited.size(), L.size());
		vector<int> p = {};
		set< PI32 > js = {};
		
		run_dp(gm, sg, W, LV, p, js, vtot);

		if(ps.find(p) != ps.end()) break;

		transcript &t = gm.transcripts[p[0]];
                transcript &t_last = gm.transcripts[p[p.size()-1]];
		if(empty_side(t) != "left" && empty_side(t) != "both" && empty_side(t_last) != "both" && empty_side(t_last) != "right")
		{
			ps.insert(p);
		}
		//else printf("path is filtered out due to imcomplete end\n");

		for(auto j : js)
		{
			visited.insert(j);
			W[j] = 0.1;

			//visitedB.insert(j.first);
			//visitedB.insert(j.second);
		}
		for(PI32 l : L)
		{
			if(visited.find(l) == visited.end()) W[l] = 10;
		}
	}

	for(auto pp : ps)
	{
		paths.push_back(pp);
		//printf("insert path #%d: ", paths.size()-1);
		//for(int i : pp) printf("%s, ", gm.transcripts[i].transcript_id.c_str());
        	//printf("\n");
	}
}

vector<int> topological_sort(vector<vector<int>> &edges)
{
	map<int, vector<int>> adj;
	map<int, int> inDegree;

	for (auto e : edges) 
	{
		adj[e[0]].push_back(e[1]);
		inDegree[e[1]]++;
	}

	queue<int> q;
	for (auto vertex : adj) {
		int u = vertex.first;
		if (inDegree[u] == 0) {
			q.push(u);
		}
	}

	vector<int> result;
	while (!q.empty()) {
		int u = q.front();
		q.pop();
		result.push_back(u);
		for (int v : adj[u]) {
			inDegree[v]--;
		    	if (inDegree[v] == 0) {
				q.push(v);
		    	}
		}
	}

	return result;
}

int find_allpaths_greedy(genome1 &gm, vector< set<int> > &subgraphs, vector< vector<int> > &paths, vector<bool> &sg_covered, vector< set<PI32> > &sgj)
{
	for(int i = 0; i < subgraphs.size(); i++)
	{
		if(sg_covered[i] == true) continue;
		if(sgj[i].size() == 0) continue;
		//printf("finding path in subgraph #%d\n", i);

		set<int> sg = subgraphs[i];
		vector< vector<int> > edges;
		set<int> vtot;

		for(int v : sg)
		{
			if(gm.transcripts[v].si.size() == 0) 
			{
				gm.transcripts[v].si.insert(-1);
				edges.push_back({-1, v});
				for(int u : gm.transcripts[v].so) edges.push_back({v, u});
				continue;
			}
			if(gm.transcripts[v].so.size() == 0)
			{
				gm.transcripts[v].so.insert(-2);
				edges.push_back({v, -2});
				vtot.insert(v);
				continue;
			}
			for(int u : gm.transcripts[v].so) edges.push_back({v, u});
		}

		vector<int> sortedsg = topological_sort(edges);

		//printf("subgraph #%d, # v = %d, # sortedv = %d, # e = %d\n", i, sg.size(), sortedsg.size(), edges.size());
		
		get_path_in_subgraph(gm, sortedsg, paths, vtot, sgj[i]);
	}

	return 0;
}

double calculate_score_abu1(genome1 &gm, transcript &curt)
{
	double score = -1;
	vector<PI32> &exons = curt.exons;

	for(int i = 0; i < exons.size() - 1; i++)
	{
		double curScore = 0;

		string h = curt.seqname;
		h.append(tostring(exons[i].second));
		h.append(tostring(exons[i+1].first));

		set<int> &tlist = gm.tjlist[h];

		for(int tidx : tlist)
		{
			//if(curt.strand != gm.transcripts[tidx].strand) continue;
			if(is_covered(gm.transcripts[tidx], curt)) curScore += gm.transcripts[tidx].coverage;
		}

		if(score == -1) score = curScore;
		else score = min(score, curScore);
	}

	return score;
}

double calculate_score_extend(genome1 &gm, transcript &curt)
{
	double score = 0;

	vector<PI32> &exons = curt.exons;

	set<int> visited;

	double times = 3;

	double times1 = 1;

	for(int i = 0; i < exons.size() - 1; i++)
	{
		string h = curt.seqname;
		h.append(tostring(exons[i].second));
		h.append(tostring(exons[i+1].first));

		set<int> &tlist = gm.tjlist[h];
		for(int tidx : tlist)
		{
			if(visited.find(tidx) != visited.end()) continue;
			visited.insert(tidx);

			transcript &t = gm.transcripts[tidx];
			
			// skip if t is inside of curt
			if(t.exons[0].second >= exons[0].second && t.exons.back().first <= exons.back().first)
			{
				/*
				if(t.exons[0].second == exons[0].second && t.exons.back().first == exons.back().first)
				{
					if(is_covered(t, curt)) score += times1 * (t.empty_style[1] + t.empty_style[2] + t.empty_style[3]);
					continue;
				}
				if(t.exons[0].second == exons[0].second && t.exons.back().first < exons.back().first) 
				{
					if(is_covered(t, curt)) score += times1 * (t.empty_style[1] + t.empty_style[2]);
					continue;
				}
				if(t.exons[0].second > exons[0].second && t.exons.back().first == exons.back().first)
				{
					if(is_covered(t, curt)) score += times1 * (t.empty_style[1] + t.empty_style[3]);
					continue;
				}
				*/
				
				continue;
			}
			
			// curt is inside of t
			if(t.exons[0].second <= exons[0].second && t.exons.back().first >= exons.back().first)
			{
				if(is_covered(curt, t)) score += times * t.coverage;
				continue;
			}

			// t can extend to curt
			if(t.exons[0].second < exons[0].second)
			{
				if(is_extendable(t, curt)) score += times * t.coverage;
				continue;
			}

			// curt can extend to t
			if(t.exons.back().first > exons.back().first)
			{
				if(is_extendable(curt, t)) score += times * t.coverage;
				continue;
			}
		}
	}

	return score;
}

double calculate_score_abu(genome1 &gm, transcript &curt, set<int> &compatible)
{
        double score = -1;

        vector<PI32> &exons = curt.exons;

        for(int i = 0; i < exons.size() - 1; i++)
        {
			double curScore = 0;

			string h = curt.seqname;
			h.append(tostring(exons[i].second));
			h.append(tostring(exons[i+1].first));

			set<int> &tlist = gm.tjlist[h];
			//printf("tlist size = %ld\n", tlist.size());

			for(int tidx : tlist)
			{
				//if(curt.strand != gm.transcripts[tidx].strand) continue;
				if(compatible.find(tidx) != compatible.end()) continue;
				
				if((gm.transcripts[tidx].exons[0].second == curt.exons[0].second) && (gm.transcripts[tidx].sb != curt.sb)) continue;
				if((gm.transcripts[tidx].exons.back().first == curt.exons.back().first) && (gm.transcripts[tidx].eb != curt.eb)) continue;

				//print_transcript_pair(gm.transcripts[tidx], curt);
				if(is_covered(gm.transcripts[tidx], curt)) 
				{
					/*
					if(gm.transcripts[tidx].exons[0].second == curt.exons[0].second) 
					{
						if(find_interval(gm.transcripts[tidx].start, gm.jblist[curt.seqname]) != find_interval(gm.transcripts[tidx].start, gm.jblist[curt.seqname])) continue;
					}
					if(gm.transcripts[tidx].exons.back().first ==  curt.exons.back().first)
					{
						if(find_interval(gm.transcripts[tidx].end, gm.jblist[curt.seqname]) != find_interval(gm.transcripts[tidx].end, gm.jblist[curt.seqname])) continue;
					}
					*/

					curScore += gm.transcripts[tidx].coverage;
				}
			}

			if(score == -1) score = curScore;
			else score = min(score, curScore);
        }

        return score;
}

double calculate_score_abu_ej(genome1 &gm, transcript &curt, set<int> &compatible)
{
        double score = -1;

        vector<PI32> &exons = curt.exons;

        for(int i = 0; i < exons.size() - 1; i++)
        {
                double curScore = 0;

                string h = curt.seqname;
                h.append(tostring(exons[i].second));
                h.append(tostring(exons[i+1].first));

                set<int> &tlist = gm.tjlist[h];
                //printf("tlist size = %ld\n", tlist.size());

                for(int tidx : tlist)
                {
			//if(curt.strand != gm.transcripts[tidx].strand) continue;
                        if(compatible.find(tidx) != compatible.end()) continue;

                        //print_transcript_pair(gm.transcripts[tidx], curt);
                        if(is_covered(gm.transcripts[tidx], curt)) curScore += gm.transcripts[tidx].coverage;
                }

                if(score == -1) score = curScore;
                else score = min(score, curScore);

        }

	for(int i = 1; i < exons.size() - 1; i++)
        {
                double curScore = 0;

                string h = curt.seqname;
                h.append(tostring(exons[i].first));
                h.append(tostring(exons[i].second));

                set<int> &tlist = gm.telist[h];
                //printf("tlist size = %ld\n", tlist.size());

                for(int tidx : tlist)
                {
                        if(compatible.find(tidx) != compatible.end()) continue;
			//if(curt.strand != gm.transcripts[tidx].strand) continue;

                        //print_transcript_pair(gm.transcripts[tidx], curt);
                        if(is_covered(gm.transcripts[tidx], curt)) curScore += gm.transcripts[tidx].coverage;
                }
		if(curScore == 0)
		{
			cout<<"tid "<<curt.transcript_id<<", exon #"<<i<<" "<<h<<"  coverage = 0"<<endl;
		}

                if(score == -1) score = curScore;
                else score = min(score, curScore);

        }

        return score;
}


double calculate_score_abu_sumcov(genome1 &gm, transcript &curt, set<int> &compatible)
{
        double score = 0;

        vector<PI32> &exons = curt.exons;

	set<int> visited;

        for(int i = 0; i < exons.size() - 1; i++)
        {
                string h = curt.seqname;
                h.append(tostring(exons[i].second));
                h.append(tostring(exons[i+1].first));

                set<int> &tlist = gm.tjlist[h];
                //printf("tlist size = %ld\n", tlist.size());

                for(int tidx : tlist)
                {
			//if(curt.strand != gm.transcripts[tidx].strand) continue;
                        if(compatible.find(tidx) != compatible.end()) continue;
			if(visited.find(tidx) != visited.end()) continue;

                        //print_transcript_pair(gm.transcripts[tidx], curt);
                        if(is_covered(gm.transcripts[tidx], curt)) score += gm.transcripts[tidx].coverage;
			visited.insert(tidx);
                }
        }

        return score;
}


double calculate_score_abuj(genome1 &gm, transcript &curt, set<int> &compatible)
{
	double score = calculate_score_abu(gm, curt, compatible);
	score *= (curt.exons.size() - 1);

    return score;
}

int calculate_score(genome1 &gm, transcript &curt)
{
	//printf("calculating score...\n");

	int score = -1;

	vector<PI32> &exons = curt.exons;
	//printf("exons size = %ld\n", exons.size());

	for(int i = 0; i < exons.size() - 1; i++)
	{
		int curScore = 0;

		string h = curt.seqname;
                h.append(tostring(exons[i].second));
                h.append(tostring(exons[i+1].first));

		set<int> &tlist = gm.tjlist[h];
		//printf("tlist size = %ld\n", tlist.size());

		for(int tidx : tlist)
		{
			//if(curt.strand != gm.transcripts[tidx].strand) continue;
			//print_transcript_pair(gm.transcripts[tidx], curt);
			if(is_covered(gm.transcripts[tidx], curt)) curScore++;
		}
		
		if(score == -1) score = curScore;
		else score = min(score, curScore);
	}

	return score;
}

void get_seeds(genome1 &gm, set<int> &sg, vector<int> &seeds, double &seedScore, set<int> &compatible)
{
	//int max = -1;
	double max = -1;

	for(int i : sg)
	{
		if(i == -1 || i == -2) continue;

		//int score = calculate_score(gm, gm.transcripts[i]);
		double score = calculate_score_abu(gm, gm.transcripts[i], compatible);
	
		if(score > max) 
		{
			seeds.clear();
			max = score;
			seeds.push_back(i);
		}
		else if(score == max) seeds.push_back(i);
	}

	seedScore = max;

	/*
	cout<<"max score = "<<max<<", seed list: ";
        for(int s : seeds) cout<<s<<" ";
        printf("\n");
	*/

}

transcript path_to_texons(genome1 &gm, vector<int> &p)
{
	transcript t;
	
	t.seqname = gm.transcripts[p[0]].seqname;
	t.start = gm.transcripts[p[0]].start+1;
	t.end = gm.transcripts[p.back()].end;
	t.strand = gm.transcripts[p[0]].strand;

	t.sb = gm.transcripts[p[0]].sb;
	t.eb = gm.transcripts[p.back()].eb;

	vector<PI32> exons;
	build_exonlist(gm, p, exons);
	t.exons = exons;


	return t;
}

void dfs_to_right(genome1 &gm, vector<int> &p, vector<double> &ps, set<int> &compatible)
{
	int cur = p.back();
	double maxScore = -1;
	int nextid = -2;

	for(int rid : gm.transcripts[cur].so)
	{
		if(rid == -2) continue;

		if(compatible.find(rid) != compatible.end()) continue;

		p.push_back(rid);
		transcript t = path_to_texons(gm, p);
		double tempScore = calculate_score_abu(gm, t, compatible);
		if(tempScore > maxScore)
		{
			nextid = rid;
			maxScore = tempScore;
		}
		p.pop_back();
	}

	if(nextid != -2)
	{
		for(double s : ps) cout<<s<<",";
                printf("\n");
	}
	double drop_threshold = 0.2;
	if(nextid != -2 && (maxScore >= drop_threshold * ps.back()))
	{
		p.push_back(nextid);
		ps.push_back(maxScore);
		//for(double s : ps) cout<<s<<",";
		//printf("\n");
		dfs_to_right(gm, p, ps, compatible);
	}
}

void dfs_to_left(genome1 &gm, vector<int> &p, vector<double> &ps, set<int> &compatible)
{
        int cur = p[0];
        double maxScore = -1;
        int nextid = -1;

        for(int lid : gm.transcripts[cur].si)
        {       
                if(lid == -1) continue;

		if(compatible.find(lid) != compatible.end()) continue;

                p.insert(p.begin(), lid);
                transcript t = path_to_texons(gm, p);
                double tempScore = calculate_score_abu(gm, t, compatible);
                if(tempScore > maxScore)
                {       
                        nextid = lid;
                        maxScore = tempScore;
                }       
                p.erase(p.begin());
        }       

	if(nextid != -1)
	{
		for(double s : ps) cout<<s<<",";
                printf("\n");
	}

	double drop_threshold = 0.2;
        if(nextid != -1 && (maxScore >= drop_threshold * ps[0]))
        {
                p.insert(p.begin(), nextid);
		ps.insert(ps.begin(), maxScore);
		//for(double s : ps) cout<<s<<",";
                //printf("\n");
		dfs_to_left(gm, p, ps, compatible);
        }
}

void get_seed_path(genome1 &gm, int seed, vector<int> &p, vector<double> &ps, set<int> &compatible)
{
	vector<int> pright = {seed};
	vector<double> psright = {ps[0]};
	vector<int> pleft = {seed};
	vector<double> psleft = {ps[0]};

	dfs_to_right(gm, pright, psright, compatible);
	dfs_to_left(gm, pleft, psleft, compatible);

	p = pleft;
	if(pright.size() > 1)
	{
		for(int i = 1; i < pright.size(); i++) p.push_back(pright[i]);
	}

	ps = psleft;
	if(psright.size() > 1)
	{
		for(int i = 1; i < psright.size(); i++) ps.push_back(psright[i]);
	}
}

int update_compatible(genome1 &gm, transcript &curt, set<int> &compatible)
{
	vector<PI32> &exons = curt.exons;

	for(int i = 0; i < exons.size() - 1; i++)
        {
                int curScore = 0;

                string h = curt.seqname;
                h.append(tostring(exons[i].second));
                h.append(tostring(exons[i+1].first));

                set<int> &tlist = gm.tjlist[h];
                //printf("tlist size = %ld\n", tlist.size());

                for(int tidx : tlist)
                {
                        //print_transcript_pair(gm.transcripts[tidx], curt);
                        if(is_covered(gm.transcripts[tidx], curt)) compatible.insert(tidx);
                }
	}

	return 0;
}

int find_allpaths_greedy_extend(genome1 &gm, vector< set<int> > &subgraphs, vector< vector<int> > &paths, vector< vector<double> > &pscores, vector<bool> &sg_covered)
{
	for(int i = 0; i < subgraphs.size(); i++)
        {
                if(sg_covered[i] == true) continue;

		//printf("checking subgraph %d...\n", i);

		set<int> sg = subgraphs[i];
		
		set<int> compatible; // all nodes that are compatible with previously determined paths
		vector<vector<int> > tempp;
        vector<vector<double> > tempps;

		while(tempp.size() < 5)
		{
			vector<int> seeds = {};
			double seedScore = 0;
			get_seeds(gm, sg, seeds, seedScore, compatible);

			if(seedScore < min_start_score) break;

			for(int seed : seeds)
			{
				vector<int> p = {};
				vector<double> ps = {seedScore};

				get_seed_path(gm, seed, p, ps, compatible);

				tempp.push_back(p);
				tempps.push_back(ps);

				transcript curt = path_to_texons(gm, p);
				update_compatible(gm, curt, compatible);
			}
		}

		for(int i = 0; i < tempp.size(); i++)
		{
			vector<int> p = tempp[i];	
			//printf("seed %d, path (", seed);
			//for(int pp : p) cout<<pp<<" ";
			//cout<<")"<<endl;

			string es = empty_side(gm.transcripts[p[0]]);
			if(es == "left" || es == "both") continue;
			es = empty_side(gm.transcripts[p.back()]);
			if(es == "right" || es == "both") continue;			

			paths.push_back(p);
			pscores.push_back(tempps[i]);
		}
	}
	return 0;
}

void dfs_to_right_abuj(genome1 &gm, vector<int> &p, vector<double> &ps, set<int> &compatible)
{
        int cur = p.back();
        double maxScore = -1;
        int nextid = -2;

        for(int rid : gm.transcripts[cur].so)
        {
                if(rid == -2) continue;

                if(compatible.find(rid) != compatible.end()) continue;

                p.push_back(rid);
                transcript t = path_to_texons(gm, p);
                double tempScore = calculate_score_abuj(gm, t, compatible);
                if(tempScore > maxScore)
                {
                        nextid = rid;
                        maxScore = tempScore;
                }
                p.pop_back();
        }

	/*
        if(nextid != -2)
        {
                for(double s : ps) cout<<s<<",";
                printf("\n");
        }
	*/

        double drop_threshold = 0.2;
        if(nextid != -2 && (maxScore >= drop_threshold * ps.back()))
        {
                p.push_back(nextid);
                ps.push_back(maxScore);
		//for(double s : ps) cout<<s<<",";
                //printf("\n");
                dfs_to_right_abuj(gm, p, ps, compatible);
        }
}

void dfs_to_left_abuj(genome1 &gm, vector<int> &p, vector<double> &ps, set<int> &compatible)
{
        int cur = p[0];
        double maxScore = -1;
        int nextid = -1;

        for(int lid : gm.transcripts[cur].si)
        {
                if(lid == -1) continue;

                if(compatible.find(lid) != compatible.end()) continue;

                p.insert(p.begin(), lid);
                transcript t = path_to_texons(gm, p);
                double tempScore = calculate_score_abuj(gm, t, compatible);
                if(tempScore > maxScore)
                {
                        nextid = lid;
                        maxScore = tempScore;
                }
                p.erase(p.begin());
        }


        double drop_threshold = 0.2;
        if(nextid != -1 && (maxScore >= drop_threshold * ps[0]))
        {
                p.insert(p.begin(), nextid);
                ps.insert(ps.begin(), maxScore);
                //for(double s : ps) cout<<s<<",";
                //printf("\n");
                dfs_to_left_abuj(gm, p, ps, compatible);
        }
}


void get_seed_path_abuj(genome1 &gm, int seed, vector<int> &p, vector<double> &ps, set<int> &compatible)
{
	vector<int> pright = {seed};
        vector<double> psright = {ps[0]};
        vector<int> pleft = {seed};
        vector<double> psleft = {ps[0]};

        dfs_to_right_abuj(gm, pright, psright, compatible);
        dfs_to_left_abuj(gm, pleft, psleft, compatible);

        p = pleft;
        if(pright.size() > 1)
        {
                for(int i = 1; i < pright.size(); i++) p.push_back(pright[i]);
        }

        ps = psleft;
        if(psright.size() > 1)
        {
                for(int i = 1; i < psright.size(); i++) ps.push_back(psright[i]);
        }
}


int find_allpaths_greedy_extend_abuj(genome1 &gm, vector< set<int> > &subgraphs, vector< vector<int> > &paths, vector< vector<double> > &pscores, vector<bool> &sg_covered)
{

	for(int i = 0; i < subgraphs.size(); i++)
	{
		if(sg_covered[i] == true) continue;

		printf("\nchecking subgraph %d...\n", i);
		printf("node list: (");
		for(int k : subgraphs[i]) cout<<k<<" ";
		cout<<")"<<endl;


		set<int> sg = subgraphs[i];

		set<int> compatible; // all nodes that are compatible with previously determined paths
		vector<vector<int> > tempp;
		vector<vector<double> > tempps;

		while(tempp.size() < 3)
		{
			vector<int> bestp;
                	double bestScore = -1;
                	vector<double> bestps;
			transcript bestt;

			set<int> visited;
			
			bool anyseed = false;

			for(int i : sg)
			{
				printf("seed node: %d\n", i);
				if(i == -1 || i == -2) continue;
				if(compatible.find(i) != compatible.end()) continue;

				if(visited.find(i) != visited.end()) continue;
				visited.insert(i);

				double seedScore = calculate_score_abuj(gm, gm.transcripts[i], compatible);
				if(seedScore < (gm.transcripts[i].exons.size()-1) * min_start_score) continue;

				anyseed = true;

				vector<int> p = {};
				vector<double> ps = {seedScore};

				get_seed_path_abuj(gm, i, p, ps, compatible);

				printf("find path: (");
                        	for(int pp : p) cout<<pp<<" ";
                        	cout<<")"<<endl;
				
				transcript t = path_to_texons(gm, p);
				double pScore = calculate_score_abuj(gm, t, compatible);
				cout<<"path score = "<<pScore<<endl;
				
				if(pScore > bestScore)
				{
					printf("find better path\n");
					bestp = p;
					bestps = ps;
					bestScore = pScore;
					bestt = t;
				}
			}

			if(!anyseed) break;

			tempp.push_back(bestp);
			tempps.push_back(bestps);
			update_compatible(gm, bestt, compatible);

			printf("finish round %ld...\n", tempp.size());

		}

		for(int i = 0; i < tempp.size(); i++)
		{
			vector<int> p = tempp[i];
			printf("path (");
			for(int pp : p) cout<<pp<<" ";
			cout<<")"<<endl;

			string es = empty_side(gm.transcripts[p[0]]);
			if(es == "left" || es == "both") continue;
			es = empty_side(gm.transcripts[p.back()]);
			if(es == "right" || es == "both") continue;

			paths.push_back(p);
			pscores.push_back(tempps[i]);
		}

	}
	return 0;
}

void get_top3_pq(priority_queue<pair<double, vector<int>>> &pq)
{
	priority_queue<pair<double, vector<int>>> temp;
	while(temp.size() < 3 && !pq.empty())
	{
		temp.push(pq.top());
		pq.pop();
	}
	pq = temp;
}

void get_top_rm_duplicate_pq(genome1 &gm, priority_queue<pair<double, vector<int>>> &pq)
{
	priority_queue<pair<double, vector<int>>> temp;
	set<double> scorelist;
	set<set<PI32>> iclist;
	while(temp.size() < num_path_per_node)
	{
		if(pq.empty()) break;
		double cur = pq.top().first;
                if(scorelist.find(cur) != scorelist.end())
                {
                        pq.pop();
                        continue;
                }

		vector<int> curp = pq.top().second;
                transcript t = path_to_texons(gm, curp);

                set<PI32> tj;
                for(int i = 0; i < t.exons.size()-1; i++)
                {
                        tj.insert({t.exons[i].second, t.exons[i+1].first});
                }

		if(iclist.find(tj) != iclist.end())
                {
                        pq.pop();
                        continue;
                }

                scorelist.insert(cur);
		iclist.insert(tj);

                temp.push(pq.top());
                pq.pop();
	}

	/*
	set<set<PI32>> icblist; // ic + {sb, eb}
	set<set<PI32>> iclist;

	while(iclist.size() < num_path_per_node)
	{	
		if(pq.empty()) break;

		vector<int> curp = pq.top().second;
		transcript t = path_to_texons(gm, curp);
		
		set<PI32> tj;
                for(int i = 0; i < t.exons.size()-1; i++)
                {
                        tj.insert({t.exons[i].second, t.exons[i+1].first});
                }

		if(iclist.find(tj) == iclist.end()) iclist.insert(tj);

		tj.insert({gm.transcripts[curp[0]].sb, gm.transcripts[curp.back()].eb});
		if(icblist.find(tj) != icblist.end())
		{
			pq.pop();
			continue;
		}
		scorelist.insert(cur);
		icblist.insert(tj);
		temp.push(pq.top());
		pq.pop();
	}

	*/

	pq = temp;
}

priority_queue<pair<double, vector<int>>> get_top_j_pq(genome1 &gm, priority_queue<pair<double, vector<int>>> pq)
{
	priority_queue<pair<double, vector<int>>> temp;
	
	/*
	set<PI32> jlist;
	map<set<PI32>, set<PI32>> jbmap;
	int num_new_j = 0;
	while(num_new_j < num_path_per_graph)
	{
		if(pq.empty()) break;

		pair<double, vector<int>> cur = pq.top();
                vector<int> &curp = cur.second;
		pq.pop();

		transcript t = path_to_texons(gm, cur.second);

                set<PI32> tj;
                for(int i = 0; i < t.exons.size()-1; i++)
                {
                        tj.insert({t.exons[i].second, t.exons[i+1].first});
                }
		
		bool any_new_j = false;
                for(auto j : tj)
                {
                        if(jlist.find(j) == jlist.end())
                        {
                                any_new_j = true;
                                break;
                        }
                }
		
		if(any_new_j)
                {
			num_new_j++;
                        temp.push(cur);
                        for(auto j : tj) jlist.insert(j);
			jbmap[tj].insert({gm.transcripts[curp[0]].sb, gm.transcripts[curp.back()].eb});
                }
		else
		{
			if(jbmap[tj].find({gm.transcripts[curp[0]].sb, gm.transcripts[curp.back()].eb}) == jbmap[tj].end())
			{
				temp.push(cur);
				jbmap[tj].insert({gm.transcripts[curp[0]].sb, gm.transcripts[curp.back()].eb});
			}
		}
	}
	*/


	set<PI32> jlist;
	while(temp.size() < num_path_per_graph)
	{

	if(pq.empty()) break;
	
	pair<double, vector<int>> cur = pq.top();
	pq.pop();
	
	transcript t = path_to_texons(gm, cur.second);
	
	set<PI32> tj;
		
	for(int i = 0; i < t.exons.size()-1; i++)
	{
		tj.insert({t.exons[i].second, t.exons[i+1].first});
	}

	bool any_new_j = false;
	for(auto j : tj)
	{
		if(jlist.find(j) == jlist.end())
		{
			any_new_j = true;
			break;
		}
	}

	if(any_new_j)
	{
				temp.push(cur);
		for(auto j : tj) jlist.insert(j);
	}
	}


	return temp;
}

priority_queue<pair<double, vector<int>>> resort_abu(genome1 &gm, priority_queue<pair<double, vector<int>>> pq)
{
	priority_queue<pair<double, vector<int>>> temp;

	while(!pq.empty())
	{
		vector<int> curp = pq.top().second;
		pq.pop();
		
		transcript t = path_to_texons(gm, curp);
		double curScore = calculate_score_abu1(gm, t);
		temp.push({curScore, curp});
	}

	return temp;
}

priority_queue<pair<double, vector<int>>> get_top_ic_pq(genome1 &gm, priority_queue<pair<double, vector<int>>> pq)
{
        priority_queue<pair<double, vector<int>>> temp;
        set<set<PI32>> iclist;

        while(temp.size() < num_path_per_graph)
        {
                if(pq.size() == 0) break;
                pair<double, vector<int>> cur = pq.top();
                pq.pop();

                transcript t = path_to_texons(gm, cur.second);

                set<PI32> tj;
                for(int i = 0; i < t.exons.size()-1; i++)
                {
                        tj.insert({t.exons[i].second, t.exons[i+1].first});
                }

                bool new_ic = false;
                if(iclist.find(tj) == iclist.end())
                {
                                new_ic = true;
                }

                if(new_ic)
                {
                        temp.push(cur);
                        iclist.insert(tj);
                }
        }

        return temp;
}

void update_related(genome1 &gm, int node, set<int> &related)
{
	related.insert(node);
	for(int v : gm.transcripts[node].so)
	{
		update_related(gm, v, related);
	}
}

void get_path_dp_abuj(genome1 &gm, vector<int> &sg, vector<vector<int>> &bestp, vector<double> &bestScore, vector<vector<double>> &bestps, set<int> &compatible)
{
	int n = sg.size();

	
	printf("sorted graph size = %d, node list = (", n);
        /*
	for(int node : sg) cout<<gm.transcripts[node].transcript_id<<", ";
        cout<<")"<<endl;
	*/

	map<int, int> vIdx;
	for(int i = 0; i < n; i++) vIdx[sg[i]] = i;

	vector< priority_queue<pair<double, vector<int>>> > dp;
	for(int i = 0; i < n; i++)
	{
		double tempScore = calculate_score_abuj(gm, gm.transcripts[sg[i]], compatible);
		priority_queue<pair<double, vector<int>>> init_pq;
		init_pq.push({tempScore, {sg[i]}});
		dp.push_back(init_pq);
	}

	
	for(int i = 0; i < n; i++)
	{
		cout<<gm.transcripts[sg[i]].transcript_id<<": init score = "<<dp[i].top().first<<endl;
	}
	

	for(int j = 1; j < n; j++) 
	{

		for(int u : gm.transcripts[sg[j]].si)
		{
			int uidx = vIdx[u];
			
			priority_queue<pair<double, vector<int>>> cur_pq = dp[uidx];

			while(!cur_pq.empty())
			{
				pair<double, vector<int>> it = cur_pq.top();
				cur_pq.pop();

				vector<int> itp = it.second;
				itp.push_back(sg[j]);
				transcript t = path_to_texons(gm, itp);
				double cur = calculate_score_abuj(gm, t, compatible);

				dp[j].push({cur, itp});
				get_top_rm_duplicate_pq(gm, dp[j]);
				//cout<<gm.transcripts[sg[j]].transcript_id<<": update score by "<<gm.transcripts[u].transcript_id<<"= "<<dp[j].top().first<<endl;
				//if(dp[j].size() > 3) get_top3_pq(dp[j]);
			}
		}
	}

	

	priority_queue<pair<double, vector<int>>> tempbest;

	for(int i = 0; i < n; i++)
	{	
		cout<<"node "<<i<<"/"<<n<<endl;
		for(int j = 0; j < num_path_per_node; j++)
		{
			if(dp[i].empty()) break;

			cout<<gm.transcripts[sg[i]].transcript_id<<": top "<<j<<" final score = "<<dp[i].top().first<<endl;
			double cur = dp[i].top().first;
			vector<int> curp = dp[i].top().second;

			dp[i].pop();

			string es = empty_side(gm.transcripts[curp[0]]);
			if(es == "left" || es == "both") continue;
			es = empty_side(gm.transcripts[curp.back()]);
			if(es == "right" || es == "both") continue;

			/*
			if(cur > bestScore)
			{
				bestp = curp;
				bestScore = cur;
			}
			*/

			tempbest.push({cur, curp});
		}
	}
	
	tempbest = get_top_j_pq(gm, tempbest);
	//tempbest = resort_abu(gm, tempbest);
	//tempbest = get_top_ic_pq(gm, tempbest);

	while(tempbest.size() > 0)
	{
		if(tempbest.top().first <= 0) break;

		vector<int> curp = tempbest.top().second;
		double curps = tempbest.top().first;

		bestp.push_back(curp);
		bestScore.push_back(curps);
		bestps.push_back({});
		
		cout<<"stored path score "<<curps<<endl;
		tempbest.pop();
	}
}

int find_allpaths_dp_abuj(genome1 &gm, vector< set<int> > &subgraphs, vector< vector<int> > &paths, vector< vector<double> > &pscores, vector<bool> &sg_covered)
{
	for(int i = 0; i < subgraphs.size(); i++)
	{
		if(sg_covered[i] == true) continue;
			
		printf("\nchecking subgraph %d...\n", i);

		set<int> sg = subgraphs[i];
		vector< vector<int> > edges;

		for(int v : sg)
		{
			if(gm.transcripts[v].si.size() == 0)
			{
				for(int u : gm.transcripts[v].so) edges.push_back({v, u});
				continue;
			}
			if(gm.transcripts[v].so.size() == 0)
			{
				continue;
			}
			for(int u : gm.transcripts[v].so) edges.push_back({v, u});
		}

		/*
		cout<<"edge list: ";
		for(vector<int> e : edges)
		{
			cout<<gm.transcripts[e[0]].transcript_id<<" -> "<<gm.transcripts[e[1]].transcript_id<<endl;
		}
		*/

		vector<int> sortedsg = topological_sort(edges);
		if(sg.size() == 1) 
		{
			for(int aa : sg) sortedsg.push_back(aa);
		}

		/*
		printf("sg size = %d/%d, node list = (", sg.size(), sortedsg.size());
		for(int node : sortedsg) cout<<gm.transcripts[node].transcript_id<<", ";
		cout<<")"<<endl;
		*/

		set<int> compatible; // all nodes that are compatible with previously determined paths
		//vector<vector<int> > tempp;
		//vector<vector<double> > tempps;

		//while(true)
		//while(tempp.size() < 1)
		//while(tempp.size() < 3)
		//{
			//vector<vector<int>> bestp;
			//vector<double> bestScore = {-1, -1, -1};
			//vector<vector<double>> bestps;
			//transcript bestt;

			//get_path_dp_abuj(gm, sortedsg, bestp, bestScore, bestps, compatible);
			//printf("round %d, find best path: (", tempp.size()+1);
			//for(int pp : bestp) cout<<gm.transcripts[pp].transcript_id<<", ";
			//cout<<"), score = "<<bestScore<<endl;

			//if(bestScore == -1) break;
			//if(bestScore == 0) break;

			//bestt = path_to_texons(gm, bestp);
			//if(bestScore >= bestt.exons.size() * min_start_score) 
			//{
			//tempp = bestp;
			//break;
			//}
			//else break;
			//update_compatible(gm, bestt, compatible);
		//}

		vector<vector<int>> bestp;
		vector<double> bestScore = {-1, -1, -1};
		vector<vector<double>> bestps;
		transcript bestt;

		get_path_dp_abuj(gm, sortedsg, bestp, bestScore, bestps, compatible);

		for(int i = 0; i < bestp.size(); i++)
		{
			vector<int> p = bestp[i];
			transcript t = path_to_texons(gm, p);

			printf("stored path (");
			for(int pp : p) cout<<gm.transcripts[pp].transcript_id<<" ";
			cout<<")"<<endl;

			string es = empty_side(gm.transcripts[p[0]]);
			if(es == "left" || es == "both") continue;
			es = empty_side(gm.transcripts[p.back()]);
			if(es == "right" || es == "both") continue;

			paths.push_back(p);
			//pscores.push_back(bestScore[i]);
		}
	}

	return 0;
}

int find_allpaths(genome1 &gm, vector< vector<int> > &paths)
{
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript t1 = gm.transcripts[i];
		
		// only check t when t has no in-edge and has out-edge, which means t would be a start node
		if(t1.si.size() > 0) continue;
		if(t1.so.size() == 0) continue;

		// skip if t is still incomplete
		if(empty_side(t1) == "left" || empty_side(t1) == "both") continue;

		set<int> visited;
		vector<int> temp;
		dfs_left(gm, paths, temp, i, visited);
	}
	
	return 0;
}

bool if_any_c_match(transcript &t1, transcript &t2)
{
	set<int> c1 = t1.clist;
	set<int> c2 = t2.clist;
	
	for(auto c : c1)
	{
		if(c2.find(c) != c2.end()) return true;
	}

	return false;
}

bool if_one_junctions(transcript &t1, transcript &t2, transcript &t3)
{
	if(t2.exons.size() <= 2) return true;
	
	int32_t s = t1.exons[t1.exons.size() - 1].first;
	int32_t e = t3.exons[0].second;

	if(t2.start >= s)
	{

	}

	return true;
}

bool if_borrow_one_junction(genome1 &gm, vector<int> &p)
{
	for(int i = 0; i < p.size() - 2 ; i+=2)
	{
		transcript t1 = gm.transcripts[p[i]];
		transcript t2 = gm.transcripts[p[i+1]];
		transcript t3 = gm.transcripts[p[i+2]];
		if(if_any_c_match(t1, t3) == false) return false;

		if(if_one_junctions(t1, t2, t3) == false) return false;
	}

	return true;
}

int allow_path_one_junction(genome1 &gm, vector< vector<int> > &paths)
{
	printf("# total paths = %ld\n", paths.size());

	vector< vector<int> > new_paths;

	for(int i = 0; i < paths.size(); i++)
	{
		if(paths[i].size() > 2) new_paths.push_back(paths[i]);
	}

	printf("# paths consists of > 2 transcripts = %ld\n", new_paths.size());
	
	paths = new_paths;
	new_paths.clear();

	for(int i = 0; i < paths.size(); i++)
        {
		if(if_borrow_one_junction(gm, paths[i]) == true) new_paths.push_back(paths[i]);	
	}

	printf("# paths allow n junction in gap = %ld\n", new_paths.size());
	paths = new_paths;
        new_paths.clear();

        return 0;
}

int build_clist(genome1 &gm, const vector<int> &p, set<int> &clist)
{
	for(int i = 0; i < p.size(); i++)
	{
		transcript t = gm.transcripts[p[i]];
		for(auto c : t.clist)
		{
			clist.insert(c);
		}
	}
	return 0;
}

// add cnum to clist if cell has any overlap on j
int build_clist_j(genome1 &gm, genome1 &gm3)
{
	map<string, set<int>> &jtocnum = gm3.jtocnum;

    for(int k = 0; k < gm.transcripts.size(); k++)
    {
		set<int> &cur_clist = gm.transcripts[k].clist;
		map<int, int> &cur_cjcount = gm.transcripts[k].cjcount; // <cell_id, #junction>

		vector<string> jlist;

		vector<PI32> &exons = gm.transcripts[k].exons;

		for(int i = 0; i < exons.size() - 1; i++)
		{
			string h = gm.transcripts[k].seqname;
			h.append(tostring(exons[i].second));
			h.append(tostring(exons[i+1].first));

			jlist.push_back(h);
		}

		for(int i = 0; i < jlist.size(); i++)
		{
			auto it = jtocnum.find(jlist[i]);
			if(it != jtocnum.end())
			{
				// TODO: comment below to get clist with a higher requirement
				//cur_clist.insert(it->second.begin(), it->second.end());
				for (const auto& item : it->second)
				{
					cur_cjcount[item] += 1;
				}
			}

		}

		// refine the clist and cjcount by adding requirements, e.g. half of j apprears
		//map<int, int> updated_cjcount;
		for (const auto& pair : cur_cjcount) 
		{
			int temp_cid = pair.first;
			int temp_jcount = pair.second;
			
			//if(temp_jcount >= 0.5 * jlist.size())
			if(temp_jcount >= 1)
			{
				//updated_cjcount[temp_cid] = temp_jcount;
				cur_clist.insert(temp_cid);
			}
		}
		//cur_cjcount = updated_cjcount;
    }
    return 0;
}


int output_transcript(ofstream &fout, genome1 &gm, genome1 &gm1, const vector<int> &p, const string tid, const vector<double> &ps)
{
	int psize = p.size();

	string chrm = gm.transcripts[p[0]].seqname.c_str();
	string source = "Beaver";
	char strand = gm.transcripts[p[0]].strand;
	string gid = gm.transcripts[p[0]].gene_id;
	double coverage = gm.transcripts[p[0]].coverage;
	
	string ttype = gm.transcripts[p[0]].transcript_type;
	string gtype = gm.transcripts[p[0]].gene_type;
	
	transcript t;

	t.seqname = chrm;
	t.source = source;
	t.feature = "transcript";
	t.gene_id = gid;
	t.transcript_id = tid;
	t.transcript_type = ttype;
	t.gene_type = gtype;
	t.start = gm.transcripts[p[0]].start+1;
	t.end = gm.transcripts[p.back()].end;
	t.strand = strand;
	t.frame = gm.transcripts[p[0]].frame;
	t.coverage = coverage;
	t.RPKM = gm.transcripts[p[0]].RPKM;
	t.FPKM = gm.transcripts[p[0]].FPKM;
	t.TPM = gm.transcripts[p[0]].TPM;
	t.hcov = gm.transcripts[p[0]].hcov;

	t.sb = gm.transcripts[p[0]].sb;
	t.eb = gm.transcripts[p.back()].eb;

	// get path
	string path = "";
	for(int i = 0; i < p.size(); i++)
	{
		path += gm.transcripts[p[i]].transcript_id;
		path += "+";
	}
	path.pop_back();

	t.path = path;

	t.fragPath = p;
	//t.scoreList = ps;

	// get exon list
	vector< pair<int32_t, int32_t> > exons;
	build_exonlist(gm, p, exons);

	if(exons.size() == 0) return 0;

	for(int i = 0; i < exons.size(); i++)
	{
		t.add_exon(exons[i].first, exons[i].second);
	}

	// get individual list
	set<int> clist;

	// only add cnum if it is any t in path
    //build_clist(gm, p, clist);
        
	// add cnum more loose, as long as any c has a same j, add it
	//vector<string> curjlist;
	//for(int k = 0; k < exons.size() - 1; k++)
	//{
		//string h = chrm;
		//int32_t q1 = exons[k].second;
        //int32_t q2 = exons[k+1].first;
		//h.append(tostring(q1));
		//h.append(tostring(q2));

		//curjlist.push_back(h);
	//}

	//map<string, set<int> > jtocnum = gm.cjlist;

	//build_clist_j(curjlist, jtocnum, clist);
	//t.clist = clist;

	gm1.add_transcript_b(t);

	return 0;
}

int output_transcripts(ofstream &fout, genome1 &gm, genome1 &gm1, const vector< vector<int> > &p, const vector< vector<double> > &ps)
{
	map<string, int> op;

	for(int i = 0; i < p.size(); i++)
	{
		cout<<"adding transcript "<<i<<endl;

		string gid = gm.transcripts[p[i][0]].gene_id;

		int idx = 1;
		if(op.find(gid) != op.end()) idx = op[gid] + 1;
		op[gid] = idx;

		string tid = gid + "T" + tostring(idx);
		//print_path(i, p[i], gm);
		output_transcript(fout, gm, gm1, p[i], tid, ps[i]);
	}
	return 0;
}

int remove_covered(genome1 &gm, genome1 &gm1)
{
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		//gm1.add_transcript(gm.transcripts[i]);
		//continue;

		bool if_covered = false;
		transcript t1 = gm.transcripts[i];
		vector< pair<int32_t, int32_t> > exons1 = t1.exons;

		for(int j = 0; j < gm.transcripts.size(); j++)
		{
			if(i == j) continue;
			transcript t2 = gm.transcripts[j];

			if(t1.seqname != t2.seqname) continue;
			if(t1.exons.size() > t2.exons.size()) continue;
			if(t1.start < t2.start || t1.end > t2.end) continue;

			// t1 must inside of t2, then check if t1 is totally covered by t2
			vector< pair<int32_t, int32_t> > exons2 = t2.exons;
			int idx1 = 0;
			int idx2 = 0;
			while(idx2 < exons2.size())
			{
				//printf("t1: %d, t2: %d\n", i, j);
				//print_transcript_pair(t1, t2);

				if(exons2[idx2].second < exons1[idx1].first)
				{
					idx2++;
					continue;
				}
				// find overlap start, t1 need to be exactly part of t2
				while(idx1 < exons1.size() && idx2 < exons2.size())
				{
					// if t1 first exon, only need to match in the end
					if(idx1 == 0)
					{
						if(exons1[idx1].second == exons2[idx2].second)
						{
							idx1++;
							idx2++;
							continue;
						}
						else break;
					}

					// if t1 last exon
					if(idx1 == exons1.size() - 1)
					{
						if(exons1[idx1].first == exons2[idx2].first)
						{
							//printf("covered!\n");
							if_covered = true;
						}
						break;
					}

					if(exons1[idx1].first == exons2[idx2].first && exons1[idx1].second == exons2[idx2].second)
					{
						idx1++;
						idx2++;
						continue;
					}
					else break;
				}
				break;
			}
			if(if_covered == true) break;
		}

		if(if_covered == false) gm1.add_transcript(t1);
	}

	return 0;
}

int get_coverage_em(genome1 &gm, genome1 &gm2)
{
	//printf("calculating coverage of merged transcripts...\n");

	// x: merged transcripts; p: transcript fragments
	map<int, set<int> > ptox;
	map<int, set<int> > xtop;

	for(int i = 0; i < gm2.transcripts.size(); i++)
	{
                transcript p = gm2.transcripts[i];
                vector< PI32 > exons1 = p.exons;

		for(int j = 0; j < gm.transcripts.size(); j++)
		{
			transcript x = gm.transcripts[j];
			set<int> xpath(x.fragPath.begin(), x.fragPath.end());

			if(xpath.find(i) != xpath.end()) 
			{
				ptox[i].insert(j);
				xtop[j].insert(i);
				continue;
			}

			if(p.seqname != x.seqname) continue;
			if(p.exons.size() > x.exons.size()) continue;
			//if(p.start < x.start || p.end > x.end) continue;

			if(is_covered(p, x)) 
			{
				ptox[i].insert(j);
				xtop[j].insert(i);
			}
		}
	}

	// init coverage of x
	vector<double> covX;
	vector<double> covX_pre;
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		set<int> pset = xtop[i];
		double temp = 0;
		for(int p : pset)
		{
			temp += gm2.transcripts[p].coverage;
		}
		temp /= pset.size();
		covX.push_back(temp);
	}
	covX_pre = covX;

	// update coverage of x
	int max_em_it_num = 100;
	int it = 0;
	double diff_pre = 0;


	while(it < max_em_it_num)
	{
		double diff = 0;
		for(int i = 0; i < gm.transcripts.size(); i++)
		{
			set<int> pset = xtop[i];
			double temp = 0;
				
			for(int p : pset)
			{
				double psum = 0;
				set<int> xset = ptox[p];
				for(int x : xset) psum += covX_pre[x];
				
				temp += (gm2.transcripts[p].coverage * (covX_pre[i] / psum) ); 
			}
			
			covX[i] = temp;
			//printf("EM it #%d: t %d pset size = %ld, coverage %f -> %f\n", it, i, pset.size(), covX_pre[i], covX[i]);
			
			if(covX[i] >= covX_pre[i]) diff += ((covX[i] - covX_pre[i]) / covX_pre[i]);
			else diff += ((covX_pre[i] - covX[i]) / covX_pre[i]);
		}

		//diff /= gm.transcripts.size();
		covX_pre = covX;
		it++;
		//printf("EM it #%d: diff %f - > %f\n", it, diff_pre, diff);
		if(diff < 1e-2) break;
		if(it != 1 && diff_pre - diff < 1e-4) break;
		diff_pre = diff;
	}

	for(int i = 0; i < gm.transcripts.size(); i++) gm.transcripts[i].coverage = covX[i];

	return 0;
}

int get_coverage_score(genome1 &gm, genome1 &gm2)
{
	set<int> compatible;
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript &t = gm.transcripts[i];
		//printf("compatible size = %ld\n", compatible.size());
		t.coverage = calculate_score_abu(gm2, t, compatible);
		//t.coverage -= calculate_score_extend(gm2, t);

	// get all info used to calculate the coverage, and store them for ML
//		
	//update_compatible(gm2, t, compatible);
	//t.coverage = calculate_score_abu1(gm2, t);
	}

	return 0;
}

int Beaver::write_individual(const string &prefix)
{
	string sgtfPath = prefix + "_sgtf";
	string cmd = "mkdir -p " + sgtfPath;
	const char *command = cmd.c_str();
	int sys = system(command);

	vector<string> ogtf;
	for(int i = 1; i <= cnum; i++)
	{
		ogtf.push_back(sgtfPath + "/" + to_string(i) + ".gtf");
	}

	gm.write_individual(ogtf);

	return 0;
}

int Beaver::write_individual_feature(const string &prefix)
{
	genome1 &meta_gm = gm;
	genome1 &pre_gm = gm2;
	genome1 &ori_gm = gm3;

	//Calculate # of output trsts for each cell
	vector<int> out_tsrt_num(cnum+1, 0);
	for(int i = 0; i < cnum; i++)
	{
		for(int j = 0; j < meta_gm.transcripts.size(); j++)
		{
			transcript &t = meta_gm.transcripts[j];
			if(t.coverage < min_transcript_coverage) continue;
			if(t.clist.find(i+1) == t.clist.end()) continue;
			out_tsrt_num[i+1]++;
		}
	}

	string sgtfPath = prefix + "_sgtf";
	//write transcript features for each cell
	for(int i = 0; i < cnum; i++)
	{
		string fo = sgtfPath + "/" + to_string(i+1) + "_feature.csv";
		ofstream fout(fo.c_str());

		for(int j = 0; j < meta_gm.transcripts.size(); j++)
		{
			transcript &t = meta_gm.transcripts[j];
			if(t.coverage < min_transcript_coverage) continue;
			if(t.clist.find(i+1) == t.clist.end()) continue;

			string tid = t.transcript_id;
			fout << tid.c_str() << ",";
			string seqname = t.seqname; //chromosome
			fout << seqname.c_str() <<",";

			//write basic cell features
			fout << ori_gm.cell_gene_num[i+1] <<","; //input #gene-loci
			fout << ori_gm.cell_trst_num[i+1] <<","; //input #transcripts
			fout << out_tsrt_num[i+1] << ","; //output #transcripts

			//write cell-specific transcript features
			vector<double> specific_features;
			calculate_cell_specific_features(ori_gm, pre_gm, t, i+1, specific_features);
			for(int k = 0; k < specific_features.size(); k++)
			{
				if(k != specific_features.size() - 1) fout<<specific_features[k]<<",";
				else fout<<specific_features[k]<<endl;
			}

			//write meta transcript features
			/*const vector<double> &f = features[j];
			for(int k = 0; k < f.size(); k++)
			{
				if(k != f.size() - 1) fout<<f[k]<<",";
				else fout<<f[k]<<endl;
			}*/
		}

	}
	return 0;
}

int Beaver::calculate_cell_specific_features(genome1 &ori_gm, genome1 &pre_gm, transcript &t, int cid, vector<double> & specific_features)
{
	printf("\nCalculating specific features for %s", t.transcript_id.c_str());
	specific_features.clear();

	//how many junction supported by cell cid
	specific_features.push_back(t.cjcount[cid]);
	specific_features.push_back(1.0*t.cjcount[cid]/(t.exons.size()-1));
	specific_features.push_back(t.exons.size()-1);

	//Calculate cell-specific compatible junction coverage, length = #junctions
	for(int i = 0; i < t.exons.size() - 1; i++)
	{
		string h = t.seqname;
		h.append(tostring(t.exons[i].second));
		h.append(tostring(t.exons[i+1].first));

		set<int> &tlist = pre_gm.tjlist[h];
		set<int> &cjlist = ori_gm.cjlist[h];

		double curScore = 0;

		printf("\ncjlist: ");
		for(auto c: cjlist) printf("%d ", c);
		printf("\n");

		if(cjlist.find(cid) != cjlist.end())
		{
			for(int tidx : tlist)
			{
				transcript & supt = pre_gm.transcripts[tidx];

				//check whether the query trst is support the cell AND compatible with current trst
				if(supt.clist.find(cid) != supt.clist.end())
				{
					if(is_covered(supt, t))
					{
						curScore += supt.coverage;
					}
				}
			}
		}
		specific_features.push_back(curScore);
	}

	//Calculate cell-specific junction coverage regardless compatibility, length = #junctions
	for(int i = 0; i < t.exons.size() - 1; i++)
	{
		string h = t.seqname;
		h.append(tostring(t.exons[i].second));
		h.append(tostring(t.exons[i+1].first));

		set<int> &tlist = pre_gm.tjlist[h];
		set<int> &cjlist = ori_gm.cjlist[h];

		if(cjlist.find(cid) == cjlist.end())
		{
			specific_features.push_back(0);
			continue;
		}

		double curScore = 0;
		for(int tidx : tlist)
		{
			transcript & supt = pre_gm.transcripts[tidx];
			printf("\nFor junction %d, Check trst %s(cov=%.4f)'s cell support: ", i+1, supt.transcript_id.c_str(), supt.coverage);
			for(auto c : supt.clist) printf("%d ", c);
			printf("\n");

			//check whether the query trst is support the cell AND compatible with current trst
			if(supt.clist.find(cid) != supt.clist.end())
			{
				curScore += supt.coverage;
			}
		}
		specific_features.push_back(curScore);
	}
	return 0;
}

int add_tjlist(genome1 &gm)
{
	map<string, set<int> > &tjlist = gm.tjlist;
	tjlist.clear();
	for(int k = 0; k < gm.transcripts.size(); k++)
	{
		vector<PI32> &exons = gm.transcripts[k].exons;

		for(int i = 0; i < exons.size() - 1; i++)
        	{
				string h = gm.transcripts[k].seqname;
				h.append(tostring(exons[i].second));
				h.append(tostring(exons[i+1].first));

				tjlist[h].insert(k);
		}
	}

	return 0;
}

int add_jtocnum(genome1 &gm)
{
	map<string, set<int> > &jtocnum = gm.jtocnum;
	for(int k = 0; k < gm.transcripts.size(); k++)
	{
		vector<PI32> &exons = gm.transcripts[k].exons;

		for(int i = 0; i < exons.size() - 1; i++)
		{
			string h = gm.transcripts[k].seqname;
			h.append(tostring(exons[i].second));
			h.append(tostring(exons[i+1].first));

			jtocnum[h].insert(gm.transcripts[k].clist.begin(), gm.transcripts[k].clist.end());
		}
	}

	return 0;
}

int add_telist(genome1 &gm)
{
        map<string, set<int> > &telist = gm.telist;
        for(int k = 0; k < gm.transcripts.size(); k++)
        {
                vector<PI32> &exons = gm.transcripts[k].exons;

                for(int i = 1; i < exons.size() - 1; i++)
                {
                        string h = gm.transcripts[k].seqname;
                        h.append(tostring(exons[i].first));
                        h.append(tostring(exons[i].second));

                        telist[h].insert(k);
                }
        }

        return 0;
}

int add_extendable_to_empty(genome1 &gm)
{
	map<string, set<int> > &tjlist = gm.tjlist;
        
	for(int k = 0; k < gm.transcripts.size(); k++)
	{
		transcript &curt = gm.transcripts[k];
		vector<PI32> &exons = curt.exons;

		set<int> visited;

		for(int i = 0; i < exons.size() - 1; i++)
		{
			string h = gm.transcripts[k].seqname;
			h.append(tostring(exons[i].second));
			h.append(tostring(exons[i+1].first));

			set<int> &tlist = gm.tjlist[h];

			for(int tidx : tlist)
			{
				if(visited.find(tidx) != visited.end()) continue;
				visited.insert(tidx);

				transcript &t = gm.transcripts[tidx];
				// skip if t is inside of curt
				if(t.exons[0].second >= exons[0].second && t.exons.back().first <= exons.back().first) continue;
				if(t.exons[0].second <= exons[0].second && t.exons.back().first >= exons.back().first)
				{
					if(!is_covered(curt, t)) continue;
					if(t.exons[0].second < exons[0].second) t.empty_vote[0][1] += t.coverage;
					if(t.exons.back().first > exons.back().first) t.empty_vote[1][1] += t.coverage;
					continue;
				}
				/*
				// curt is inside of t
				if(t.exons[0].second <= exons[0].second && t.exons.back().first >= exons.back().first)
				{
					if(!is_covered(curt, t)) continue;
					if(t.exons[0].second < exons[0].second) t.empty_vote[0][1] += t.coverage;
					if(t.exons.back().first > exons.back().first) t.empty_vote[1][1] += t.coverage;
					continue;
				}

				// t can extend to curt
                        	if(t.exons[0].second < exons[0].second)
                        	{
                                	if(is_extendable(t, curt)) t.empty_vote[0][1] += t.coverage;
					continue;
                        	}

				// curt can extend to t
                        	if(t.exons.back().first > exons.back().first)
                        	{
                                	if(is_extendable(curt, t)) t.empty_vote[1][1] += t.coverage;
					continue;
                        	}
				*/

			}
		}
        }

	return 0;
}

int addback_ft_coverage(genome1 &gm, vector<transcript> &ft)
{
	for(int i = 0; i < ft.size(); i++)
	{
		string s = gm.compute_intron_hashing(ft[i]);
		if(gm.intron_hashing.find(s) != gm.intron_hashing.end())
		{
			transcript &t = gm.transcripts[gm.intron_hashing[s]];
			t.coverage += ft[i].coverage;
		}
	}

	return 0;
}

int vote_empty(genome1 &gm)
{
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript &t = gm.transcripts[i];
		vector< vector<double> > &ev = t.empty_vote;

		double threshold = 1;
		
		if(ev[0][0] >= threshold * ev[0][1]) t.left_empty = false;
		else t.left_empty = true;

		if(ev[1][0] >= threshold * ev[1][1]) t.right_empty = false;
		else t.right_empty = true;

		/*
		//printf("# %d t: (%f, %f), (%f, %f)\n", i, ev[0][0], ev[0][1], ev[1][0], ev[1][1]);
		if(ev[0][0] > 0 && ev[0][1] > 0) 
		{
			printf("left!!! # %d t: (%f, %f), (%f, %f)\n", i, ev[0][0], ev[0][1], ev[1][0], ev[1][1]);
		}
		if(ev[1][0] > 0 && ev[1][1] > 0)
                {
                        printf("right!!! # %d t: (%f, %f), (%f, %f)\n", i, ev[0][0], ev[0][1], ev[1][0], ev[1][1]);
                }
		*/
	}
	return 0;
}

int store_jb(genome1 &gm)
{
	map<string, vector<int32_t>> &jb = gm.jblist;

	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript &t = gm.transcripts[i];
		if(jb.find(t.seqname) == jb.end()) jb[t.seqname] = {};

		for(int j = 0; j < t.exons.size() - 1; j++)
		{
			jb[t.seqname].push_back(t.exons[j].second);
			jb[t.seqname].push_back(t.exons[j+1].first);
		}
	}

	for(auto c : jb)
	{
		vector<int32_t> &list = c.second;
		sort(list.begin(), list.end());
	}

	return 0;
}

int regroup_transcripts(genome1 &gm, genome1 &gm1)
{
	// store all junc boundaries by chr
	store_jb(gm);
	map<string, vector<int32_t>> &jb = gm.jblist;

	// regroup t
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript &t = gm.transcripts[i];
		vector<transcript> &tlist = t.tlist;
		
		for(transcript &curt : tlist)
		{
			int32_t s = find_interval(curt.start, jb[curt.seqname]);
			int32_t e = find_interval(curt.end, jb[curt.seqname]);
			curt.sb = s;
			curt.eb = e;
			gm1.add_transcript_b(curt);
		}
	}

	return 0;
}

int store_features(genome1 &meta_gm, genome1 &gm, vector<vector<double>> &features)
{
	features.clear();

	for(int i = 0; i < meta_gm.transcripts.size(); i++)
    {
		vector<double> cur_feature;

        transcript &curt = meta_gm.transcripts[i];
		vector<PI32> &exons = curt.exons;

		// 1st, 2nd feature, bottleneck coverage and higest coverage of fragments
		cur_feature.push_back(curt.coverage);
		cur_feature.push_back(curt.hcov);

		// 3rd feature: extendable score; field[0]
		set<int> visited;
		double score = 0;

		for(int i = 0; i < exons.size() - 1; i++)
		{
			string h = curt.seqname;
			h.append(tostring(exons[i].second));
			h.append(tostring(exons[i+1].first));

			set<int> &tlist = gm.tjlist[h];
			for(int tidx : tlist)
			{
				if(visited.find(tidx) != visited.end()) continue;
				visited.insert(tidx);

				transcript &t = gm.transcripts[tidx];

				// skip if t is inside of curt
				if(t.exons[0].second >= exons[0].second && t.exons.back().first <= exons.back().first) continue;
				
				if(t.exons[0].second <= exons[0].second && t.exons.back().first >= exons.back().first)
				{
					if(is_covered(curt, t)) score += t.coverage;
					continue;
				}

				// t can extend to curt
				if(t.exons[0].second < exons[0].second)
				{
					if(is_extendable(t, curt)) score += t.coverage;
					continue;
				}

				// curt can extend to t
				if(t.exons.back().first > exons.back().first)
				{
					if(is_extendable(curt, t)) score += t.coverage;
					continue;
				}
			}
		}
		cur_feature.push_back(score);

		// 2nd feature: # junctions;field[1]
		cur_feature.push_back(exons.size() - 1);

		// other features: score for each junction; field[2, #junction+1]
		for(int i = 0; i < exons.size() - 1; i++)
		{
			string h = curt.seqname;
			h.append(tostring(exons[i].second));
			h.append(tostring(exons[i+1].first));

			set<int> &tlist = gm.tjlist[h];

			double curScore = 0;

			for(int tidx : tlist)
			{
				
				// print exon lists of two t
				/*	
				for(int m = 0; m < exons.size(); m++)
				{
					cout<<"("<<exons[m].first<<","<<exons[m].second<<") ";
				}
				printf("\n");
				vector<PI32> &exons1 = gm.transcripts[tidx].exons;
				for(int n = 0; n < exons1.size(); n++)
				{
						cout<<"("<<exons1[n].first<<","<<exons1[n].second<<") ";
				}
				printf("\n");
				*/
				
				//if((gm.transcripts[tidx].exons[0].second == curt.exons[0].second) && (gm.transcripts[tidx].sb != curt.sb)) continue;
				//if((gm.transcripts[tidx].exons.back().first == curt.exons.back().first) && (gm.transcripts[tidx].eb != curt.eb)) continue;

				if(is_covered(gm.transcripts[tidx], curt))
				{
					//printf("covered!\n");
					/*
					if(gm.transcripts[tidx].exons[0].second == curt.exons[0].second)
					{
							if(find_interval(gm.transcripts[tidx].start, gm.jblist[curt.seqname]) != find_interval(gm.transcripts[tidx].start, gm.jblist[curt.seqname])) continue;
					}
					if(gm.transcripts[tidx].exons.back().first ==  curt.exons.back().first)
					{
							if(find_interval(gm.transcripts[tidx].end, gm.jblist[curt.seqname]) != find_interval(gm.transcripts[tidx].end, gm.jblist[curt.seqname])) continue;
					}
					*/

					curScore += gm.transcripts[tidx].coverage;
				}
				//printf("curScore = %f\n", curScore);

			}
			cur_feature.push_back(curScore);
		}

		// support coverage for each junction regardless compatibility; field[#junction+2, 2*#junction+1]
		for(int i = 0; i < exons.size() - 1; i++)
		{
			string h = curt.seqname;
			h.append(tostring(exons[i].second));
			h.append(tostring(exons[i+1].first));

			set<int> &tlist = gm.tjlist[h];

			double curAllScore = 0;

			for(int tidx : tlist)
			{
				transcript &t = gm.transcripts[tidx];
				curAllScore += t.coverage;
			}

			cur_feature.push_back(curAllScore);
		}

		// how junction supported by cells; field[2*#junction+2, 2*#junction+2+2*#cell_support]
		map<int,int> &cur_cjcount = curt.cjcount;
		
		cur_feature.push_back(cur_cjcount.size());
		for (const auto& pair : cur_cjcount) {
			cur_feature.push_back(pair.first);  // Push the cell id
			cur_feature.push_back(pair.second); // Push the overlapping j count
		}

		// how meta-trst supported by fragments
		printf("Meta transcript %s consist of %ld fragments: \n", curt.transcript_id.c_str(), curt.fragPath.size());
		if(curt.fragPath.size() == 0) //additional trst from aletsch
		{
			cur_feature.push_back(1);
			cur_feature.push_back(curt.coverage);
		}
		else
		{
			cur_feature.push_back(curt.fragPath.size());
			for(int j = 0; j < curt.fragPath.size(); j++)
			{
				transcript &frag = gm.transcripts[curt.fragPath[j]];
				printf("%d -- %s.\n", j+1, frag.transcript_id.c_str());
				cur_feature.push_back(frag.coverage);
			}
		}

		features.push_back(cur_feature);
	}
	return 0;
}

int write_meta_features(genome1 &gm, vector<vector<double>> &features, const string &prefix)
{
	string fo = prefix + "_feature.csv";
	ofstream fout(fo.c_str());

	for(int i = 0; i < features.size(); i++)
	{
		vector<double> &f = features[i];
		transcript &t = gm.transcripts[i];

		string tid = t.transcript_id;
		fout<<tid.c_str()<<",";
		string seqname = t.seqname;
		fout<<seqname.c_str()<<",";

		for(int k = 0; k < f.size(); k++)
		{
			if(k != f.size() - 1) fout<<f[k]<<",";
			else fout<<f[k]<<endl;
		}
	}

	return 0;
}

int addback_all_individuals(genome1& gm1, genome1& ori_gm)
{
	for(int i = 0; i < ori_gm.transcripts.size(); i++)
	{
		transcript &t = ori_gm.transcripts[i];
		if(t.exons.size() > 1 && t.coverage > confident_input_coverage) //high confident
		{
			gm1.add_transcript_b(t);
		}
	}
	return 0;
}

int Beaver::link_merge(const string &prefix)
{
	printf("# transcripts loaded by Beaver = %ld\n", gm.transcripts.size());
	if(gm.transcripts.size() == 0) return 0;

	/*
	// remove covered transcripts, store removed set in gm1
	gm1.clear();
	remove_covered(gm, gm1);
	gm = gm1;
	gm1.clear();
	printf("# transcripts after removing covered = %ld\n", gm.transcripts.size());

	remove_ir(gm, gm1);
	gm = gm1;
	gm1.clear();
	*/

	//add_transcript_cw(gm);
	//add_transcript_cwe(gm);
	//printf("constructing junction-node list...\n");
	//add_tjlist(gm);
	//add_telist(gm);
	//printf("voting empty...\n");
	//share_empty_info(gm);
	//add_extendable_to_empty(gm),
	//vote_empty(gm);	


	//string fo = prefix + "_all.gtf";
	//gm.write(fo);

	gm3 = gm; //original genome
	add_jtocnum(gm3);
	store_jb(gm3);
	// store all junction boundaries in each chr, then group t
	regroup_transcripts(gm, gm1);
	gm = gm1;
	gm1.clear();
	add_tjlist(gm);
	add_tjlist(gm3);
	printf("# transcripts after regroup = %ld\n", gm.transcripts.size());

	printf("adding edges...\n");
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		for(int j = i + 1; j < gm.transcripts.size(); j++)
		{
			transcript &t1 = gm.transcripts[i];
			transcript &t2 = gm.transcripts[j];
				
			if(t1.seqname != t2.seqname) continue;
			if(t1.strand != t2.strand) continue;

			if(t1.exons[0].first > t2.exons.back().first && t1.exons[0].first < t2.exons.back().second && t1.exons[0].second > t2.exons.back().second) 
			{
				if(t1.exons[0].second - t2.exons.back().first > 200) continue;
			}
			if(t2.exons[0].first > t1.exons.back().first && t2.exons[0].first < t1.exons.back().second && t2.exons[0].second > t1.exons.back().second)
			{
				if(t2.exons[0].second - t1.exons.back().first > 200) continue;
			}
		
			if(t1.end < t2.start)
			{
				if(t1.end + max_gap_bp < t2.start) continue;
				else 
				{
					int w = get_edge_weight(t1, t2, 0);
					t1.so.insert(j);
					t1.sow[j] = w;
					t2.si.insert(i);
					t2.siw[i] = w;
				}
			}
			if(t1.start > t2.end) 
			{
				if(t1.start > t2.end + max_gap_bp) continue;
				else
				{
					int w = get_edge_weight(t1, t2, 0);
					t2.so.insert(i);
					t2.sow[i] = w;
					t1.si.insert(j);
					t1.siw[j] = w;

				}
			}

			int if_link = 0;
			
			if(t1.start <= t2.start) if_link = add_in_and_out_edge(t1, t2, i, j);
			else if_link = add_in_and_out_edge(t2, t1, j, i);
			//if(t1.exons[0].second <= t2.exons[0].second) if_link = add_in_and_out_edge(t1, t2, i, j);
			//else if_link = add_in_and_out_edge(t2, t1, j, i);
		}
	}

	/*
	for(int i = 0; i < gm.transcripts.size(); i++)
        {
			print_transcript(gm.transcripts[i], i);
        }
	*/

	//store_ft(gm, ft);
	gm2 = gm;

	assign_gid(gm, subgraphs);

	//print_graph(gm, prefix);
	
	rm_duplicate_subgraph(gm, subgraphs, sg_covered, sgj);

	printf("finding paths...\n");
	
	//find_allpaths_greedy(gm, subgraphs, paths, sg_covered, sgj); // dp + greedy
	
	//find_allpaths_greedy_extend(gm, subgraphs, paths, pscores, sg_covered); // greedy extend from node with max score

	//find_allpaths_greedy_extend_abuj(gm, subgraphs, paths, pscores, sg_covered); // greedy extend from all nodes, dp store top 3 faths, final score = score*#junctions

	find_allpaths_dp_abuj(gm, subgraphs, paths, pscores, sg_covered);

	//printf("printing all intput transcripts...\n");	
	//string fo = prefix + "_all.gtf";
        //gm.write(fo);

	//find_allpaths(gm, paths); // junction weight approach
	
	// one edge - one junction
	//allow_path_one_junction(gm, paths);

	/*
	printf("# paths = %ld\n", paths.size());
	for(int i = 0; i < paths.size(); i++)
	{
		print_path(i, paths[i], gm);
	}
	*/
	

	string ofile = prefix + ".gtf";
	ofstream fout(ofile.c_str());
        //if(fout.fail()) return 0;
	gm1.clear();
	// store merged transcripts in gm1
	printf("paths to transcripts....\n");
	output_transcripts(fout, gm, gm1, paths, pscores);
	//gm1.write(ofile);
	//fout.close();
	
	printf("# transcripts from individual: %ld\n", gm2.transcripts.size());	
	// add back all individuals to the Beaverd set
	addback_all_individuals(gm1, gm2);
	printf("# transcripts after adding back all individuals: %ld\n", gm1.transcripts.size());
	
	// store remove covered transcripts in gm
	gm.clear();
	//addback_ft(gm1, ft);
	//printf("remove covered...\n");
	//remove_covered(gm1, gm);
	gm = gm1;
	
	printf("get coverage using final score...\n");
	// set coverage of merged transcripts
	//get_coverage_em(gm, gm2);
	
	get_coverage_score(gm, gm2);
	//get_coverage_score(gm, gm3);

	// get clist for each transcript
	build_clist_j(gm, gm3);

	//addback_ft_coverage(gm, ft);

	printf("writing meta-level gtf...\n");
	gm.write(ofile);

	store_features(gm, gm2, features);
	printf("writing features...\n");
	//store_features(gm, gm3, features);
	write_meta_features(gm, features, prefix);

	// write individual gtf
	
	/*
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		print_transcript(gm.transcripts[i], i);
	}
	*/
	write_individual(prefix);
	//gm -- meta_gm; gm2 -- original trsts with compatible edges; gm3 -- original trsts
	write_individual_feature(prefix);
	
	return 0;
}

int Beaver::rm_cover()
{

	remove_covered(gm, gm1);
	return 0;
}

int Beaver::meta_vote(const string &prefix)
{
	printf("# transcripts loaded by meta-vote = %ld\n", gm.transcripts.size());

        string ofile = prefix;

        gm1.clear();

        for(int i = 0; i < gm.transcripts.size(); i++)
        {
                transcript &t1 = gm.transcripts[i];
                if( t1.clist.size() >= 3)
                {
                        //printf("t.clist size = %ld\n", t1.clist.size());
                        gm1.add_transcript(t1);
                }
        }

        gm1.write(ofile);

        return 0;
}

int Beaver::group_vote(const string &glist, const string &prefix)
{
	build_union1(glist);
	string ofile = prefix;

        gm1.clear();

        for(int i = 0; i < gm.transcripts.size(); i++)
        {
                transcript &t1 = gm.transcripts[i];

                if(t1.clist.size() >= 3)
                {
                        //printf("t.clist size = %ld\n", t1.clist.size());
                        gm1.add_transcript(t1);
                }
        }
        
	gm1.write(ofile);
	return 0;
}

void print_vector(vector<vector<double> > two_D_vector)
{
    for (auto it : two_D_vector) {
        for (auto ij : it) {
            cout << ij << " ";
        }
        cout << endl;
    }
}


int Beaver::filter_junction(const string &fi, const string &fo)
{
	load_genome_all(fi.c_str(), &gm, 0);

	map< string, vector<int> > genes; // transcript idx list in each gene_id;
	for(int i = 0; i < gm.transcripts.size(); i++) genes[gm.transcripts[i].gene_id].push_back(i);

	for(auto g : genes)
	{
		vector<int> &tlist = g.second;
		
		vector< vector<double> > tclist; // ((tcov, tidx))
		for(int i : tlist) tclist.push_back({gm.transcripts[i].coverage, (double)i});
		
		sort(tclist.begin(), tclist.end());
		reverse(tclist.begin(), tclist.end());
		//print_vector(tclist);

		vector<int> tl;
		for(auto i : tclist) tl.push_back((int)i[1]);

		set<PI32> jlist;
		for(int i : tl)
		{
			transcript &t = gm.transcripts[i];
			
			vector<PI32> curj;
			for(int j = 0; j < t.exons.size() - 1; j++)
			{
				curj.push_back({t.exons[j].second, t.exons[j+1].first});
			}

			bool if_new_junction = false;
			
			for(auto junc : curj)
			{
				if(jlist.find(junc) == jlist.end())
				{
					if_new_junction = true;
					break;
				}
			}

			if(if_new_junction) 
			{
				for(auto junc : curj) jlist.insert(junc);
				gm1.add_transcript(t);
			}
			else printf("transcript %d filtered...\n", i);
		}
	}
	
	gm1.write(fo);

	return 0;
}

bool if_left_chain_match(vector<PI32> exons1, vector<PI32> exons2)
{
	bool if_match = false;
	
	if(exons2.size() > exons1.size()) return if_match;

	int i = 0;

	for(int j = 0; j < exons2.size(); j++)
	{
		if(j == 0)
		{
			if(exons1[i].second != exons2[j].second) break;
		}	
		else if(j > 0 && j < exons2.size() - 1)
		{
			if(exons1[i].first != exons2[j].first || exons1[i].second != exons2[j].second) break;
		}
		else if(j == exons2.size() - 1)
		{
			if(exons1[i].first == exons2[j].first) if_match = true; 
		}
		
		i++;
	}

	return if_match;
}

bool if_right_chain_match(vector<PI32> exons1, vector<PI32> exons2)
{
        bool if_match = false;

	if(exons2.size() > exons1.size()) return if_match;

        int i = exons1.size() - 1;

        for(int j = exons2.size() - 1; j >= 0; j--)
        {
                if(j == 0)
                {
                        if(exons1[i].second == exons2[j].second) if_match = true;
                }
                else if(j > 0 && j < exons2.size() - 1)
                {
                        if(exons1[i].first != exons2[j].first || exons1[i].second != exons2[j].second) break;
                }
                else if(j == exons2.size() - 1)
                {
                        if(exons1[i].first != exons2[j].first) break;
                }
		
		i--;
        }

        return if_match;
}


int Beaver::split_empty(const string &fi)
{
	gm1 = gm;
	gm.clear();
	load_genome_all(fi.c_str(), &gm, 0);

	vector<bool> if_share_empty(gm.transcripts.size(), false);

	vector<string> empty_side_list;
	for(int i = 0; i < gm1.transcripts.size(); i++) empty_side_list.push_back(empty_side(gm1.transcripts[i]));

	// gm: full-length; gm1: non-full-length
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		transcript &ft = gm.transcripts[i];

		for(int j = 0; j < gm1.transcripts.size(); j++)
		{
			transcript &nft = gm1.transcripts[j];

			if(ft.seqname != nft.seqname) continue;
			if(ft.strand != nft.strand) continue;

			if(ft.exons[0].second == nft.exons[0].second)
			{
				if(empty_side_list[j] == "left" || empty_side_list[j] == "both")
				{
					if(if_left_chain_match(nft.exons, ft.exons))
					{
						if_share_empty[i] = true;
						break;
					}
				}
			}

			if(ft.exons[ft.exons.size()-1].first == nft.exons[nft.exons.size()-1].first)
			{
				if(empty_side_list[j] == "right" || empty_side_list[j] == "both")
				{
					if(if_right_chain_match(nft.exons, ft.exons))
					{
						if_share_empty[i] = true;
                                                break;
					}
				}
			}
		}
	}

	string fo1 = fi + ".e"; // full-length transcripts may contain non-full point
	string fo2 = fi + ".f"; // full-length transcripts still full

	gm1.clear();
	gm2.clear();
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		if(if_share_empty[i]) gm1.add_transcript(gm.transcripts[i]);
		else gm2.add_transcript(gm.transcripts[i]);

	}
	gm1.write(fo1);
	gm2.write(fo2);

	return 0;
}

int Beaver::split_support_cell(const string &prefix)
{
	map<int, vector<int> > X;

	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		X[(int)gm.transcripts[i].clist.size()].push_back(i);
	}

	for(auto x : X)
	{
		string fo = prefix + "." + to_string(x.first) + ".cellsupport.gtf";
		
		for(int i : x.second) gm1.add_transcript(gm.transcripts[i]);

		gm1.write(fo);
		gm1.clear();

	}

	return 0;
}

int Beaver::split_match(const string &fi)
{
	load_genome_all(fi.c_str(), &gm1, 0);

	vector<int> match_set;
	vector<int> unmatch_set;

	//printf("# full = %d/%d, # non-full = %d/%d\n", gm.transcripts.size(), gm.intron_hashing.size(), gm1.transcripts.size(), gm1.intron_hashing.size());
	for(auto m : gm1.intron_hashing)
	{
		if(gm.intron_hashing.find(m.first) != gm.intron_hashing.end())
		{
			match_set.push_back(gm1.intron_hashing[m.first]);
		}
		else unmatch_set.push_back(gm1.intron_hashing[m.first]);
	}

	sort(match_set.begin(), match_set.end());
	sort(unmatch_set.begin(), unmatch_set.end());

	string fo = fi + ".f";
	for(int i = 0; i < match_set.size(); i++) gm2.add_transcript(gm1.transcripts[match_set[i]]);
	gm2.write(fo);
	
        gm2.clear();
	fo = fi + ".nf";
	for(int i = 0; i < unmatch_set.size(); i++) gm2.add_transcript(gm1.transcripts[unmatch_set[i]]);
	gm2.write(fo);
	
	return 0;
}

void build_blist(genome1 &gm, vector<set<int32_t>> &blist)
{

	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		vector<pair<int32_t, int32_t>> &exons = gm.transcripts[i].exons;
		set<int32_t> temp;
		
		if(exons.size() == 1)
		{
			blist.push_back(temp);
			continue;
		}

		for(int j = 0; j < exons.size(); j++)
		{
			if(j != 0) temp.insert(exons[j].first);
			if(j != exons.size()-1) temp.insert(exons[j].second);
		}

		blist.push_back(temp);
	}
}

bool is_boundary_covered(set<int32_t> &bl1, set<int32_t> &bl2)
{
	bool if_covered = true;

	for(int32_t b : bl1)
	{
		if(bl2.find(b) == bl2.end()) 
		{
			if_covered = false;
			break;
		}
	}

	return if_covered;
}

int Beaver::filter_boundary(const string &fi, const string &fo)
{
	load_genome_all(fi.c_str(), &gm, 0);

	vector<set<int32_t>> blist;
	build_blist(gm, blist);

	for(int i = 0; i < gm.transcripts.size(); i++)
        {
		bool add = true;
		transcript &t1 = gm.transcripts[i];

		for(int j = 0; j < gm.transcripts.size(); j++)
		{
			transcript &t2 = gm.transcripts[j];

			if(i == j) continue;
			if(t1.seqname != t2.seqname) continue;
			if(t1.end < t2.start) continue;
			if(t1.start > t2.end) continue;
			if(blist[i].size() >= blist[j].size()) continue;

			if(is_boundary_covered(blist[i], blist[j])) 
			{
				add = false;
				break;
			}
		}

		if(add) gm1.add_transcript(t1);
	}

	gm1.write(fo);

	return 0;
}

int Beaver::write_union(const string &fo)
{
	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		gm.transcripts[i].clist = {};
	}
	gm.write(fo);
	return 0;
}

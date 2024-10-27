#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

double min_transcript_coverage = -10000;
bool merge_coverage_as_counts = false;
bool merge_coverage_log = false;
int num_threads = 10;
int max_shift_bp = 2;
int max_gap_bp = 0;
int max_unmatch_exon = 3;//3;
float max_unmatch_exon_perc = 2;
int max_path = 2;
double min_start_score = 2;
int num_path_per_node = 15;
int num_path_per_graph = 100;
double min_input_transcript_coverage = 0;
double confident_input_coverage = 0.5;

int parse_parameters(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		if(string(argv[i]) == "-c")
		{
			min_transcript_coverage = atof(argv[i + 1]);
			i++;
		}
		if(string(argv[i]) == "-ic")
		{
				min_input_transcript_coverage = atof(argv[i + 1]);
				i++;
		}
		else if(string(argv[i]) == "-t")
		{
			num_threads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-n")
		{
			merge_coverage_as_counts = true;
		}
		else if(string(argv[i]) == "-e")
		{
			merge_coverage_log = true;
		}
		else if(string(argv[i]) == "-pn")
		{
			num_path_per_node = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-pg")
		{
			num_path_per_graph = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-sbp")
		{
			max_shift_bp = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-mue")
		{
			max_unmatch_exon = atoi(argv[i + 1]);
			i++;
		} 
	}

	return 0;
}

int print_help()
{
	printf("\n");
	printf("usage: scAletsch <input-gtf-list> <output-prefix> [options]\n");
	printf("\n");
	printf("options:\n");
	printf(" %-14s  %s\n", "-t <integer>",  "number of threads");
	return 0;	
}

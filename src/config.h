#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

extern double min_transcript_coverage;
extern bool merge_coverage_as_counts;
extern bool merge_coverage_log;
extern int num_threads;
extern int max_shift_bp;
extern int max_gap_bp;
extern int max_unmatch_exon;
extern float max_unmatch_exon_perc;
extern int max_path;
extern double min_start_score;
extern int num_path_per_node;
extern int num_path_per_graph;
extern double min_input_transcript_coverage;
extern double confident_input_coverage;


int parse_parameters(int argc, const char ** argv);
int print_help();

#endif

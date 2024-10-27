#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "scAletsch.h"
#include "config.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc == 1)
	{
		print_help();
		return 0;
	}

	parse_parameters(argc, argv);

        assert(argc >= 3);
        scAletsch prog;
        prog.build_union1(argv[1]);
        prog.link_merge(argv[2]);

    return 0;
}

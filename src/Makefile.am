bin_PROGRAMS=scAletsch

scAletsch_CPPFLAGS = -pthread -I$(top_srcdir)/lib/gtf -std=c++11
scAletsch_LDFLAGS = -L$(top_builddir)/lib/gtf -pthread
scAletsch_LDADD = -lgtf

scAletsch_SOURCES = genome1.h genome1.cc scAletsch.h scAletsch.cc config.h config.cc main.cc

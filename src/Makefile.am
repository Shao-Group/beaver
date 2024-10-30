bin_PROGRAMS=Beaver

Beaver_CPPFLAGS = -pthread -I$(top_srcdir)/lib/gtf -std=c++11
Beaver_LDFLAGS = -L$(top_builddir)/lib/gtf -pthread
Beaver_LDADD = -lgtf

Beaver_SOURCES = genome1.h genome1.cc beaver.h beaver.cc config.h config.cc main.cc

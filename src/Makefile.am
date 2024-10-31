bin_PROGRAMS=beaver

beaver_CPPFLAGS = -pthread -I$(top_srcdir)/lib/gtf -std=c++11
beaver_LDFLAGS = -L$(top_builddir)/lib/gtf -pthread
beaver_LDADD = -lgtf

beaver_SOURCES = genome1.h genome1.cc beaver.h beaver.cc config.h config.cc main.cc

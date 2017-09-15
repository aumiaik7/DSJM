#DSJM_FLAGS = $(EXTRA_FLAGS)  -O3
#AM_CPPFLAGS =  $(DSJM_FLAGS) 
#AM_LDFLAGS = $(DSJM_FLAGS)
#AM_CXXFLAGS= -O3
CORE_SOURCES = mmio.cpp Matrix.cpp  CLI.cpp Utility.cc Timer.cpp  IMatrix.cc  Configuration.cc  IRowColumnDS.cc 
INCLUDEdir = -I$(top_srcdir)/src
INCLUDE_HEADERS = BucketPQ.hh CLI.h Definitions.h Configuration.hh IMatrix.hh IRowColumnDS.hh Matrix.hh mmio.h NNZTag.hh OrderingMethod.hh detail/PriorityQueue.hh Result.hh RLF.hh Timer.h Utility.h
## Profiling Flags
## CPP_PROFILING_FLAGS = -g3 -O0 -pg
## CPP_PROFILING_FLAGS = -g -O0
## CPP_PROFILING_FLAGS = -O0

## Executables to Build.
noinst_LIBRARIES = libmatrix.a
libmatrix_a_SOURCES = $(CORE_SOURCES)





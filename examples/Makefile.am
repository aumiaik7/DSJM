#DSJM_FLAGS = $(EXTRA_FLAGS)  -O3
#AM_CPPFLAGS =  $(DSJM_FLAGS) 
#AM_LDFLAGS = $(DSJM_FLAGS)
#AM_CXXFLAGS= -O3
CORE_SOURCE_FILES = $(top_srcdir)/src/mmio.cpp $(top_srcdir)/src/Matrix.cpp  $(top_srcdir)/src/CLI.cpp $(top_srcdir)/src/Utility.cc $(top_srcdir)/src/Timer.cpp  $(top_srcdir)/src/IMatrix.cc  $(top_srcdir)/src/Configuration.cc  $(top_srcdir)/src/IRowColumnDS.cc 
INCLUDEdir = -I$(top_srcdir)/example
INCLUDE_HEADERS = driver.hh Eye.hh FD.hh gcolor.hh StructuredJacobian.hh
# Include Directory for compilation
INCLUDES 	= -I$(top_srcdir)/src
#LIB = $(top_srcdir)/libs/libcsparse.a

FORWARDDIFFERENCE_SOURCE_FILES = ForwardDifference.cc  StructuredJacobian.cc

## Programs to Build
bin_PROGRAMS = $(FORWARDDIFFERENCE) $(VERIFY) $(GCOLOR)

## Mappings of the executables. / Dictionary-Like
FORWARDDIFFERENCE = ForwardDifference
VERIFY = dsjmverify

ForwardDifference_SOURCES =$(FORWARDDIFFERENCE_SOURCE_FILES) $(CORE_SOURCE_FILES)
ForwardDifference_CXXFLAGS = $(CPP_PROFILING_FLAGS)
#ForwardDifference_LDADD = $(LIB)

dsjmverify_SOURCES = dsjmverify.cpp

GCOLOR = gcolor
gcolor_SOURCES = gcolor.cpp driver.cc $(CORE_SOURCE_FILES)
gcolor_CXXFLAGS = $(CPP_PROFILING_FLAGS)

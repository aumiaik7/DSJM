#DSJM_FLAGS = $(EXTRA_FLAGS)  -O3
#AM_CPPFLAGS =  $(DSJM_FLAGS) 
#AM_LDFLAGS = $(DSJM_FLAGS)
#AM_CXXFLAGS= -O3
EXTRA_DIST = MatrixTest.h ApplicationTest.h HeapPQTest.hh

CORE_SOURCE_FILES = $(top_srcdir)/src/mmio.cpp $(top_srcdir)/src/Matrix.cpp  $(top_srcdir)/src/CLI.cpp $(top_srcdir)/src/Utility.cc $(top_srcdir)/src/Timer.cpp  $(top_srcdir)/src/IMatrix.cc  $(top_srcdir)/src/Configuration.cc  $(top_srcdir)/src/IRowColumnDS.cc 

# Include Directory for compilation
INCLUDES 	= -I$(top_srcdir)/src
#INCLUDE_DIR = $(top_builddir)/src
#INCLUDE_HEADERS = mmio.h Matrix.h CLI.h

#MAIN_SOURCE_FILES = $(top_srcdir)/src/mmio.cpp $(top_srcdir)/src/Matrix.cpp $(top_srcdir)/src/CLI.cpp $(top_srcdir)/src/Utility.cc $(top_srcdir)/src/IMatrix.cc $(top_srcdir)/src/NNZTag.cc

TEST_SOURCE_FILES = MatrixTest.cpp unittest.cpp ApplicationTest.cpp HeapPQTest.cc BucketPQTest.cc

## Programs to Build
bin_PROGRAMS = $(TUTORIAL) $(C_SRTDAT)  $(APPLICATION)  $(VERIFY) $(TRIMONE)
## Test Programs to Build
check_PROGRAMS = $(TESTS) $(APPLICATION_TESTS)

## Mappings of the executables. / Dictionary-Like
TESTS = unittest
APPLICATION_TESTS = ApplicationApp
TUTORIAL = MatrixMarket

## APPLICATION_TESTS / ApplicationApp
ApplicationApp_SOURCES = ApplicationApp.cpp  ApplicationTest.cpp $(CORE_SOURCE_FILES)
ApplicationApp_CXXFLAGS = $(CPPUNIT_CFLAGS) $(LOG4CPP_CFLAGS)
ApplicationApp_LDFLAGS = $(CPPUNIT_LIBS) $(LOG4CPP_LIBS)

## TESTS/ unittest
unittest_SOURCES = $(CORE_SOURCE_FILES) $(TEST_SOURCE_FILES)
unittest_CXXFLAGS = $(CPPUNIT_CFLAGS) $(LOG4CPP_CFLAGS)
unittest_LDFLAGS = $(CPPUNIT_LIBS) $(LOG4CPP_LIBS)

## TUTORIAL / MatrixMarket
MatrixMarket_SOURCES =$(top_builddir)/src/mmio.cpp  example_read.cpp
MatrixMarket_CXXFLAGS = $(CPP_PROFILING_FLAGS)


## Verification Program
VERIFY = Verify
Verify_SOURCES = verify.cpp
Verify_CXXFLAGS=$(CPP_PROFILING_FLAGS)

## Output all but First Line : Program
TRIMONE = TrimOne
TrimOne_SOURCES = trimone.cpp
#TrimOne_CXXFLAGS=$(CPP_PROFILING_FLAGS)

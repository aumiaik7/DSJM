#ifndef MATRIXTEST_H
#define MATRIXTEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "Matrix.hh"


/** Log4CPP Disabled
#include "log4cpp/Portability.hh"
#ifdef LOG4CPP_HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <iostream>
#include "log4cpp/Category.hh"
#include "log4cpp/Appender.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/OstreamAppender.hh"
#ifdef LOG4CPP_HAVE_SYSLOG
#include "log4cpp/SyslogAppender.hh"
#endif
#include "log4cpp/Layout.hh"
#include "log4cpp/BasicLayout.hh"
#include "log4cpp/Priority.hh"
#include "log4cpp/NDC.hh"

*/


class MatrixTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MatrixTest );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testMin );
  CPPUNIT_TEST( testMax );
  CPPUNIT_TEST( testSLO );
  CPPUNIT_TEST( testIDO );
  CPPUNIT_TEST( testSrtdat );
  CPPUNIT_TEST( testSrtdat_3 );
  CPPUNIT_TEST( testSetr );
  CPPUNIT_TEST( testSetr_1 );
  // CPPUNIT_TEST( testSeq );
  // CPPUNIT_TEST( testDsm );
  // CPPUNIT_TEST( testNumSort );

  CPPUNIT_TEST(test_getColMax);
  CPPUNIT_TEST(testGetColMin);
  CPPUNIT_TEST(testGetColAvg);

  CPPUNIT_TEST(testGetRowMax);
  CPPUNIT_TEST(testGetRowMin);
  CPPUNIT_TEST(testGetRowAvg);

  CPPUNIT_TEST(testDegr);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp();
  void tearDown();

  void testConstructor();
  void testMin();
  void testMax();
  void testSLO();
  void testIDO();
  void testNumSort();
  void testSrtdat();
  void testSrtdat_3();
  void testSetr();
  void testSetr_1();
  void testSeq();
  void testDsm();
  void testDegr();


  void test_getColMax();
  void testGetColMin();
  void testGetColAvg();
  void testGetRowMax();
  void testGetRowMin();
  void testGetRowAvg();

  //private:
  //log4cpp::Appender* appender;
  //log4cpp::Category& log;

 private:
  Matrix* matrix;

};

#endif // MATRIXTEST_H

#ifndef APPLICATIONTEST_H
#define APPLICATIONTEST_H

#include <cppunit/extensions/HelperMacros.h> 

class ApplicationTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ApplicationTest ) ; 
  CPPUNIT_TEST ( testParseMethod ) ; 
  CPPUNIT_TEST_SUITE_END();

 public:
  void testParseMethod();
}; 

#endif // APPLICATIONTEST_H

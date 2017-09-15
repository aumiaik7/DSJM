#include "ApplicationTest.h"
#include "CLI.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( ApplicationTest);

void ApplicationTest::testParseMethod()
{
  cout << "ApplicationTest::testParseMethod()" << endl;
  CLI::parseMethod(NULL);
  CPPUNIT_ASSERT_EQUAL(CLI::parseMethod("lfo") ,CLI::LFO);
  CPPUNIT_ASSERT_EQUAL(CLI::parseMethod("slo") ,CLI::SLO);
  CPPUNIT_ASSERT_EQUAL(CLI::parseMethod("ido") ,CLI::IDO);
  CPPUNIT_ASSERT_EQUAL(CLI::parseMethod("rlf") ,CLI::RLF);
  CPPUNIT_ASSERT_EQUAL(CLI::parseMethod("sdf") ,CLI::UNKNOWN);
  CPPUNIT_ASSERT_EQUAL(CLI::parseMethod(NULL) ,CLI::UNKNOWN);
}

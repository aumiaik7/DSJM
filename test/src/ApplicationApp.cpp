#include <cstdio> 
#include <cppunit/CompilerOutputter.h> 
#include <cppunit/extensions/TestFactoryRegistry.h> 
#include <cppunit/ui/text/TestRunner.h>
#include <iostream> 
#include "ApplicationTest.h"
#include <cppunit/TestResult.h>

using namespace std; 

int main(int argc , char **argv)
{
  /* CppUnit::TestSuite suite; 
  CppUnit::TestResult result; 

  suite.addTest(new CppUnit::TestCaller<ApplicationTest>(
							 "testParseMethod",
							 &ApplicationTest::testParseMethod)); 
  */
  // Adds the test to the list of test to run. 
  CppUnit::TextUi::TestRunner runner; 
  // runner.addTest(&suite);
  runner.addTest( ApplicationTest::suite() ); 

  // Change the default outputter to a compiler error format outputter 
  runner.setOutputter( new CppUnit::CompilerOutputter ( &runner.result(), 
							std::cerr) ) ; 

  // Run the tests. 
  bool wasSuccessfull = runner.run(); 

  // Return error code 1 if one of the test failed. 
  return wasSuccessfull ? 0 : 1; 
} 

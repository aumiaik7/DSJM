#include <stdio.h>
#include <cppunit/CompilerOutputter.h> 
#include <cppunit/extensions/TestFactoryRegistry.h> 
#include <cppunit/ui/text/TestRunner.h> 
#include <cppunit/TestResult.h>


int main( int argc, char* argv[])
{

  // Get the top level suite from the registry 
  CppUnit::Test *suite = CppUnit::TestFactoryRegistry::getRegistry().makeTest();

  // Adds the test to the list of test to run. 
  CppUnit::TextUi::TestRunner runner; 
  runner.addTest(suite);
  
  /**
   * If I want to run a specific tests then I go with the following 
   */ 

  // ----------------------------------------
  /* CppUnit::TestSuite adhocSuite; 
  CppUnit::TestResult result; 

  adhocSuite.addTest(new CppUnit::TestCaller<ApplicationTest>(
							 "testParseMethod",
							 &ApplicationTest::testParseMethod)); 
							 runner.addTest(&adhocSuite); 
  */
  // ----------------------------------------

  // Change the default outputter to a compiler error format outputter
  runner.setOutputter( new CppUnit::CompilerOutputter( &runner.result(), 
						       std::cerr) ); 

  // Run the tests. 
  bool wasSuccessful = runner.run(); 

  // Return error code 1 if the one of test failed. 
  return wasSuccessful ? 0 : 1; 
}

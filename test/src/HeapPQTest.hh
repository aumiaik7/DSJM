#ifndef HEAPPQTEST_HH
#define HEAPPQTEST_HH


#include <cppunit/extensions/HelperMacros.h>

class HeapPQTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(HeapPQTest);
    CPPUNIT_TEST(testPriorityQueue);
    CPPUNIT_TEST_SUITE_END();
public:
    void testPriorityQueue();
};


#endif

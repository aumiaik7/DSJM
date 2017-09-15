#ifndef BUCKETPQTEST_HH
#define BUCKETPQTEST_HH


#include <cppunit/extensions/HelperMacros.h>

class BucketPQTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(BucketPQTest);

    CPPUNIT_TEST(test_pop);
    CPPUNIT_TEST(test_top);

    CPPUNIT_TEST_SUITE_END();
public:

    void setUp();
    void tearDown();

    void test_top();
    void test_pop();

};


#endif

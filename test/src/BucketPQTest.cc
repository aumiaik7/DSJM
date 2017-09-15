#include "BucketPQTest.hh"

#include "BucketPQ.hh"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( BucketPQTest );

/**
 * @Description : This functions
 * @Creation Date: 17.09.2009
 */
void
BucketPQTest::setUp()
{

}
/* setUp() ENDS*/

/**
 * @Description : This functions
 * @Creation Date: 17.09.2009
 */
void
BucketPQTest::tearDown()
{

}
/* tearDown() ENDS*/



/**
 * @Description :
 * @Creation Date: 17.09.2009
 */
void
BucketPQTest::test_top()
{

    std::cout << "BucketPQTest::test_top() " << std::endl;
    BucketPQ<MaxQueue> p_queue(5,6);

    p_queue.insert(1,5);
    p_queue.insert(2,5);
    p_queue.insert(3,5);
    p_queue.insert(4,5);
    p_queue.insert(5,5);
    p_queue.insert(6,5);

    while(!p_queue.empty())
    {
        std::cout <<"Top " << p_queue.top().index << std::endl;
        p_queue.pop();
    }
    CPPUNIT_FAIL("Test Case not Implemented");

    std::cout << "BucketPQTest::test_top() ends " << std::endl;
}
/* test_top() ENDS*/

/**
 * @Description :
 * @Creation Date: 17.09.2009
 */
void
BucketPQTest::test_pop()
{
    CPPUNIT_FAIL("Test Case not Implemented");
}
/* test_pop() ENDS*/



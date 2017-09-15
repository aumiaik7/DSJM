#include "HeapPQTest.hh"

#include "HeapPQ.hh"
#include "BucketPQ.hh"



// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( HeapPQTest );


template <template <typename > class T_PriorityQueue>
void run_PriorityQueueBasicTest()
{
    T_PriorityQueue<MinQueue> priority_queue(10,10);


    std::vector<Item> vertices;
    std::vector<Item>::iterator v_iter;

    vertices.push_back(Item(2,2));
    vertices.push_back(Item(1,3));
    vertices.push_back(Item(3,1));


    priority_queue.insert(2,2);
    priority_queue.insert(1,3);
    priority_queue.insert(3,1);

    while( !priority_queue.empty() )
    {
        Item item = priority_queue.top();
        // std::cout << "(" << item.first << "," << item.second << ")" << std::endl;
        priority_queue.pop();
    }

    T_PriorityQueue<MinQueue> p_queue_min(10,10);
    getPriorityQueue(vertices,p_queue_min);


    Item item;

    item = p_queue_min.top();
    p_queue_min.pop();
    CPPUNIT_ASSERT_EQUAL(item.index,3);

    item = p_queue_min.top();
    p_queue_min.pop();
    CPPUNIT_ASSERT_EQUAL(item.index,2);

    item = p_queue_min.top();
    p_queue_min.pop();
    CPPUNIT_ASSERT_EQUAL(item.index,1);



    T_PriorityQueue<MaxQueue> p_queue_max(10,10);
    getPriorityQueue(vertices,p_queue_max);

    item = p_queue_max.top();
    p_queue_max.pop();
    CPPUNIT_ASSERT_EQUAL(item.index,1);

    item = p_queue_max.top();
    p_queue_max.pop();
    CPPUNIT_ASSERT_EQUAL(item.index,2);

    item = p_queue_max.top();
    p_queue_max.pop();
    CPPUNIT_ASSERT_EQUAL(item.index,3);


    CPPUNIT_ASSERT_EQUAL(p_queue_max.empty(), true);


    T_PriorityQueue<MaxQueue> p_queue(10,10);
    getPriorityQueue(vertices,p_queue);

    p_queue.remove(1);

    item = p_queue.top();
    p_queue.pop();
    CPPUNIT_ASSERT_EQUAL(item.index,2);

    item = p_queue_max.top();
    p_queue_max.pop();
    CPPUNIT_ASSERT_EQUAL(item.index,3);

    CPPUNIT_FAIL("Test Case Not Implemented");
}

template<typename ComparisonPolicy, template <typename> class T_PriorityQueue>
void getPriorityQueue(std::vector<Item>& vertices,T_PriorityQueue<ComparisonPolicy>& p_queue)
{
    p_queue.make_empty();
    std::vector<Item>::iterator  v;
    for (v = vertices.begin(); v != vertices.end() ; v++)
    {
        p_queue.insert((*v).index,(*v).priority);
    }
}

/**
 * @Description :
 * @Creation Date: 16.09.2009
 */
void
HeapPQTest::testPriorityQueue()
{
    run_PriorityQueueBasicTest< HeapPQ >();
    run_PriorityQueueBasicTest< BucketPQ >();
}
/* testPriorityQueue() ENDS*/

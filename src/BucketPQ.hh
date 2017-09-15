#ifndef BUCKETPQ_HH
#define BUCKETPQ_HH

#include "detail/PriorityQueue.hh"
#include <vector>
#include <cassert>
#include <iostream>

template <typename MinMaxQueue>
class BucketPQ
{
private:

    const int _invalid;

    std::vector<int> head;
    std::vector<int> next;
    std::vector<int> prev;

    std::vector<Priority> v_priority;
    int maxBucket;
    int numberOfItems;
    Priority optimum;
    int size;

    MinMaxQueue cmp_fun;


public:
    BucketPQ(int _maxBucket, int _numberOfItems)
        :maxBucket(_maxBucket), numberOfItems(_numberOfItems), size(0),_invalid(-1)
    {
        head.reserve(_maxBucket+1);
        next.reserve(_numberOfItems + 1);
        prev.reserve(_numberOfItems + 1);
        v_priority.reserve(_numberOfItems + 1);
        optimum = cmp_fun.opt_init_value(maxBucket);

        make_empty(); // Setting head, next, and prev to their default
                      // value.
    }

    bool empty() const
    {
        return size == 0;
    }

    void insert(Index vertex_id, Priority priority)
    {
        size++;
        optimum = cmp_fun.opt(priority,optimum);
        v_priority[vertex_id] = priority;

        prev[vertex_id] = _invalid;
        next[vertex_id] = head[priority];
        if( head[priority] > _invalid)
        {
            prev[head[priority]] = vertex_id;
        }
        head[priority] = vertex_id;
    }


    Item top()
    {
        int jcol;
        do
        {
            jcol = head[optimum];
            if(jcol > _invalid) // BUG LOCATION: RISK,
            {             // As we are using zero based indexing in
                          // our PriorityQueue Implementation. So ,be
                          // careful about it.
                break;
            }

            optimum = cmp_fun.next_op(optimum);
            assert(optimum > _invalid);
        }while(1);
        assert(v_priority[jcol] == optimum);
        return Item( jcol, v_priority[jcol] );
    }

    void pop()
    {
        remove(head[optimum]);
    }

    Item get(Index vertex_id)
    {
        return Item(vertex_id, v_priority[vertex_id]);
    }

    void increase(Index vertex_id)
    {
        optimum = cmp_fun.opt(optimum, v_priority[vertex_id] + 1);
        change(vertex_id, v_priority[vertex_id]+1);
    }

    void decrease(Index vertex_id)
    {
        optimum = cmp_fun.opt(optimum, v_priority[vertex_id] -1);
        change(vertex_id, v_priority[vertex_id]-1);
    }

    void change(Index vertex_id, Priority priority)
    {
        remove(vertex_id);
        insert(vertex_id, priority);
    }

    void remove(Index vertex_id)
    {
        size--;
        Priority priority = v_priority[vertex_id];
        if(prev[vertex_id] == _invalid)
        {
            head[priority] = next[vertex_id];
        }
        else
        {
            next[prev[vertex_id]] = next[vertex_id];
        }

        if(next[vertex_id] > _invalid)
        {
            prev[next[vertex_id]] = prev[vertex_id];
        }
    }

    void make_empty()
    {
        size = 0;
        for(int i = 0; i < numberOfItems; i++)
        {
            if(i <= maxBucket)
                head[i] = _invalid;
            next[i] = _invalid;
            prev[i] = _invalid;
        }
    }

    // Needed for IDO.
    Item betterTop(int jp, const int maxlst, const int *ndeg)
    {
        jp = jp;
        int numwgt = -1;
        int jcol;
        for(int numlst = 1; numlst <= maxlst; numlst++)
        {
            assert(ndeg[jp] > _invalid);
            if (ndeg[jp] > numwgt)
            {
                numwgt = ndeg[jp];
                jcol = jp;
            }
            jp = next[jp];
            if (jp <= 0)
                break;
        }
        return Item(jcol, v_priority[jcol]);
    }
};

#endif

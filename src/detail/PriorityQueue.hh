#ifndef PRIORITYQUEUE_HH
#define PRIORITYQUEUE_HH

#define startIndex 1

#include <algorithm>

typedef int Priority;
typedef int Index;

template<typename T>
inline T PQ_STEP_UP(T t)
{
    return t + startIndex;
}

template<typename T>
inline T PQ_STEP_DOWN(T t)
{
    return t - startIndex;
}

typedef class PriorityQueueItem
{
public:
    Index index;
    Priority priority;

    PriorityQueueItem(Index i, Priority p): index(i), priority(p){}
    PriorityQueueItem() {}

}Item;

struct MinQueue
{
    bool operator()(Item &a, Item &b)
    {
        return a.priority > b.priority;
    }

    int opt_init_value(int maxBucket)
    {
        return maxBucket;
        // return (std::numeric_limits<size_t>::max)();
    }

    int opt(Priority priority, Priority _min)
    {
        if ( priority < _min)
            return priority;
        else
            return _min;
    }

    int next_op(int optimum)
    {
        return optimum +1 ;
    }
};

struct MaxQueue
{
    bool operator() (Item &a, Item &b)
    {
        return a.priority < b.priority;
    }

    int opt_init_value(int maxBucket)
    {
        return -1;
    }

    int opt(Priority priority, Priority _max)
    {
        if ( priority > _max)
            return priority;
        else return _max;
    }

    int next_op(int optimum)
    {
        return optimum - 1;
    }
};



// typedef std::pair<Index, Priority> Item;


#endif

#ifndef RLF_HH
#define RLF_HH

#include "Utility.h"
#include "OrderingMethod.hh"

using namespace std;

class RLF
{
public:
    const OrderingMethod::ordering_direction ord_direction;
private:
    int M;
    int N;
    int *jpntr;
    int *indRow;
    int *ipntr;
    int *indCol;
    int *ndeg;
    int *list;
    int *ngrp;
    int numord;
    BucketPQ<MaxQueue> priority_queue;
    BucketPQ<MaxQueue> u_queue;
    std::vector<int> tag;

    int *blackList;
    int *u_list;
    bool *inU;
    int *u_tag;
    bool newColorClass;
    int u_maxdeg ;
    int u_numdeg;
    int maxdeg;

    int howMany;

    int countV;
    int countV_prime;
    int countC;
    int q;

    int columnRemaining;
    int countColored;

public:
    // Constructor .
    template <typename PreviousOrderingMethod>
    RLF(int M_, int N_, int *ndeg_, int *jpntr_,
        int *indRow_, int *ipntr_, int *indCol_, int *list_, int *ngrp_,
        int maxBucket_, PreviousOrderingMethod& previousOrderingMethod, int howMany_)
        :ord_direction(OrderingMethod::FORWARD),
         M(M_),
         N(N_),
         jpntr(jpntr_),
         indRow(indRow_),
         ipntr(ipntr_),
         indCol(indCol_),
         ndeg(ndeg_),
         list(list_),
         ngrp(ngrp_),
         priority_queue(maxBucket_, N_),
         u_queue(maxBucket_, N_),
         howMany(howMany_),
         countColored(0)
    {

        // Initialize Numord.
        if(previousOrderingMethod.ord_direction == this->ord_direction)
            numord = previousOrderingMethod.getNumOrd();
        else
            numord = 1;

        tag.reserve(N+1);
        columnRemaining = 0;
        for (int jcol = 1; jcol <= N; jcol++)
        {
            if(previousOrderingMethod.isTagged(jcol) )
            {
                // We tag it too.
                tagit(jcol);
            }
            else
            {
                ++columnRemaining;
                untagit(jcol);
            }
        }
        blackList = new int[M+1];
        for ( int i = 1 ; i <= M; i++)
        {
            blackList[i] = 0;
        }

        u_list = new int[N+1];
        inU = new bool[N+1];
        u_tag = new int[N+1];

        newColorClass = true;

        u_maxdeg = 0;
        u_numdeg = 0;

        maxdeg = 0;
        for(int jp = 1; jp <=N ; jp++)
        {
            inU[jp] = false;
            u_tag[jp] = 0;
        }

        for (int jp = 1; jp <= N; jp++)
        {
            if(!previousOrderingMethod.isTagged(jp))
            {
                priority_queue.insert(jp, previousOrderingMethod.getDeg(jp));
                u_queue.insert(jp,0);
            }
        }

    }
    // Destructor.
    ~RLF()
    {

    }
    // The actual Work.
    int work()
    {

        newColorClass = true;

        countV = columnRemaining; // We used to have N here, but it is
                                  // dependent on the remaining column
                                  // actually, when we are thinking
                                  // about the size of N.
        countV_prime = 0;
        countC = 0;
        q = 1;
        int count = 0;
        while(true)
        {
            int jcol;
            if (newColorClass == true)
            {
                newColorClass = false;
                /**
                 * Choose a column jcol of maximal degree mindeg.
                 */
                //printf("Choosing JCol from the List of MaxDegree such\N that it has the Maximum Degree in V\N" );
                Item item  = priority_queue.top();
                jcol = item.index;
                maxdeg = item.priority;

            }
            else
            {
                /**
                 * Choose a column jcol that has maximal degree in
                 * set U.
                 */
                //printf("Choosing JCol from the list of V such that \nit has the maximum edges to U\N");
                Item item = u_queue.top();
                jcol = item.index;
                u_maxdeg = item.priority;
            }
            countV = countV - 1;
            countC = countC + 1;

            list[numord] = jcol;
            // printf("list[%d] = %d\n",numord,jcol);

            priority_queue.remove(jcol);
            u_queue.remove(jcol);

            /**
             * Tag the column, as it has been colored.
             */
            tagit(jcol);
            ++countColored;

            numord = numord + 1;
            if(numord > N || countColored >= howMany)
            {
                // De allocate Memory.

                delete[] u_list;
                delete[] inU;
                delete[] u_tag;

                delete[] blackList;
                return q;
            }


            // TODO: Termination Test for MixUp things/ Based on Count
            // of colored vertices through this ordering method.


            /**
             * Find All Columns adjacent to column jcol.
             * Determine all positions (ir,jcol) which correspond
             * to non-zeroes in the matrix.
             */
            for(int jp = jpntr[jcol] ; jp < jpntr[jcol+1] ; jp++)
            {
                /*
                 * We do not process any of the columns found on this
                 * column, while we are updating the u_degree for each of the
                 * vertices.
                 */
                int ir = indRow[jp];
                blackList[ir] = q;
            }


            for ( int jp = jpntr[jcol] ; jp < jpntr[jcol+1]  ; jp++)
            {
                int ir = indRow[jp];

                // blackList[ir] = 1; // Make this row a black one.


                /**
                 * For each row ir, determine all positions (ir,ic)
                 * which correspond to non-zeroes in the matrix.
                 */
                for ( int ip = ipntr[ir]; ip < ipntr[ir+1]; ip++)
                {
                    int ic= indCol[ip];

                    /** Array tag marks columns whch are adjacent
                     * to column jcol
                     */

                    if(tag[ic] < numord)
                    {
                        tag[ic] = numord;
                        priority_queue.decrease(ic);

                        if (inU[ic] == false)
                        {
                            //printf("ADD ic = %d to  V_Prime\N",ic);
                            inU[ic] = true;
                            countV_prime = countV_prime + 1;
                            countV = countV - 1;

                            u_queue.remove(ic);

                            // Update the degrees of the adjacent vertices.

                            u_maxdeg = pq_updateDegreesToUVertices(N,ic,u_maxdeg, jpntr,indRow,ipntr,indCol,inU,
                                                                   tag,u_tag,u_queue,blackList,q);

                        }
                    }
                }
            }

            if ( ( countV_prime + countC) == columnRemaining)
            {
                if ( countV_prime == 0 )
                {
                    // We found the solution with total coloring q

                    // We actually never come to this portion ??? !!! ALERT
                    // This comment is misleading, not verified.

                    // De allocate Memory.
                    delete[] u_list;
                    delete[] inU;
                    delete[] u_tag;
                    delete[] blackList;

                    return q;
                    /**
                     * We found the solution
                     * */
                }
                else
                {
                    /* ----------------------------------------
                       Start a new Color Class.
                       ---------------------------------------- */
                    q = q + 1;

                    // We are reinitializing

                    /**
                     * We are supposed to initialize V, U
                     */
                    newColorClass = true;
                    countV =  countV_prime;
                    countV_prime = 0;

                    u_numdeg = 0;

                    u_maxdeg = 0;

                    pq_initializeDegreesToUVertices(N,tag,u_queue,inU,u_tag);
                }
            }
        }
        // It does not reach This point.
    }

    inline int getDeg(int jcol)
    {
        return priority_queue.get(jcol).priority;
    }

    inline int getNumOrd()
    {
        return numord;
    }

    inline void tagit(int jcol)
    {
        tagit(jcol,N);
    }

    inline void tagit(int jcol, int value)
    {
        tag[jcol] = value;
    }

    inline void untagit(int jcol)
    {
        tagit(jcol,0);
    }

    inline bool isTagged(int jcol)
    {
        return tag[jcol] == N;
    }


    template <typename T_PriorityQueue, typename Tag>
    static int pq_updateDegreesToUVertices(int n, int jcol,int maxdeg, int *jpntr,int *indRow,
                                           int *ipntr,int *indCol, bool * inU, Tag& tag, int *u_tag,
                                           T_PriorityQueue& u_queue, int *blackList,const int q)
    {
        // Update the degrees of the adjacent vertices.

        int numdeg;
        bool zero = false;

        for( int jp = jpntr[jcol] ; jp < jpntr[jcol+1] ; jp++)
        {
            int ir = indRow[jp];
            if(blackList[ir] == q)
                continue;
            for(int ip = ipntr[ir]; ip < ipntr[ir+1] ; ip++)
            {
                int ic = indCol[ip];

                if(inU[ic] == false && tag[ic] < n && u_tag[ic] != jcol)
                {
                    u_tag[ic] = jcol;
                    u_queue.increase(ic);
                    maxdeg = max(u_queue.get(ic).priority, maxdeg);

                }
            }
        }
        return maxdeg;
    }

    template <typename T_PriorityQueue, typename Tag>
    static void pq_initializeDegreesToUVertices(int n,Tag& tag,T_PriorityQueue& u_queue,
                                                bool *inU, int *u_tag)
    {

        for (int jp = 1;jp <= n; jp++)
        {
            inU[jp] = false;
            u_tag[jp] = 0;
        }

        u_queue.make_empty();
        for (int jp = 1; jp <= n ; jp++)
        {
            if(tag[jp] < n)
            {
                u_queue.insert(jp,0);
            }
        }
    }

    inline int getNumOrd() const
    {
        return numord;
    }

};

#endif

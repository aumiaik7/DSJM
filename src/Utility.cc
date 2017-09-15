//
// Utility.cc
//
// Made by Mahmudul Hasan
// Login   <hasanm@hasanm-laptop>
//
// Started on  Tue Feb 17 23:51:11 2009 Mahmudul Hasan
// Started on  Tue Feb 17 23:51:11 2009 Mahmudul Hasan
//

#include <iostream>
#include "Utility.h"
#include "mmio.h"

using namespace std;


int MatrixUtility::loadMatrix(string fileName,Matrix* matrix)
{
    // Copy The Setup Method From MatrixTest.h
    //

    MM_typecode matcode;
    int M,N,nz,NZ;
    int ret_code;

    FILE *f;

    if ( (f = fopen(fileName.c_str(),"r") ) == NULL)
    {
        // TODO: Give Feedback / Exception Error Condition
        return -1;
    }

    if ( mm_read_banner (f , &matcode ) != 0)
    {
        // TODO: Give Feedback / Exception Error Condition
        return -1;
    }

    if( mm_is_complex( matcode) && mm_is_matrix( matcode) &&
        mm_is_sparse(matcode))
    {
        // TODO: Check these Functions Definitions, I don't have enough
        // ideas on them.

        // TODO: Give Feedback/ Exception
        return -1;
    }

    if( (ret_code = mm_read_mtx_crd_size ( f, &M, &N, &nz ) ) != 0)
    {
        // TODO: Error
        return -1;
    }

    NZ = nz;

    if( mm_is_symmetric(matcode))
    {
        NZ = 2 *  nz;
    }

    matrix = new Matrix(M,N,NZ, false);

    double value;
    int row, col;

    for (int i ; i <= nz; i++ )
    {
        if ( mm_is_pattern( matcode) )
        {
            fscanf(f, "%d %d\n", &row, &col);
        }
        else
        {
            fscanf(f, "%d %d %lg\n", &row, &col, &value);
        }

        matrix->setIndRowEntry(i,row);
        matrix->setIndColEntry(i,col);

        if(mm_is_symmetric(matcode))
        {
            matrix->setIndRowEntry(i+nz,col);
            matrix->setIndColEntry(i+nz,row);
        }
    }

    if ( f!= stdin )
        fclose(f);

    matrix->computeCCS();
    int nnz = matrix->compress();
    matrix->computeCRS();

    return nnz;

}

void MatrixUtility::indexsort(int n,int nmax,int *num,int mode,int *index,int *last,int *next)
{

    //  subroutine indexsort(n,nmax,num,mode,index,last,next)
    //  integer n,nmax,mode
    //  integer num(n),index(n),last(0:nmax),next(n)
    /*     **********.
     *
     *     subroutine numsrt
     *
     *     given a sequence of integers, this subroutine groups
     *     together those indices with the same sequence value
     *     and, optionally, sorts the sequence into either
     *     ascending or descending order.
     *
     *     the sequence of integers is defined by the array num,
     *     and it is assumed that the integers are each from the set
     *     0,1,...,nmax. on output the indices k such that num(k) = l
     *     for any l = 0,1,...,nmax can be obtained from the arrays
     *     last and next as follows.
     *
     *           k = last(l)
     *           while (k != 0) k = next(k)
     *
     *     optionally, the subroutine produces an array index so that
     *     the sequence num(index(i)), i = 1,2,...,n is sorted.
     *
     *     the subroutine statement is
     *
     *       subroutine numsrt(n,nmax,num,mode,index,last,next)
     *
     *     where
     *
     *       n is a positive integer input variable.
     *
     *       nmax is a positive integer input variable.
     *
     *       num is an input array of length n which contains the
     *         sequence of integers to be grouped and sorted. it
     *         is assumed that the integers are each from the set
     *         0,1,...,nmax.
     *
     *       mode is an integer input variable. the sequence num is
     *         sorted in ascending order if mode is positive and in
     *         descending order if mode is negative. if mode is 0,
     *         no sorting is done.
     *
     *       index is an integer output array of length n set so
     *         that the sequence
     *
     *               num(index(i)), i = 1,2,...,n
     *
     *         is sorted according to the setting of mode. if mode
     *         is 0, index is not referenced.
     *
     *       last is an integer output array of length nmax + 1. the
     *         index of num for the last occurrence of l is last(l)
     *         for any l = 0,1,...,nmax unless last(l) = 0. in
     *         this case l does not appear in num.
     *
     *       next is an integer output array of length n. if
     *         num(k) = l, then the index of num for the previous
     *         occurrence of l is next(k) for any l = 0,1,...,nmax
     *         unless next(k) = 0. in this case there is no previous
     *         occurrence of l in num.
     *
     *     argonne national laboratory. minpack project. july 1983.
     *     thomas f. coleman, burton s. garbow, jorge j. more'
     *
     */
    int i,j,jinc,jl,ju,k,l;
    /*
     *     determine the arrays next and last.
     */

    for(int i = 0; i <=  nmax; i++)
    {
        last[i] = 0;
    }

    for(int  k = 1; k <=  n ; k++)
    {
        l = num[k];
        next[k] = last[l];
        last[l] = k;
    }
    if (mode == 0)
        return;
    /*
     *     store the pointers to the sorted array in index.
     */
    i = 1;
    if (mode > 0)
    {
        jl = 0;
        ju = nmax;
        jinc = 1;
    }
    else
    {
        jl = nmax;
        ju = 0;
        jinc = -1;
    }

    for(int j = jl; j*jinc <= ju*jinc; j = j + jinc)
    {
        k = last[j];
        while(k!=0)
        {
        	index[i] = k;
        	i = i + 1;
        	k = next[k];
        }
    }

    return;
}

#ifndef UTILITY_H
#define UTILITY_H

/*
** Utility.h
**
** Made by Mahmudul Hasan
** Login   <hasanm@hasanm-laptop>
**
** Started on  Tue Feb 17 23:51:03 2009 Mahmudul Hasan
** Last update Tue Feb 17 23:51:03 2009 Mahmudul Hasan
*/

#ifndef   	UTILITY_H_
# define   	UTILITY_H_

#include "Matrix.hh"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>


namespace MatrixUtility
{
    /**

     * Int Return Value for error code
     */
    int loadMatrix(std::string fileName, Matrix* matrix);
	inline std::string createString(char* arg)
    {
        std::string s(arg);
        return s;
    }

    inline int square(int a)
    {
        return a * a;
    }
    /**
     * Method Numsrt()
     * Originally developed in Coleman's Implementation.
     */
    void indexsort(int n,int nmax,int *num,int mode,int *index,int *last,int *next);



}

namespace DSJM
{
    class MemoryError
    {
    };

}
#endif 	    /* !UTILITY_H_ */


#endif

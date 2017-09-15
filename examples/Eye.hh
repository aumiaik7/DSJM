#ifndef EYE_HH
#define EYE_HH

#include "Matrix.hh"

class Eye: public Matrix
{
public:
    Eye(int n):Matrix(n,n,n, false)
    {
        for(int i = 1; i<=n; i++)
        {
            setIndRowEntry(i,i);
            setIndColEntry(i,i);
        }

        computeCCS();
        compress();
        computeCRS();
    }
    virtual ~Eye()
    {

    }
};

#endif

#ifndef FD_HH
#define FD_HH

#include <cassert>
#include <memory>
#include <iostream>

using namespace std;

class FD
{
public:
    FD(){}
    virtual ~FD(){}

    template<class Function,class SeedMatrix>
    static Matrix* finiteDifference(Function f,const vector<double> &x,SeedMatrix* seed)
    {
        const double h = 0.0000005;
        Matrix *mT;
        int m, n, nz;
        m = f.m; n = f.n; nz = f.getNzmax(m);
        mT = new Matrix(m,seed->getN(),nz,true);

        int count = 0;
        auto_ptr<vector<double> > y =  f.eval(x);


        for(int j = 1; j <= seed->getN(); j++)
        {

            vector<double> x1 = x;
            for(int ic = seed->getJpntrEntry(j); ic < seed->getJpntrEntry(j+1);ic++)
            {
                int ir = seed->getIndRowEntry(ic);
                x1[ir] = x1[ir] + h;
            }

            auto_ptr<vector<double> > y1 = f.eval(x1);


            for(int i = 1; i <= f.m; i++)
            {
                double val = ((*y1)[i] -  (*y)[i])/h;

                if(fabs(val) > 0.000000000001)
                {
                    mT->entry(i,j,val);
                    count++;
                }

            }
        }
        mT->computeCCS();
        mT->compress();
        mT->computeCRS();
        return mT;
    }
};

#endif

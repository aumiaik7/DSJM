#ifndef MATRIX_H
#define MATRIX_H

#include "IMatrix.hh"
#include "NNZTag.hh"
#include "BucketPQ.hh"

#include <iostream>
#include <cmath>
#include "CLI.h"
#include <boost/dynamic_bitset.hpp>



///
///
/// \todo
//	This class need better namings for the variables.
class Matrix : public IMatrix
{

public:
    ~Matrix();
    Matrix(int M, int N, int nz, bool value);

    int updateDegreesToUVertices(int n, int ic, int u_maxdeg, int *jpntr,int *indRow,int *ipntr,int *indCol,
                                 bool * f_added, int *tag, int *f_tag, int *u_list,
                                 int *u_head, int *u_next, int *u_previous,int *list,int *blackList,const int q);

    // ========================================
    // IMatrix Methods
    bool slo(int *order);

    /**
     * Method ido()
     **/
    bool ido(int *order);
    bool idoDsatur(int *order, int *clique);	
    bool lfo(int *order);


    
    /*int getColumn();
    //int branchColor(int order,int currentColor, int *color);	
    int branchColor(int order,int currentColor);    
    //bool colorAvailable(int jcol, int colorNo, int *color);	
    bool colorAvailable(int jcol, int colorNo); 
    void satDegInc(int jcol,int order, int colorNo);
    void satDegDec(int jcol,int order, int colorNo);*/		



    bool computedegree();
    int greedycolor(int *list, int *ngrp);
    int rlf(int *ngrp);
    int sdo(int *ngrp);
    int sdo2(int *ngrp);
    int slo_rlf(int *list, int *ngrp);
    //int dsatur(int *ngrp, int *clique, int UB);
    //int dsatur(int *clique, int UB, int tbChoice);    
    int dsatur(int UB, int tbChoice);    

    Matrix* getSeedMatrix(int *ngrp);

    int getNumberOfColors() const;

    // IMatrix Methods
    // ========================================

private:
    /**
     * Private Constructor
     */
    Matrix();

    int *ndeg;
    int maxi, maxj;
    int rho_max;	

    int numberOfColors; // TODO: Remove

    double *val; // This is supposed to hold the values,
    // Currently we are not using it anywhere.
};

#endif

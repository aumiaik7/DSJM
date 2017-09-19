#ifndef IMATRIX_HH
#define IMATRIX_HH

#include "IRowColumnDS.hh"
#include "BucketPQ.hh"

#include <string>


class IMatrix : public IRowColumnDS
{
public:
    IMatrix(int M,int N, int nz, bool value);
    virtual ~IMatrix();


    /**
     * Method slo()
     */
    virtual bool slo(int *order) = 0;

    /**
     * Method ido()
     **/
    virtual bool ido(int *order) = 0;
    virtual bool idoDsatur(int *order,int *clique) = 0;	
    virtual bool lfo(int *order) = 0 ;

    virtual bool computedegree() = 0 ;
    virtual int greedycolor(int *list, int *ngrp) = 0 ;
    virtual int rlf(int *ngrp) = 0;
    virtual int dsatur(int UB, int tbChoice) = 0;
    virtual int slo_rlf(int *list, int *ngrp) = 0 ;	  

    void setVerify(bool v);
    bool getVerify();



    void buildPriorityQueue(int n, int *ndeg, int *head, int *next, int *previous);

    void deleteColumn(int *head,int *next,int *previous,int numdeg,int jcol);
    //void deleteColumn2(int **head,int *next,int *previous,int* headj,int numdeg,int jcol); 	
    void deleteColumn2(int **head,int *next,int *previous,int ideg,int numdeg,int jcol);  
    void addColumn(int *head,int *next,int *previous,int numdeg,int jcol);
    //void addColumn2(int **head,int *next,int *previous,int* headj, int* inducedDeg,int numdeg,int jcol);
    void addColumn2(int **head,int *next,int *previous, int ideg,int numdeg,int jcol);

    void initializeDegreesToUVertices(int n,int *tag,int *u_head,int *u_next,int *u_previous,
                                      int *u_list, bool *f_added, int *f_tag);


    //dsat related methods
    int getColumn();
    //int branchColor(int order,int currentColor, int *color);  
    int branchColor(int order,int currentColor);    
    //bool colorAvailable(int jcol, int colorNo, int *color);   
    bool colorAvailable(int jcol, int colorNo); 
    void satDegInc(int jcol,int order, int colorNo);
    void satDegDec(int jcol,int order, int colorNo);


    int *satDegDsat;
int *headDsat;
int *nextDsat;
int *previousDsat;
int *tagDsat;
int *seqTagDsat;
int maxgrpDsat;
int *inducedDegDsat; 
int **colorTracker;

int LB;//lower bound
int UB;//upper bound  
int maxsatDsat;
int *colorDsat;
int maxlstDsat; 
bool *handled;
int startTime,currentTime; //to keep track of timing 
double subProblems;
int tbChoice;

std::string name;




    // ========================================
    // RowColumnDS
protected:


    bool shouldVerify;
    int maxdeg;
};

#endif

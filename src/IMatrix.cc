#include "IMatrix.hh"
#include <iostream>
#include <cmath>
#include "IRowColumnDS.hh"

#define CUT_OFF_TIME 60 //how many minutes dsatur will run if it doesn't terminate

using namespace std;

/**
 * @Description : This functions
 * @Creation Date: 20.06.2009
 */
IMatrix::IMatrix(int M, int N, int nz, bool value)
    : IRowColumnDS(M,N,nz, value),
      shouldVerify(false)
{

}
/* IMatrix() ENDS*/

IMatrix::~IMatrix()
{

}

void IMatrix::buildPriorityQueue(int n, int *ndeg, int *head, int *next, int *previous)
{
    int numdeg;
    for (int jp = 1; jp <= n ; jp++)
    {
        numdeg = ndeg[jp];

        previous[jp] = 0;

        next[jp] = head[numdeg];

        if ( head[numdeg] > 0 )
        {
            previous[head[numdeg]] = jp;
        }
        head[numdeg] = jp;
    }
}

/**
 * Constant Operation
 */
void IMatrix::deleteColumn(int *head,int *next,int *previous,int numdeg,int jcol)
{
    if ( previous[jcol] == 0)
    {
        head[numdeg] = next[jcol];
    }
    else
    {
        next[previous[jcol]] = next[jcol];
    }

    if(next[jcol] > 0)
    {
        previous[next[jcol]] = previous[jcol];
    }
     next[jcol] = previous[jcol] = 0;	
}

//unlike deletecolumn method this method deletes a column from a bucket with 2d head array
void IMatrix::deleteColumn2(int **head,int *next,int *previous,int ideg, int numdeg,int jcol)
{
    if ( previous[jcol] == 0)
    {
        head[numdeg][ideg] = next[jcol];
    }
    else
    {
        next[previous[jcol]] = next[jcol];
    }

    if(next[jcol] > 0)
    {
        previous[next[jcol]] = previous[jcol];
    }
     next[jcol] = previous[jcol] = 0;

}

/**
 * Constant Operation
 */
void IMatrix::addColumn(int *head,int *next,int *previous, int numdeg, int jcol)
{

    previous[jcol] = 0;
    next[jcol] = head[numdeg];
    if(head[numdeg] > 0)
    {
        previous[head[numdeg]] = jcol;
    }
    head[numdeg] = jcol;
}

//unlike addcolumn method this method adds a column to a bucket with 2d head array
void IMatrix::addColumn2(int **head,int *next,int *previous, int ideg, int numdeg, int jcol)
{
   	previous[jcol] = 0;
        next[jcol] = head[numdeg][ideg];
        if(head[numdeg][ideg]>0)
        {
            previous[head[numdeg][ideg]] = jcol;
        }
        head[numdeg][ideg] = jcol;
}

void IMatrix::initializeDegreesToUVertices(int n,int *tag,int *u_head,int *u_next,int *u_previous, int *u_list, bool *inU, int *u_tag)
{

    for (int jp = 1;jp <= n; jp++)
    {
        u_list[jp] = 0;
        u_head[jp-1] = 0;
        u_next[jp] = 0;
        u_previous[jp] = 0;
        inU[jp] = false;
        u_tag[jp] = 0;
    }
    for (int jp = 1; jp <= n ; jp++)
    {
        if(tag[jp] < n)
        {
            u_previous[jp] = 0;
            u_next[jp] = u_head[0] ;
            if ( u_head[0] > 0)
            {
                u_previous[u_head[0]] = jp;
            }
            u_head[0] = jp;
        }
    }
}

void IMatrix::setVerify(bool v)
{
    shouldVerify = v;
}

bool IMatrix::getVerify()
{
    return shouldVerify;
}


//methods used for DSATUR based exact coloring alagorithms
//this method is used to select a column in each iteration 
//based on user's tie-breking choice(1-4) 
int IMatrix::getColumn()
{
    int column;
    int maxsatLevel = maxsatDsat;
    //it simply retruns the head of maxsat
    if(tbChoice==1)
    {   
        
        while(true)
        {
            column = headDsat[maxsatLevel];
            //jp = headDsat[maxsatDsat];
            if(column>0)
                break;
            maxsatLevel--;
            //maxsatDsat--;
        }
        return column;
    }
    //maximum common unused colors-1
    else if(tbChoice==2)
    {
        while(true)
        {

            column = headDsat[maxsatLevel];
            //jp = headDsat[maxsatDsat];
            if(column>0)
            {
               int maxCol = column; 
               int maxCount = 0;    
               while(true)
               {
                    int count = 0;
                    bool *tagCol = new bool[N+1]();
                    tagCol[column] = true;
                    for (int jp = jpntr[column] ; jp < jpntr[column+1] ; jp++)
                    {
                        int ir = row_ind[jp];

                        for( int ip = ipntr[ir] ; ip < ipntr[ir+1] ; ip++)
                        {
                            int ic = col_ind[ip];
                            if(!tagCol[ic] && !handled[ic])
                            {   
                                //it is now tagged so will not be processed in any next iteration
                                tagCol[ic] = true;  
                                for(int color=1; color<=maxsatDsat; color++)
                                {
                                     if(colorTracker[color][column] == 0 && colorTracker[color][ic] == 0)
                                     {
                                        count++;
                                        if(count>maxCount)
                                        {
                                            maxCount = count;
                                            maxCol = column;
                                        }
                                     }      
                                }
                            }
                        }
                    }
                    if(tagCol) delete [] tagCol;
                    if(nextDsat[column]<=0)
                    {
                        return maxCol;
                    }
                    column = nextDsat[column];
               
               }

            }    
            maxsatLevel--;
            //maxsatDsat--;
        }
    }
    //maximum common unused colors-3
    else if(tbChoice==3)
    {
        while(true)
        {

            column = headDsat[maxsatLevel];
            //jp = headDsat[maxsatDsat];
            if(column>0)
            {
               int maxCol = column; 
               int maxCount = 0;    
               while(true)
               {
                    int count = 0;
                    bool *tagCol = new bool[N+1]();
                    tagCol[column] = true;
                    for (int jp = jpntr[column] ; jp < jpntr[column+1] ; jp++)
                    {
                        int ir = row_ind[jp];

                        for( int ip = ipntr[ir] ; ip < ipntr[ir+1] ; ip++)
                        {
                            int ic = col_ind[ip];
                            if(satDegDsat[ic] != satDegDsat[column])
                                continue;
                            if(!tagCol[ic] && !handled[ic] )
                            {   
                                //it is now tagged so will not be processed in any next iteration
                                tagCol[ic] = true;  
                                for(int color=1; color<=maxsatDsat; color++)
                                {
                                     if(colorTracker[color][column] == 0 && colorTracker[color][ic] == 0)
                                     {
                                        count++;
                                        if(count>maxCount)
                                        {
                                            maxCount = count;
                                            maxCol = column;
                                        }
                                     }      
                                }
                            }
                        }
                    }
                    if(tagCol) delete [] tagCol;
                    if(nextDsat[column]<=0)
                    {
                        return maxCol;
                    }
                    column = nextDsat[column];
               
               }

            }    
            maxsatLevel--;
            //maxsatDsat--;
        }
    }
    //maximum common unused colors-4
    else if(tbChoice==4)
    {
        while(true)
        {

            column = headDsat[maxsatLevel];
            //jp = headDsat[maxsatDsat];
            if(column>0)
            {
               int maxCol = column; 
               int maxCount = 0;    
               while(true)
               {
                    int count = 0;
                    bool *tagCol = new bool[N+1]();
                    tagCol[column] = true;
                    for (int jp = jpntr[column] ; jp < jpntr[column+1] ; jp++)
                    {
                        int ir = row_ind[jp];

                        for( int ip = ipntr[ir] ; ip < ipntr[ir+1] ; ip++)
                        {
                            int ic = col_ind[ip];
                            if(satDegDsat[ic]<1)
                                continue;
                            if(!tagCol[ic] && !handled[ic] )
                            {   
                                //it is now tagged so will not be processed in any next iteration
                                tagCol[ic] = true;  
                                for(int color=1; color<=maxsatDsat; color++)
                                {
                                     if(colorTracker[color][column] == 0 && colorTracker[color][ic] == 0)
                                     {
                                        count++;
                                        if(count>maxCount)
                                        {
                                            maxCount = count;
                                            maxCol = column;
                                        }
                                     }      
                                }
                            }
                        }
                    }
                    if(tagCol) delete [] tagCol;
                    if(nextDsat[column]<=0)
                    {
                        return maxCol;
                    }
                    column = nextDsat[column];
               
               }

            }    
            maxsatLevel--;
            //maxsatDsat--;
        }
    }
}   

//After coloring a colunm we look for all of its adjacent columns
//and update their saturation degree where applicable
void IMatrix::satDegInc(int jcol,int order,int colorNo)
{
    //this array is for simple tagging. we select a nonzero element of jcol
    //if we process a neighboring column of jcol we dont need to process that
    //column in any next iteration. so we tag that column and don't procees that  
    bool *tagCol = new bool[N+1]();
    tagCol[jcol] = true;
    for (int jp = jpntr[jcol] ; jp < jpntr[jcol+1] ; jp++)
    {
        int ir = row_ind[jp];

        for( int ip = ipntr[ir] ; ip < ipntr[ir+1] ; ip++)
        {
            int ic = col_ind[ip];
           
         
            //column ic has not been processed
            if(!tagCol[ic])
            {   
                //it is now tagged so will not be processed in any next iteration
                tagCol[ic] = true;  
                
                //noumber of "colorNO" colored neighbors in ordered graph
                int prevColorCount = colorTracker[colorNo][ic];
                                        
                //no "colorNo" colored neighbors in ordered graph so we can increase the saturation degree of ic                        
                if(prevColorCount==0)// && toggleFlag[ic] ==0)
                {
                    
                    satDegDsat[ic]++;
                    
                    // update the maxsat.
                    maxsatDsat = max(maxsatDsat,satDegDsat[ic]);
                   
                    //this means the column isn't handled 
                    //and in the bucket so we can update its
                    //degree list (sat degree) in bucket 
                    if(!handled[ic])
                    {
                        deleteColumn(headDsat,nextDsat,previousDsat,satDegDsat[ic]-1,ic);
                        addColumn(headDsat,nextDsat,previousDsat,satDegDsat[ic],ic);
                    }
                                        
                }

                //increase the "colorNo" colored neighbor(s) in ordered graph of ic 
                colorTracker[colorNo][ic]++;
                //decrease ic's degree in unordered graph
                inducedDegDsat[ic] = inducedDegDsat[ic] - 1;
                
            }
        }
    }
    if(tagCol) delete[] tagCol;
}

//After removing a color of a colunm we look for all of its adjacent 
//columns and update their saturation degree where applicable
void IMatrix::satDegDec(int jcol,int order,int colorNo)
{
    //this array is for simple tagging. we select a nonzero element of jcol
    //if we process a neighboring column of jcol we dont need to process that
    //column in any next iteration. so we tag that column and don't procees that  
    bool *tagCol = new bool[N+1]();
    tagCol[jcol] = true;
    for (int jp = jpntr[jcol] ; jp < jpntr[jcol+1] ; jp++)
    {
        int ir = row_ind[jp];

        for( int ip = ipntr[ir] ; ip < ipntr[ir+1] ; ip++)
        {
            int ic = col_ind[ip];
            
            //column ic has not been processed
            if(!tagCol[ic])
            {
                //it is now tagged so will not be processed in any next iteration
                tagCol[ic] = true;
                //noumber of "colorNO" colored neighbors in ordered graph
                int prevColorCount = colorTracker[colorNo][ic];
                //one "colorNo" colored neighbor(which will be deleted) in ordered graph so we can decrease the saturation degree of ic     
                if(prevColorCount==1)
                {
                    satDegDsat[ic]--;
                    maxsatDsat = max(maxsatDsat,satDegDsat[ic]);
                    //this means the column isn't handled 
                    //and in the bucket so we can update its
                    //degree list (sat degree) in bucket 
                    if(!handled[ic])
                    {
                        deleteColumn(headDsat,nextDsat,previousDsat,satDegDsat[ic]+1,ic);
                        addColumn(headDsat,nextDsat,previousDsat,satDegDsat[ic],ic);
                    }
                }
                //decrease the "colorNo" colored neighbor(s) in ordered graph of ic 
                colorTracker[colorNo][ic]--;
                if(colorTracker[colorNo][ic]<0)
                    cout<<"Onnoooooo "<<ic<<" jcol "<<jcol<<" colorTracker[colorNo][ic] "<<colorTracker[colorNo][ic]<<endl; 
                //increase ic's degree in unordered graph
                inducedDegDsat[ic] = inducedDegDsat[ic] + 1;
            }
        }
    }
    if(tagCol) delete[] tagCol;
}

//this function takes a column, a color no and color matrix as input
//it check the colorNo is available for it or no by checking the
//colors of its neighbors
bool IMatrix::colorAvailable(int jcol, int colorNo)
{

    for (int jp = jpntr[jcol] ; jp < jpntr[jcol+1] ; jp++)
        {
            int ir = row_ind[jp];
             for( int ip = ipntr[ir] ; ip < ipntr[ir+1] ; ip++)
             {
                int ic = col_ind[ip];
                if(colorDsat[ic] == colorNo)
                    return false;
             }
        
    }

    return true;
}   

//main branch and bound coloring method for dsatur exact
int IMatrix::branchColor(int order,int colorBoundary)
{

    int updatedColoring,jcol,colorNo;
    subProblems++;
    
    if(colorBoundary >= UB) 
    {
        return colorBoundary;
    }
    if(UB<=LB)
    {
        return UB;
    }
    if(order >= N)
    {
        return colorBoundary;
    }

    //if time exceeds the limit we will stop and return
    if(( clock() - startTime ) / (double) CLOCKS_PER_SEC > CUT_OFF_TIME*60) 
    {
       return UB;
    }
    
    
    //find the vertex with maximum saturation degree. 
    //if there is more than one the we take the head to 
    //break the tie
    jcol = getColumn();
    handled[jcol] = true;
    

    //jcol is colored by all available colors (one at a time) and then 
    //recursively colors rest of the unordered columns to find a new coloring 
    for(colorNo=1; colorNo<=colorBoundary;colorNo++)
    {
        //from color 1 to colorboundary availability of those colors are checked for jcol
        if(colorAvailable(jcol,colorNo))    
        {
            tagDsat[jcol] = N;  
            
            //delete jcol from unorderd list 
            deleteColumn(headDsat,nextDsat,previousDsat,satDegDsat[jcol],jcol);
            
            colorDsat[jcol] = colorNo; // numord is new color for each clique member    
            //available color is>maxsat 
            if(colorNo>maxsatDsat)  
            {
                maxsatDsat++;
                colorTracker[maxsatDsat] = new int[N+1]();  
            }
            //increase saturation degrees of jcol's neighbors 
            satDegInc(jcol,order,colorDsat[jcol]);
            //recursuvely color the rest of the graph
            updatedColoring = branchColor(order+1,colorBoundary);
            //new improved coloring found
            if(updatedColoring<UB)
            {
                //this is our new decreased upperbound
                UB = updatedColoring;
                cout<<endl<<"Got current best coloring at time : "<<( clock() - startTime ) / (double) CLOCKS_PER_SEC<<endl;
                //coloring
                cout<<"New coloring: "<<UB<<" Subproblems: "<<subProblems<<endl;
		//uncomment this part to see which column gets which color	
                /*for(int i=1;i<=N;i++)
                    cout<<i<<"-"<<colorDsat[i]<<";";
                cout<<endl; */
            }
            
            //Remove color of jcol part

            //decrease saturation degrees of its neighbors
            satDegDec(jcol,order,colorDsat[jcol]);  
            //removing color from jcol
            colorDsat[jcol] = N;
            //add jcol back to unordered list
            addColumn(headDsat,nextDsat,previousDsat,satDegDsat[jcol],jcol);

            tagDsat[jcol] = 0;
            //new upperbound is <= color boundary so return from here
            if(UB <= colorBoundary)
            {
                handled[jcol] = false;
                return UB;
            }
        
        }   
    }
    //if colorBoundary+1 is still < the UB the we will color jcol
    //using colorBoundary+1 and recursively color the rest of the 
    //graoh to fin a new improved upper bound
    if(colorBoundary+1 < UB)
    {
        
        colorDsat[jcol] = colorBoundary+1;
        tagDsat[jcol] = N;  
        
        //delete jcol from unorderd list 
        deleteColumn(headDsat,nextDsat,previousDsat,satDegDsat[jcol],jcol);

        //color+1 is>maxsat 
        if(colorBoundary+1>maxsatDsat)  
        {
            maxsatDsat++;
            colorTracker[maxsatDsat] = new int[N+1]();  
        }
        //increase saturation degrees of jcol's neighbors 
        satDegInc(jcol,order,colorDsat[jcol]);

        //recursuvely color the rest of the graph using new color boundary
        updatedColoring = branchColor(order+1,colorBoundary+1);
        //new improved coloring found
        if(updatedColoring<UB)
        {
            //this is our new decreased upperbound
            UB = updatedColoring;
            cout<<endl<<"Got current best coloring at time : "<<( clock() - startTime ) / (double) CLOCKS_PER_SEC<<endl;
            //coloring
            cout<<"New coloring: "<<UB<<" Subproblems: "<<subProblems<<endl;
            //uncomment this part to see which column gets which color
	    /*for(int i=1;i<=N;i++)
                cout<<i<<"-"<<colorDsat[i]<<";";
            cout<<endl; */
        }
        //Remove color of jcol part

        //decrease saturation degrees of its neighbors
        satDegDec(jcol,order,colorDsat[jcol]);  
        //adding it back to uncolored list of columns   
        addColumn(headDsat,nextDsat,previousDsat,satDegDsat[jcol],jcol);
        //remove color
        colorDsat[jcol] = N;
        tagDsat[jcol] = 0; 
        
    }       
    handled[jcol] = false;
    return UB;
}



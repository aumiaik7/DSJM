#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sys/times.h>
#include <sys/resource.h>
#include <cstring>

#include <string>
#include <sstream>

#include <vector>

#include "Timer.h"

#include "driver.hh"

#include "CLI.h"
#include "gcolor.hh"
#include "mmio.h"
#include "Matrix.hh"
#include "Result.hh"



#define CLK_TCK 100
#define MAXLINE 100

class RunningTimeInfo
{
public:
    int color;
    double ordering_time;
    double coloring_time;
    int *ngrp;

    RunningTimeInfo() { }
    RunningTimeInfo(int c, double o_t, double c_t, int *grp) : color(c), ordering_time(o_t), coloring_time(c_t), ngrp(grp) { }
};

template <typename T>
void readData(Configuration *configuration, FILE* f, T *matrix)
{
    // As we are not going to use the value, we are simply using
    // this placeholder value so that we can read each line
    // from the input matrix.
    double value;


    for ( int i = 1,row,col ; i <= configuration->nz; i++ )
    {
        if (configuration->is_pattern)
        {
            fscanf (f, "%d %d\n", &row, &col);
        }
        else
        {
            fscanf (f, "%d %d %lg\n", &row, &col, &value);
        }

        matrix->setIndRowEntry(i,row);
        matrix->setIndColEntry(i,col);

        if(configuration->is_symmetric)
        {
            matrix->setIndRowEntry(i + configuration->nz,col);
            matrix->setIndColEntry(i + configuration->nz,row);
        }

    }

    if (f != stdin)
        fclose(f);
}

template <typename T>
void printStatistics(Configuration* configuration, T* matrix)
{
  cout << configuration->inputfile << " " << configuration->M << " & " << configuration->N << " & " << matrix->getNNZ() << " & " << matrix->getRowMax() << " & " << matrix->getRowAvg() << " & " << matrix->getRowMin() << " & " << matrix->getColMax() << " & " << matrix->getColAvg() << " & " << matrix->getColMin() << endl;
}

template <typename T, typename T_Timer>
RunningTimeInfo runOrderingAndColoringAlgorithm(Configuration* configuration, T* matrix,T_Timer &timer )
{

    int maxgrp;
    double ordering_time = -1;
    double coloring_time = -1;
    bool success = false;
    int *ngrp;



    if(configuration->is_cseg != true)
    {
        if (configuration->oMethod == CLI::SLO )
        {
            timer.Start();
            int *list;
            list = new int[configuration->N+1];
            ngrp = new int[configuration->N+1];
            success = matrix->slo(list);
            timer.Stop();

            ordering_time = timer.GetWallTime();

            timer.Start();
            maxgrp = matrix->greedycolor(list,ngrp);
            timer.Stop();

            coloring_time = timer.GetWallTime();
            delete[] list;

        }
        else if( configuration->oMethod == CLI::IDO )
        {
             timer.Start();
            int *list = new int[configuration->N+1];
            ngrp = new int[configuration->N+1];
            success = matrix->ido(list);
            timer.Stop();

            ordering_time = timer.GetWallTime();

            timer.Start();
            maxgrp = matrix->greedycolor(list,ngrp);
            timer.Stop();

            coloring_time = timer.GetWallTime();
            delete[] list;
        }
        else if ( configuration->oMethod == CLI::LFO)
        {
            timer.Start();
            int *list = new int[configuration->N+1];
            ngrp = new int[configuration->N+1];
            success = matrix->lfo(list);
            timer.Stop();

            ordering_time = timer.GetWallTime();

            timer.Start();
            maxgrp = matrix->greedycolor(list,ngrp);
            timer.Stop();

            coloring_time = timer.GetWallTime();

            delete[] list;
        }
        else if ( configuration->oMethod == CLI::RLF)
        {
            timer.Start();
            ngrp = new int[configuration->N+1];
            maxgrp =  matrix->rlf(ngrp);
            timer.Stop();

            coloring_time = timer.GetWallTime();
        }
        else if ( configuration->oMethod == CLI::SDO)
        {
            /*if(configuration->is_cseg )
            {
                cout << "DSatur method has not been implemented for CSegMatrix Yet!" << endl;
                exit(1);
            }*/
	          timer.Start();
            ngrp = new int[configuration->N+1];
            Matrix *nmatrix = dynamic_cast<Matrix*>(matrix);
            maxgrp = nmatrix->sdo(ngrp);
            timer.Stop();

            // maxgrp = matrix->seq();

            coloring_time = timer.GetWallTime();
        }
	else if(configuration->oMethod == CLI::SLO_RLF)
	{
 		timer.Start();int *list;
            	list = new int[configuration->N+1];
		ngrp = new int[configuration->N+1];
		success = matrix->slo(list);
		maxgrp =  matrix->slo_rlf(list,ngrp);
            	timer.Stop();
			
	}
	else if ( configuration->oMethod == CLI::EXACT)
        {
            //maxgrp is either no of columns or coloring we get from any heuristic
	          maxgrp = configuration->N;
            
            //We calculate lower bound in 2 ways. 
            //No-1: rho_max is lower bound. It is calculated in computedegree method
            //No-2: lower bound is calculated using slo_exact method the following code
            //block is used for that
            timer.Start();
            int *list;
            list = new int[configuration->N+1];
            // ngrp = new int[configuration->N+1];
            int *clique = new int[configuration->N+1]();
            success = matrix->slo_exact(list,clique);
            timer.Stop();

          
            ordering_time = timer.GetWallTime();

	         // here we set an Upper bound from the coloring we get from ido and coloring	
	         /*timer.Start();
            int *list = new int[configuration->N+1];
	          int *clique = new int[configuration->N+1];	
            ngrp = new int[configuration->N+1];
            success = matrix->idoDsatur(list,clique);
            timer.Stop();

            ordering_time = timer.GetWallTime();

            timer.Start();
            maxgrp = matrix->greedycolor(list,ngrp);
            timer.Stop();
	         */ //uncomment this part to set ido ordering+coloring as upperbound	

	          timer.Start();
            ngrp = new int[configuration->N+1];
            Matrix *nmatrix = dynamic_cast<Matrix*>(matrix);
            maxgrp = nmatrix->exact(maxgrp,clique,1,configuration->tieBrkDsat); // Lower bound from rho_max 
            //maxgrp = nmatrix->exact(maxgrp,clique,2,configuration->tieBrkDsat); // Lower bound we get from slo_exact method     
            timer.Stop();	

	    	
            coloring_time = timer.GetWallTime();
	     	
            //delete[] list;
        }
        else
        {
            cerr << "Unknown Method for Ordering" << endl;
            exit(EXIT_FAILURE);
        }
    }
    return RunningTimeInfo(maxgrp, ordering_time, coloring_time, ngrp);
}

template <typename T>
void writeNormalColoring(Configuration* configuration, T* matrix, Result result, int *ngrp)
{
    FILE *cfp;
    if((cfp = fopen("coloring.txt","w")) == NULL)
    {
        fprintf(stderr,"We could not open the coloring file, so exiting\n");
        exit(1);
    }

    fprintf(cfp,"%d\n", result.totalColors);
    for (int i = 1 ; i < configuration->N+1 ; i++)
    {
        fprintf(cfp,"%d %d\n",i,ngrp[i]);
        // assert(matrix->getNgrpEntry(i) != configuration->N ); This
        // assertion is not right, when we have to color each column
        // with its own color. So, I am removing it. But we had a
        // glitch in the verification system from the very
        // beginning. We do not check, whether we have colored every
        // one of them, or not.
    }
    fclose(cfp);
}

int gcolor_getrusage(struct rusage &runtime)
{
  getrusage (RUSAGE_SELF, &runtime);
  return (100 * runtime.ru_utime.tv_sec) + (runtime.ru_utime.tv_usec / 10000);
}

template <typename T>
void writeGraphDescription(Configuration *configuration, T *matrix)
{
    int *tag;
    tag = new int[configuration->N+1];
    FILE *gfp = fopen("graph.description","w");
    if(gfp == NULL)
    {
        fprintf(stderr,"File Error for graph.description\n");
        exit(1);
    }

    fprintf(gfp,"%d\n",configuration->N);

    for (int i = 1; i <= configuration->N; i++)
    {
        tag[i] = 0;
    }

    for (int jcol = 1; jcol <= configuration->N; jcol++)
    {
        tag[jcol] = configuration->N;
        for (int jp = matrix->getJpntrEntry(jcol) ; jp <= matrix->getJpntrEntry(jcol+1) - 1; jp++)
        {
            int ir = matrix->getIndRowEntry(jp);
            for( int ip = matrix->getIpntrEntry(ir) ; ip <= matrix->getIpntrEntry(ir+1) -1; ip++)
            {
                int ic = matrix->getIndColEntry(ip);

                if(tag[ic] < jcol)
                {
                    tag[ic] = jcol;
                    fprintf(gfp,"%d %d\n",jcol,ic);
                }
            }
        }
    }
    delete[] tag;

    fclose(gfp);
}

bool canReadPartitionFile(Configuration* configuration)
{
  FILE *pf;
  if(( pf = fopen(configuration->partitionFile, "r")) == NULL)
    {
      return false;
    }
  else
    {
      fclose(pf);
      return true;
    }
}

bool canReadFile(Configuration *configuration)
{
  FILE *f;
  if( configuration->inputfile_given == false)
    {
      fprintf (stderr, "Usage: %s -i [martix-market-filename] -m [slo/ido/rlf/lfo/sdo/exact]\n", configuration->exeName.c_str());
      fprintf (stderr, "Use an input file with -i $filename\n");
      return false;
    }
  if ((f = fopen (configuration->inputfile, "r")) == NULL)
    {
      cerr << configuration->inputfile << " cannot be found " << endl;
      return false;
    }
  else
    {
      fclose(f);
    }
  return true;
}

FILE* readBanner(Configuration* configuration)
{
  /**
   * ----------------------------------------------------------------------
   * Let's Try to Read the Banner of the Matrix Market Format
   * FIle.
   *
   * */

  /**
   * Matrix Market or mmio.h Specific Data Types' Declaration
   * ----------------------------------------------------------------------
   **/
  MM_typecode matcode;
  int ret_code;

  FILE *f;
  f  = fopen(configuration->inputfile, "r");
  if (mm_read_banner (f, &matcode) != 0)
    {
      fprintf (stderr,"%s -> Could not process Matrix Market banner.\n",configuration->inputfile);
      exit (1);
    }

  /**
   * ----------------------------------------------------------------------
   * This is how one can screen matrix types if their applicaiton
   * only supports a subset of the Matrix Market data types.
   **/

  if (mm_is_complex (matcode) && mm_is_matrix (matcode) &&
      mm_is_sparse (matcode))
    {
      printf ("Sorry, this application does not support ");
      printf ("Market Market type: [%s]\n", mm_typecode_to_str (matcode));
      exit (1);
    }

  /**
  *
   * Find out the size of the sparse matrix
   **/

  if ((ret_code = mm_read_mtx_crd_size (f, &(configuration->M), &(configuration->N), &(configuration->nz))) != 0)
    exit (1);


  configuration->NZ = configuration->nz;

  configuration->is_symmetric = mm_is_symmetric(matcode);
  configuration->is_pattern = mm_is_pattern(matcode);
  if(configuration->is_symmetric)
    configuration->NZ = 2 * configuration->nz;

  return f;
}


vector<string> split(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

int run(Configuration *configuration)
{
    /*
     *
     *  Declaration Block
     *  ----------------------------------------
     *
     */

    FILE *pf;
    char line[MAXLINE+1];

    struct rusage runtime;
    int t1, t2,t3,t4;


    /**
     * ----------------------------------------------------------------------
     * In this portion we are going to handle whether we can
     * open a file and read its input.
     * So we check the following things,
     *      1. We check whether we are given at least two
     *         arguments in our console or not. The second
     *         one obviously the name of the file, while
     *         the first one becoming the name of our program.
     *
     *      2. Then we try to open a handle for the file. If
     *         we are successfull finding the file and has,
     *         then we proceed on next iteration.
     * ----------------------------------------------------------------------
     */

    Result result;
    result.configuration = configuration;




    if( !canReadFile(configuration) )
    {
        std::cerr << "Error, File cannot be found!" << std::endl;
        exit(EXIT_FAILURE);
    }

    result.inputFile = configuration->inputfile;
    result.is_Cseg = configuration->is_cseg;
    result.method = configuration->oMethod;
    result.printInfo();

    FILE *f = readBanner(configuration);



    /**
     * Now we need to allocate some memory for the things we are going to read
     * ----------------------------------------------------------------------
     */
    int *tag;
    IMatrix *matrix = NULL;
    try
    {


    	Matrix* _matrix = new Matrix(configuration->M,configuration->N,configuration->NZ, false);
		matrix = _matrix;
		matrix->setVerify(configuration->verify);


        /** Time Related Variables -- CPack */
        double ordering_time = -1;
        double coloring_time = -1;

        ColPack::Timer timer;



        // This is the point where we have taken the input from the mtx
        // file.
        readData(configuration,f,matrix);



        //timer.Start();
        matrix->computeCCS();


        int nnz = matrix->compress();

        result.matrix_M = configuration->M;
        result.matrix_N = configuration->N;
        result.matrix_NNZ = nnz;

        if(!configuration->io_at_end)
            result.printMatrixInfo();


        matrix->computeCRS();



        if ( configuration->statistics_only == true)
        {
            // Print the information related to the Matrix and exit.
            printStatistics(configuration, matrix);
            return 0;
        }
 
      
        matrix->computedegree();
        // maxgrp, coloring_time, ordering_time, kind of a tuple.
        int maxgrp;
        RunningTimeInfo runningTimeInfo = runOrderingAndColoringAlgorithm(configuration, matrix,timer);
        maxgrp = runningTimeInfo.color;
        ordering_time = runningTimeInfo.ordering_time;
        coloring_time = runningTimeInfo.coloring_time;


        result.totalColors = maxgrp;
        result.ordering_time = ordering_time;
        result.coloring_time = coloring_time;

	//cout<<configuration->inputfile<<","<<configuration->M<<","<<configuration->printInfo<<","<<nnz<<",";

        if(configuration->io_at_end)
            result.printMatrixInfo();
        result.printResult();

        // Verification Block
        // ------------------------------------------------------------
        if(configuration->verify)
        {
            writeGraphDescription(configuration,matrix);

            writeNormalColoring(configuration, matrix,result, runningTimeInfo.ngrp);
        }
        // Verification Block Ends
        // ------------------------------------------------------------
        if(matrix)
            delete matrix;

        if(configuration->inputfile)
            delete[] configuration->inputfile;
    }

    catch(std::bad_alloc)
    {
        std::cerr << "Memory Exhausteddd\n";
        delete[] tag;

        if(matrix)
            delete matrix;

        exit(1);
    }
}

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <iostream>
#include "mmio.h"
#include "Matrix.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "CLI.h"
#include <cstring>

using namespace std;


/**
 * The following struct is going to capture all the command line parameters
 * we are supplying 
 */ 

struct globalArgs_t
{
  CLI::ordering_method oMethod; /* Method to use -- > LFO, RLF, IDO, SLO */
  char *inputfile; /* Matrix Description in Matrix Market format */ 

}globalArgs; 


/**
 * 
 * The following optString is an important part of the Command Line Interface. 
 * -m lfo -h -? -i inputfile.mtx
 */ 
static const char *optString = "m:i:h?";


/**
 * This function is called when the application is called with the display help
 * command Line options. 
 */ 
void display_usage(void) 
{ 
  puts("Graph Coloring Driver Program"); 
  /* ... */ 
  exit(EXIT_FAILURE); 
}


int main (int argc, char *argv[])
{
  int opt  = 0; 

  /**
   * Initialize the globalArgs before we get to work. */ 
  globalArgs.oMethod = CLI::SLO; 
  
  /* We are starting parsing the command Line arguments with getopt 
   * method 
   */ 

  opt = getopt(argc,argv, optString); 
  while ( opt != -1 ) 
    { 
      switch( opt ) 
	{ 
	case 'm' : 
	  // Select the method for ordering 
	  // between SLO, IDO and LFO
	  // RLF is handled by a separate program. 
	  
	  // But we might merge them in a single program. in near future. 
	  // globalArgs.oMethod = optarg; 
	  globalArgs.oMethod = CLI::parseMethod(optarg);
	  if(globalArgs.oMethod == CLI::LFO) 
	    // cout << "You have selected LFO" << endl; 
	    ;
	  else if( globalArgs.oMethod == CLI::SLO ) 
	    // cout << "You have selected SLO" << endl; 
	    ;
	  else if ( globalArgs.oMethod == CLI::IDO ) 
	    // cout << "You have selected IDO " << endl; 
	    ;
	  else if( globalArgs.oMethod == CLI::RLF  )
	    // cout << "You have selected RLF" << endl; 
	    ;
	  else 
	    cout << "Unknown method selected for ordering " << endl; 
	  break; 
	case 'i': 
	  /** 
	   * \todo Should we change the struct globalArgs to a class 
	   * so that it can better serve in case dealing with memory handling ? 
	   * So we got something with the -i option.  
	   * This is a character string and most probably the name 
	   * of the file containing the matrix description in matrix
	   * market format. We should store this in our globalArgs struct. 
	   */ 
	  if(optarg != NULL && strlen(optarg ) > 0) 
	    {
	      delete[] globalArgs.inputfile ; 
	      globalArgs.inputfile = new char[strlen(optarg) + 1]; 
	      strcpy(globalArgs.inputfile, optarg); 
	      // cout << "We have copied it to the globalArgs - " <<
	      // globalArgs.inputfile << endl;
	    } 
	  	   
	  break;
	case 'h': /* fall-through is intentional */ 
	case '?': 
	  display_usage();
	  break; 
	default: 
	  break;
	}

      opt = getopt(argc,argv,optString); 
    }

  /* 
   * I think we pretty much gathered all the informations from the
   * command line as you can see, now we should make changes to the
   * code so that it actually honors the code presented here. 
   */ 

 

	  
  
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, nz;
  int t1, t2,t3,t4;


  // int *indRow,*jpntr,*indCol, *ipntr, maxClique, *iwa;

  double *val;

  if (argc < 2)
    {
      fprintf (stderr, "Usage: %s -i [martix-market-filename] -m [lfo/ido/slo/rlf]\n", argv[0]);
      exit (1);
    }
  else
    {
      // if ((f = fopen (argv[1], "r")) == NULL)
      if ((f = fopen (globalArgs.inputfile, "r")) == NULL)	
	exit (1);
    }

  if (mm_read_banner (f, &matcode) != 0)
    {
      printf ("Could not process Matrix Market banner.\n");
      exit (1);
    }


  /* This is how one can screen matrix types if their application */
  /* only supports a subset of the Matrix Market data types.  */

  if (mm_is_complex (matcode) && mm_is_matrix (matcode) &&
      mm_is_sparse (matcode))
    {
      printf ("Sorry, this application does not support ");
      printf ("Market Market type: [%s]\n", mm_typecode_to_str (matcode));
      exit (1);
    }

  /* find out size of sparse matrix .... */

  if ((ret_code = mm_read_mtx_crd_size (f, &M, &N, &nz)) != 0)
    exit (1);
  /**
   * We are actually starting the block where we should or 
   * we can mimick subroutine DSM 
   */




  /* reseve memory for matrices */


  /** 
  std::cerr << sizeof (int);
  std::cerr << INT_MIN << ":" << INT_MAX << "\n";
  std::cerr << "(M,N) (" << M << "," << N << ")\n";
  */

  /*  if ((indRow = (int *) malloc((nz+1) * sizeof(int))) == NULL)
      exit(1);
      if ((indCol = (int *) malloc((nz+1) * sizeof(int))) == NULL)
      exit(1);
      if ((jpntr = (int *) malloc((N+2) * sizeof(int))) == NULL)
      exit(1);
      if ((ipntr = (int *) malloc((M+2) * sizeof(int))) == NULL)
      exit(1);
      if ((iwa = (int *) malloc((6*N+1) * sizeof(int))) == NULL)
      exit(1);
  */

  int maxClique;
  try
    {
      int *iwa = new int[6 * N + 1];

      val = (double *) malloc ((nz + 1) * sizeof (double));
      maxClique = nz;      // \todo isnt maxclq is an output of slo and ido.
      /** \todo
       * I just allocated the objects, 
       * But didn't do the actually reading of the datas. 
       *
       * \todo 
       * Do the same task as the above for Fortran Code. 
       *
       * Declare the object matrix; 
       **/
      Matrix matrix(M,N,nz);

      for (int i = 1,row,col; i <= nz; i++)
	{
	  fscanf (f, "%d %d %lg\n", &row, &col, val + i);

	  matrix.setIndRowEntry(i,row); 
	  matrix.setIndColEntry(i,col);
	  /** printf("%d %d %lg\n",*(indRow+i),*(indCol+i),val+i);
	      if (i == 10) 
	      exit(0);
	  **/
	}
      if (f != stdin)
	fclose (f);

      //  for(int i = 1 ; i <= nz; i++)
      //std::cout << "(" << *(indRow + i)  << "," << *(indCol + i)  << ")\n";


      
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // ,----										     //
    // | Start Counting Time. So We do not want any further output until we stop counting time.  //
    // `----										     //
    ///////////////////////////////////////////////////////////////////////////////////////////////


      struct rusage runtime;
      getrusage (RUSAGE_SELF, &runtime);
      t1 = (100 * runtime.ru_utime.tv_sec) + (runtime.ru_utime.tv_usec / 10000);

      matrix.srtdat();

      /** Complress the data and determine the number of 
       * non-zero elements of A
       */

      /*

      int nnz = nz; 
      for (int i = 1; i <= M; i++)
	{
	  iwa[i] = 0;
	}
      nnz = 1;
      int k;
      int ir;
      for (int j = 1; j <= N; ++j)
	{
	  k = nnz;
	  for (int jp = jpntr[j]; jp <=(jpntr[j + 1] - 1); ++jp)
	    {
	      ir = indRow[jp];
	      if (iwa[ir] != j)
		{
		  indRow[nnz] = ir;
		  ++nnz;
		  iwa[ir] = j;
		}
	    }
	  jpntr[j] = k;
	}
      jpntr[N + 1] = nnz;

      */    

      int nnz = matrix.compress();

      matrix.setr();

      getrusage (RUSAGE_SELF, &runtime);
      t2 = (100 * runtime.ru_utime.tv_sec) + (runtime.ru_utime.tv_usec / 10000);


      /*       _\|/_
               (o o)
       +----oOO-{_}-OOo---------------------------------------------+
       |                                                            |
       |Write the graph informations on a file, so that we can run a|
       |verifier on the coloring informations.                      |
       +-----------------------------------------------------------*/
      


      FILE *gfp = fopen("graph.description","w"); 

      int *tag = new int[N+1];
      if(gfp == NULL) 
	{ 
	  fprintf(stderr,"File Error for graph.description\n"); 
	  exit(1); 
	} 

      fprintf(gfp,"%d\n",N);

      for (int i = 1; i <= N; i++)
	{
	  tag[i] = 0; 
	}

      int s = 1;
      for (int jcol = 1; jcol <= N; jcol++,s++) 
	{
	  tag[jcol] = N; 
	  for (int jp = matrix.getJpntrEntry(jcol) ; jp <= matrix.getJpntrEntry(jcol+1) - 1; jp++)
	    { 
	      int ir = matrix.getIndRowEntry(jp); 
	      for( int ip = matrix.getIpntrEntry(ir) ; ip <= matrix.getIpntrEntry(ir+1) -1; ip++) 
		{
		  int ic = matrix.getIndColEntry(ip); 

		  if(tag[ic] < jcol) 
		    {
		      tag[ic] = jcol; 
		      fprintf(gfp,"%d %d\n",jcol,ic);
		    } 
		} 
	    } 
	} 

      fclose(gfp); 

      getrusage (RUSAGE_SELF, &runtime);
      t3 = (100 * runtime.ru_utime.tv_sec) + (runtime.ru_utime.tv_usec / 10000);


    


      matrix.degr ();


      /**
       * Following commented line used to verify the codes. Kind of test codes. 
       */ 
      /* FILE *fp = fopen("Application.out","w");
	 if (fp == NULL) 
	 exit(1); 
	 fprintf(fp,"%d %d %d %d\n",M,N,nz,nnz);
	 for(int i = 1; i <= nnz ;i++)
	 {
	 fprintf(fp,"(%d %d)\n",*(indRow+i),*(indCol+i));
	 }

	 for (int i = 1; i <= N + 1 ; i++) 
	 {
	 fprintf(fp,"%d\n",*(jpntr+i));
	 }

	 for (int i = 1; i <= M + 1 ; i++) 
	 {
	 fprintf(fp,"%d\n",*(ipntr+i));
	 }

	 for (int i = 1; i <= N ; i++) 
	 {
	 fprintf(fp,"%d\n",*(iwa + 5*N + i));
	 }
	 fclose(fp);
  
	 return 0;
      */


        /////////////////////////
        // ,----	       //
        // | Find out ordering //
        // `----	       //
        /////////////////////////

      /** \class Matrix
       * \todo
       * We have to use some kind of code generation or any other framework so that, 
       * we can have three different executables or we can have command line parameters 
       * to have three different functionalities from a single executable 
       */ 

      /** \class Matrix
       * \todo
       * We have to make lots of unit testing so that we can make sure that all the time 
       * that the coded procedures are still working in the changing world and changing 
       * changing requirements. 
       */

      if (globalArgs.oMethod == CLI::SLO ) 
	{ 
	  matrix.slo(); 
	}
      else if( globalArgs.oMethod == CLI::IDO ) 
	{
	matrix.ido(); 
	}
      else if ( globalArgs.oMethod == CLI::LFO) 
	{
	matrix.lfo();
	}
      else if ( globalArgs.oMethod == CLI::RLF) 
	{	
	// matrix.rlf(); 
	  cerr << "Not Implemented in this Executable, use ./RecursiveLargestFirst [matrix-file]" << endl;
	  return 1;
	}
      else
	{
	cerr << "Unknown Method for Ordering" << endl; 
	return 1;
	}
      /* matrix.slo (); */ 


      /* matrix.ido(); */

      /* matrix.lfo(); */ 
      

      int maxgrp;

      /***************************/
      /*  D o   c o l o r i n g  */
      /***************************/

      maxgrp = matrix.seq();
      //matrix.seq(N, indRow,jpntr,indCol,ipntr,iwa+5*N,iwa +2 * N ,&maxgrp,iwa+1);
        
      /*
       *  Time Count Stopped
       */

      getrusage (RUSAGE_SELF, &runtime);
      t4 = (100 * runtime.ru_utime.tv_sec) + (runtime.ru_utime.tv_usec / 10000);



      /* We are calculating the time taken for the whole coloring algorithm. */


      // cerr<<"Ordering done. " << endl;
      cout<<globalArgs.inputfile  << "\t" << globalArgs.oMethod<< "\tK = " << maxgrp<<  " Took "<<(t4-t1-(t3-t2))/100<<"."<<((t4-t1-(t3-t2))%100)/10<<" seconds"<<endl;


      /** 
       * 
       * We are going to stop the output of the whole coloring to the 
       * standard for a while
       
       * 
       */


      
      /* for (int i = 1; i < N + 1; i++)
	{
	  std::cout << " LIST; " << i << "; " << *(iwa + 4 * N + i) << "\n";
	}

      for ( int i = 1; i < N + 1 ; i++) 
	{ 
	  std::cout << "COLORING; " << i << "; " <<  *(iwa + 2*N + i )  << "\n" ; 
	} 
      */

      



      /*
                __   _,--="=--,_   __
               /  \."    .-.    "./  \
              /  ,/  _   : :   _  \/` \
              \  `| /o\  :_:  /o\ |\__/
               `-'| :="~` _ `~"=: |
                  \`     (_)     `/
           .-"-.   \      |      /   .-"-.
      .---{     }--|  /,.-'-.,\  |--{     }---.
       )  (_)_)_)  \_/`~-===-~`\_/  (_(_(_)  (
      (  Write Coloring output to file        )
       )                                     (
      '---------------------------------------'
      */


      FILE *fp = fopen("coloring.txt","w"); 
      if (fp == NULL) 
	{
	  std::cout << "File Error for coloring informations\n"; 
	  exit(1);
	}

      for ( int i = 1; i < N + 1 ; i++) 
	{ 
	  // fprintf(fp,"%d %d\n",i , *(iwa + 2*N + i ));
	  fprintf(fp,"%d %d\n",i , matrix.getNgrpEntry(i));
	} 

      fclose(fp);


    
    }
  catch (std::bad_alloc)
    {
      std::cerr << "Memory Exhausted\n";
      exit (1);
    }

  return 0;
}
 

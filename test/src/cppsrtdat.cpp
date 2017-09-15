#include <cstdio>
#include <cstdlib>
#include <memory>
#include <iostream>
#include "mmio.h"
#include "Matrix.h"
#include <climits>

int main (int argc, char *argv[])
{
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, nz;   
  // int *indRow,*jpntr,*indCol, *ipntr, maxClique, *iwa;

  double *val;

  if (argc < 2)
    {
      fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
      exit(1);
    }
  else    
    { 
      if ((f = fopen(argv[1], "r")) == NULL) 
	exit(1);
    }

  if (mm_read_banner(f, &matcode) != 0)
    {
      printf("Could not process Matrix Market banner.\n");
      exit(1);
    }


  /* This is how one can screen matrix types if their application */
  /* only supports a subset of the Matrix Market data types.  */

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
      mm_is_sparse(matcode) )
    {
      printf("Sorry, this application does not support ");
      printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
      exit(1);
    }

  /* find out size of sparse matrix .... */

  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
    exit(1);
  /**
   * We are actually starting the block where we should or 
   * we can mimick subroutine DSM 
   */

  /* reseve memory for matrices */
  std::cout << sizeof(int);
  std::cout << INT_MIN << ":" << INT_MAX << "\n";
  std::cout << "(M,N) (" << M << "," << N << ")\n";

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
      int *indRow =new int[nz+1]; 
      int *indCol = new int[nz+1]; 
      int *jpntr = new int[N+2]; 
      int *ipntr = new int[N+2]; 
      int *iwa( new int[6*N+1]); 

      val = (double *) malloc((nz+1) * sizeof(double));
      maxClique = nz; // \todo isnt maxclq is an output of slo and ido.
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

      for (int i = 1; i <= nz; i++)
	{
	  fscanf(f, "%d %d %lg\n",indRow+i, indCol+i,val+i);
	  /** printf("%d %d %lg\n",*(indRow+i),*(indCol+i),val+i);
	      if (i == 10) 
	      exit(0);
	  **/
	}
      if (f !=stdin) fclose(f);

      //  for(int i = 1 ; i <= nz; i++)
      //std::cout << "(" << *(indRow + i)  << "," << *(indCol + i)  << ")\n";

      matrix.srtdat();

      int nnz = nz;
      /** Complress the data and determine the number of 
       * non-zero elements of A
       */ 
      for (int i = 1 ; i <= M; i++)
	{
	  iwa[i] = 0; 
	}
      nnz = 1; 
      int k;
      int ir;
      for (int j=1 ; j <= N; ++j)
	{
	  k = nnz; 
	  for (int jp = jpntr[j]; jp <= (jpntr[j+1] - 1); ++jp)
	    {
	      ir = indRow[jp]; 
	      if(iwa[ir] != j)
		{
		  indRow[nnz] = ir; 
		  ++nnz; 
		  iwa[ir] = j;
		}
	    }
	  jpntr[j] = k; 
	}
      jpntr[N+1] = nnz;
      std::cerr << "NNZ " << nnz << "\n";

      matrix.setr();
      matrix.degr();

      *(indRow+nnz) = *(indCol + nnz ) = 1;

      matrix.slo();

      
      for (int i = 1; i <= nnz ; i++) 
	std::cout << " " << *(indRow+i) << ", " << *(indCol+i)<< "\n";
      for (int i = 1 ; i <= N+1 ; i++)
       	std::cout << " " << *(jpntr+i) << ", "<< *(ipntr+i)<< ", " << *(iwa+5*N+1+i)<< "\n";
      // for (int i = 1 ; i <= N+1 ; i++)
      // 	std::cout << " " << *(jpntr+i) << ", "<< *(ipntr+i)<< "\n" ; // ", " << *(iwa+5*N+1+i)<< "\n";
    }
  catch(std::bad_alloc)
    {
      std::cerr << "Memory Exhausted\n";
      exit(1);
    }

  return 0;
}

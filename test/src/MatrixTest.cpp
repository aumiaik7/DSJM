#include "MatrixTest.h"
#include <iostream>
#include "mmio.h"
#include <unistd.h>
#include <cstdio>
#include <cstdlib>

using namespace std;


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MatrixTest );

// void GetCurrentPath( char* buffer);

void
MatrixTest::setUp()
{
  /*
    if(log4cpp::Appender::getAppender("FileAppender") == NULL)
    appender = new log4cpp::FileAppender("FileAppender", "logs/test.log");
    else
    appender = log4cpp::Appender::getAppender("FileAppender");
    appender->setLayout(new log4cpp::BasicLayout());
    // log4cpp::Category& log = log4cpp::Category::getInstance(std::string("log"));
    log.setAdditivity(false);
    log.setAppender(appender);

    log.setPriority(log4cpp::Priority::INFO);
  */


  cout << "CSegMatrix::setUp() " << endl;
  // char PATH[1000];
  // GetCurrentPath(PATH);
  // cout << "CurrentPath " << PATH << endl;

  char fileName[] = "abb313.mtx" ;
  MM_typecode matcode;

  FILE *f;
  if (( f = fopen(fileName, "r")) == NULL)
    {
      CPPUNIT_FAIL("File Not Found");
    }

  if (mm_read_banner (f,&matcode) != 0)
    {
      CPPUNIT_FAIL("Could not processs Matrix Market banner.\n");
    }

  int M,N,nz;
  int ret_code;
  if ( (ret_code = mm_read_mtx_crd_size(f,&M,&N,&nz)) != 0)
    CPPUNIT_FAIL("I don't know why");
  CPPUNIT_ASSERT_EQUAL(313,M);
  CPPUNIT_ASSERT_EQUAL(176,N);
  CPPUNIT_ASSERT_EQUAL(1557,nz);


  int NZ = nz;

  if(mm_is_symmetric(matcode))
    NZ = 2 * nz;

  matrix = new Matrix(M,N,NZ,false);

  double value;


  // We are reading co-ordinates from the files.
  for ( int i = 1, row, col ; i <= nz; i++ )
    {
      if ( mm_is_pattern(matcode))
	{
	  fscanf(f, "%d %d\n", &row, &col);
	  // cout << row << " " << col << endl;
	}
      else
	{
	  fscanf(f, "%d %d %lg\n", &row, &col, &value);
	}

      matrix->setIndRowEntry(i,row);
      matrix->setIndColEntry(i,col);

      if(mm_is_symmetric(matcode))
	{
	  matrix->setIndRowEntry(i + nz,col);
	  matrix->setIndColEntry(i + nz, row);
	}

    }

  // File IS Being CLosed Here
  if ( f != stdin )
    fclose(f);

  matrix->srtdat();
  matrix->compress();
  matrix->setr();

}

void MatrixTest::tearDown()
{
  cout << "CSegMatrixTest::tearDown() " << endl;
  delete matrix;
}




void
MatrixTest::testConstructor()
{
  // CPPUNIT_FAIL( "not Implemented" );
}

void
MatrixTest::testSLO()
{
  // Process
  int n;
  int indRow[11] = {0,1,2,1,3,2,4,3,4,5,5};
  int jpntr[7] = {0,1,3,5,7,10,11};
  int indCol[11] = {0,1,2,1,3,2,4,3,4,4,5};
  int ipntr[7] = { 0,1,3,5,7,9,11};
  int ndeg[6] = {0,2,2,2,3,1};
  int *list;
  int maxClique;
  int *iwa1;
  int *iwa2;
  int *iwa3;
  int *iwa4;


  n =5;
  int nnz = 10;
  int m = 5;

  Matrix matrix(m,n,nnz,false);


  // indRow = new int[nnz];
  // jpntr = new int[n+1];
  //indCol = new int[nnz];
  // ipntr = new int[m+1];
  // ndeg = new int[n];
  list = new int[n+1];
  maxClique = 10;
  iwa1 = new int[n+1];
  iwa2 = new int[n+1];
  iwa3 = new int[n+1];
  iwa4 = new int[n+1];
  // matrix.slo();

  // Check
  CPPUNIT_ASSERT(true);
  /* for(int i = 0 ;i < 6 ; i++)
    {
      std::cout << "list[" << i << "] = " << list[i] << "\n";
    }
  */

}

void
MatrixTest::testMin()
{
    // Set up
    Matrix m(5,5,5,false);
    int s,t;
    int result;

    // Process & Check

    // give s and t same. answer should be same to s.
    s = t = 4;
    result = m.min(s,t);
    CPPUNIT_ASSERT_EQUAL(result,s);

    // give s larger. answer should be t.
    s = 5;
    t = 4;
    result = m.min(s,t);
    CPPUNIT_ASSERT_EQUAL(result,t);

    // give s smaller. answer should be t.
    s = 3;
    t = 4;
    result = m.min(s,t);
    CPPUNIT_ASSERT_EQUAL(result,s);
}

void
MatrixTest::testMax()
{
    // Set up
    Matrix m(5,5,5,false);
    int s,t;
    int result;

    // Process & Check
    // give s and t same. answer should be same to s.
    s = t = 4;
    result = m.max(s,t);
    CPPUNIT_ASSERT_EQUAL(result,s);

    // give s larger. answer should be t.
    s = 5;
    t = 4;
    result = m.max(s,t);
    CPPUNIT_ASSERT_EQUAL(result,s);

    // give s smaller. answer should be t.
    s = 3;
    t = 4;
    result = m.max(s,t);
    CPPUNIT_ASSERT_EQUAL(result,t);
}



void
MatrixTest::testIDO()
{

  // Process
  int n;
  int indRow[11] = {0,1,2,1,3,2,4,3,4,5,5};
  int jpntr[7] = {0,1,3,5,7,10,11};
  int indCol[11] = {0,1,2,1,3,2,4,3,4,4,5};
  int ipntr[7] = { 0,1,3,5,7,9,11};
  int ndeg[6] = {0,2,2,2,3,1};
  int *list;
  int maxClique;
  int *iwa1;
  int *iwa2;
  int *iwa3;
  int *iwa4;


  n =5;
  int nnz = 10;
  int m = 5;

  // Set up
  Matrix matrix(m,n,nnz,false);

  // indRow = new int[nnz];
  // jpntr = new int[n+1];
  //indCol = new int[nnz];
  // ipntr = new int[m+1];
  // ndeg = new int[n];
  list = new int[n+1];
  maxClique = 10;
  iwa1 = new int[n+1];
  iwa2 = new int[n+1];
  iwa3 = new int[n+1];
  iwa4 = new int[n+1];
  // matrix.ido();

  // Check
  // CPPUNIT_ASSERT(true);
//   for(int i = 0 ;i < 6 ; i++)
//     {
//       std::cout << "list[" << i << "] = " << list[i] << "\n";
//     }

}

void
MatrixTest::testNumSort()
{

  // Process
  int n;
  int indRow[11] = {0,1,2,1,3,2,4,3,4,5,5};
  int jpntr[7] = {0,1,3,5,7,10,11};
  int indCol[11] = {0,1,2,1,3,2,4,3,4,4,5};
  int ipntr[7] = { 0,1,3,5,7,9,11};
  int ndeg[6] = {0,2,2,2,3,1};
  int *list;
  int maxClique;
  int *iwa1;
  int *iwa2;
  int *iwa3;
  int *iwa4;


  n =5;
  int nnz = 10;
  int m = 5;

  // Set up
  Matrix matrix(m,n,nnz,false);

  // indRow = new int[nnz];
  // jpntr = new int[n+1];
  //indCol = new int[nnz];
  // ipntr = new int[m+1];
  // ndeg = new int[n];
  list = new int[n+1];
  maxClique = 10;
  iwa1 = new int[n+1];
  iwa2 = new int[n+1];
  iwa3 = new int[n+1];
  iwa4 = new int[n+1];


  // matrix.numsrt(n,n-1,ndeg,-1,iwa4,iwa2,iwa3);

  for(int i = 1 ;i < 6 ; i++)
    {
      std::cout << "ndeg[" << i << "] = " << ndeg[i] << "\n";
    }

  for(int i = 1 ;i < 6 ; i++)
    {
      std::cout << "iwa4[" << i << "] = " << iwa4[i] << "\n";
    }

  for(int i = 1 ;i < 6 ; i++)
    {
      std::cout << "iwa2[" << i << "] = " << iwa2[i] << "\n";
    }

  for(int i = 1 ;i < 6 ; i++)
    {
      std::cout << "iwa3[" << i << "] = " << iwa3[i] << "\n";
    }



}


void MatrixTest::testSrtdat()
{




  int n,nz;
  n = 5;
  nz = 10;

  Matrix matrix(n,n,nz,false);
  int *iwa1;

  iwa1 = new int[n+1];
  int indRow[11] = {0,1,1,2,2,3,3,4,4,5,5};
  int indCol[11] = {0,1,2,1,3,2,4,3,4,4,5};

  int jpntr[7];
  int ipntr[7] = { 0,1,3,5,7,9,11};

  int outputRow[11] = { 0,1,2,1,3,2,4,3,4,5,5};
  int outputCol[11] = { 0,1,1,2,2,3,3,4,4,4,5};
  int outputJpntr[7] = {0,1,3,5,7,10,11};
  int outputIpntr[7] = { 0,1,3,5,7,9,11};


  CPPUNIT_FAIL("Test Case need to be reimplemented");
  // matrix.srtdat();

  for(int i = 1 ; i <= 10; i++)
    {
      CPPUNIT_ASSERT_EQUAL(outputRow[i],indRow[i]);
      CPPUNIT_ASSERT_EQUAL(outputCol[i],indCol[i]);

    }

  for (int i = 1; i <= 6 ; i++)
    {
      CPPUNIT_ASSERT_EQUAL(outputJpntr[i],jpntr[i]);
    }

  // Second Test.
  // matrix.srtdat();


  int secondOutputRow[11] = { 0,1,1,2,2,3,3,4,4,5,5};
  int secondOutputCol[11] = { 0,1,2,1,3,2,4,3,4,4,5};

  for(int i = 1 ; i <= 10; i++)
    {
      CPPUNIT_ASSERT_EQUAL(secondOutputRow[i],indRow[i]);
      CPPUNIT_ASSERT_EQUAL(secondOutputCol[i],indCol[i]);

    }

  for (int i = 1; i <= 6 ; i++)
    {
      CPPUNIT_ASSERT_EQUAL(outputIpntr[i],ipntr[i]);
    }


}

/**
 * This is the third test for srtdat
 *
 */

void MatrixTest::testSrtdat_3()
{
  /**
   ****************************************
   * Third Test
   ****************************************
   **/
  /**
   * Here we have done only a static test case, but we should start
   * developing a method so that we can test our methods by just
   * replacing files.
   */

  /**
   * Suppose we have two file, one is srtdat.in and the other is
   * srtdat.out. Then we can test those.
   */

  /**
   * Let's make the srtdat.in file.
   **/

  /**
   * Let's read some data from that.
   **/
  FILE* f;
  int m,n,nz,*indRow,*indCol;
  double *val;

  int ret_code;
  MM_typecode matcode;

  if((f = fopen("srtdat.in","r")) == NULL)
    CPPUNIT_FAIL("strtdat.in file cannot be opened");

  if(mm_read_banner(f,&matcode) != 0)
    CPPUNIT_FAIL("Could not process Matrix market baner.\n");

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
      mm_is_sparse(matcode))
    {
      CPPUNIT_FAIL("We cannot handle MATRIP of this type");
    }

  if((ret_code = mm_read_mtx_crd_size(f,&m,&n,&nz) != 0) )
    CPPUNIT_FAIL("Cannot Read the size information from the file");




//   for (int i = 1; i <= nz; i++)
//     {
//       fscanf(f, "%d %d %lg]n", indRow+i, indCol+i, val+i);
//     }

  CPPUNIT_FAIL("Test Case Not Implemented");
}



/**
 * This is another test Method.
 *
 */

void MatrixTest::testSetr()
{


  int n,nz;
  int m = 5;
  n = 5;
  nz = 10;

  Matrix matrix(m,n,nz,false);
  int *iwa1;

  iwa1 = new int[n+1];
  int indRow[11] = {0,1,1,2,2,3,3,4,4,5,5};
  int indCol[11] = {0,1,2,1,3,2,4,3,4,4,5};
  int jpntr[7];
  int ipntr[7] = { 0,1,3,5,7,9,11};

  int outputRow[11] = { 0,1,2,1,3,2,4,3,4,5,5};
  int outputCol[11] = { 0,1,2,1,3,2,4,3,4,4,5};
  int outputJpntr[7] = {0,1,3,5,7,10,11};
  int outputIpntr[7] = { 0,1,3,5,7,9,11};

  // matrix.srtdat();

  // matrix.setr();
  CPPUNIT_FAIL("Test Case need to be reimplemented");

  for(int i = 1 ; i <= 10; i++)
    {
      // std::cout << outputRow[i] << " " << indRow[i] << "\n";
      CPPUNIT_ASSERT_EQUAL(outputRow[i],indRow[i]);
      CPPUNIT_ASSERT_EQUAL(outputCol[i],indCol[i]);

    }

  for (int i = 1; i <= 6 ; i++)
    {
      CPPUNIT_ASSERT_EQUAL(outputJpntr[i],jpntr[i]);
    }

  for (int i = 1; i <= 6 ; i++)
    {
      CPPUNIT_ASSERT_EQUAL(outputIpntr[i],ipntr[i]);
    }
}
void MatrixTest::testSetr_1()
{


  int n,nz;
  int m = 5;
  n = 5;
  nz = 10;

  Matrix matrix(m,n,nz,false);
  int *iwa1;

  iwa1 = new int[n+1];
  int indRow[11] = {0,1,1,2,2,3,3,4,4,5,5};
  int indCol[11] = {0,1,5,1,3,2,4,3,4,4,5};
  int jpntr[7];
  int ipntr[7] = { 0,1,3,5,7,9,11};

  int outputRow[11] = { 0,1,2,3,4,2,5,3,4,1,5};
  int outputCol[11] = { 0,1,5,1,3,2,4,3,4,4,5};
  int outputJpntr[7] = {0,1,3,4,6,9,11};
  int outputIpntr[7] = { 0,1,3,5,7,9,11};

  // matrix.srtdat();

  // matrix.setr();
  CPPUNIT_FAIL("Test Case need to be reimplemented");

  for(int i = 1 ; i <= 10; i++)
    {
      // std::cout << outputRow[i] << " " << indRow[i] << "\n";
      CPPUNIT_ASSERT_EQUAL(outputRow[i],indRow[i]);
      CPPUNIT_ASSERT_EQUAL(outputCol[i],indCol[i]);

    }
  for (int i = 1; i <= 6 ; i++)
    {
      CPPUNIT_ASSERT_EQUAL(outputJpntr[i],jpntr[i]);
    }

  for (int i = 1; i <= 6 ; i++)
    {
      CPPUNIT_ASSERT_EQUAL(outputIpntr[i],ipntr[i]);
    }
}


/**
 * This is another test Method.
 *
 */

void MatrixTest::testSeq()
{
  CPPUNIT_FAIL("Test Case not Implemented!");
}


/**
 * This is another test Method.
 *
 */

void MatrixTest::testDsm()
{
  CPPUNIT_FAIL("Test Case not Implemented!");
}


void MatrixTest::testDegr()
{
  CPPUNIT_FAIL("Test Case not Implemented!");
}


void MatrixTest::test_getColMax()
{
  CPPUNIT_ASSERT_EQUAL(26,matrix->getColMax());
}


void MatrixTest::testGetColMin()
{
  CPPUNIT_ASSERT_EQUAL(2,matrix->getColMin());
}

void MatrixTest::testGetColAvg()
{
  CPPUNIT_ASSERT_EQUAL(8,matrix->getColAvg());
}

void MatrixTest::testGetRowMax()
{
  CPPUNIT_ASSERT_EQUAL(6,matrix->getRowMax());
}

void MatrixTest::testGetRowMin()
{
  CPPUNIT_ASSERT_EQUAL(1,matrix->getRowMin());
}

void MatrixTest::testGetRowAvg()
{
  CPPUNIT_ASSERT_EQUAL(4,matrix->getRowAvg());
}


/* void GetCurrentPath( char* buffer)
{
  getcwd( buffer, PATH_MAX);
  }
*/

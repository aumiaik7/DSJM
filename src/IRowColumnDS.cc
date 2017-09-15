#include "IRowColumnDS.hh"
#include <iostream>
#include <cmath>
// (setq fill-column 100)


using namespace std;

/*
 * Purpose:         Constructor for IRowColumnDS. It allocates the necessary arrays for storing the
 *                  sparse matrix in CCS and CRS format.
 *
 * Pre-condition: 	None.
 *
 * Post-condition:  IRowColumnDS object is constructed with memory allocated to <id:row_ind>,
 *                  <id:col_ind>, <id:jpntr>, <id:ipntr> and <id:x>. <id:x> is is only allocated
 *                  memory if member variable <id:value> evaluates to true.
 */

IRowColumnDS::IRowColumnDS(int M, int N, int p_nz, bool p_value)
    :  M(M),
       N(N),
       nz(p_nz),
       value(p_value),
       entry_index(1)
{
    row_ind = new int[nz + 1];
    col_ind = new int[nz + 1];
    jpntr = new int[N + 2];
    ipntr = new int[M + 2];
    if(value)
        x = new double[nz+1];
}
/* IRowColumnDS() ENDS*/

/**
 * @Description : Destructor
 * @Creation Date: 17.08.2009
 */
IRowColumnDS::~IRowColumnDS()
{
    delete[] row_ind;
    delete[] col_ind;
    delete[] jpntr;
    delete[] ipntr;
    if(value)
        delete[] x;
}
/* ~IRowColumnDS() ENDS*/

/*
 * Purpose:         Computes Compressed Column Storage(CCS) format of the sparse matrix. The CCS
 *          		format stores the columns of matrix A in three member arrays in IRowColumnDS
 *          		object: <id:jpntr>, <id:row_ind> and <id:x>. Data member <id:x> is empty if
 *          		<id:value>, a boolean member variable evaluates to false.
 *
 * Pre-condition: 	Assumes that the matrix definition is stored in co-ordinate format in
 *          		<id:row_ind> and <id:col_ind> integer array. For every non-zero position in the
 *          		sparse matrix there is two entry: row_ind[i] and col_ind[i] holding the row, and
 *          		column coordinate of the nonzero entry. If we <id:value> is true then x[i]
 *          		stores the corresponding nonzero item.
 *
 * Post-condition: 	Column-oriented definition of the sparse matrix is stored in the two array
 *          		<id:jpntr> and <id:row_ind>. If value of the nonzero items are being stored ,
 *          		then <id:x> is also organized in column oriented definition.
 *
 * Return values:   Returns true when the function is executed successfully, otherwise returns
 *          		false.
 *
 */
bool IRowColumnDS::computeCCS()
{
    int *tag;
    try
    {
        tag = new int[N+1]; //  Temporary working array of size N+1.

        // Initialize the working array <id:w>
        for(int j=1; j <=N ; j++)
        {
            tag[j] = 0;
        }

        // Store the number of nonzeroe entries in w[i] for each
        // column i
        //std::cout << "printing the value of nz: " << std::endl;
        //std::cout << nz << std::endl;
        for(int k = 1; k <= nz; k++)
        {
            tag[col_ind[k]] = tag[col_ind[k]] + 1;
        }

        // Compute jpntr such that jpntr[j] points to the beginning of
        // row indices for column j in <id:row_ind> integer array.
        jpntr[1] = 1;
        for(int j = 1 ; j <= N; j++)
        {
            jpntr[j+1] = jpntr[j] + tag[j];
            tag[j] = jpntr[j];

        }

        // Perform In Place Sort.
        int k = 1;
        do
        {
            int j = col_ind[k];
            if (k >= jpntr[j])
            {
                // Current element is already in right position. Examine the next element of the
                // first un-sorted element in the j-th group.
                k = max(k+1,tag[j]);
            }
            else
            {
                // Current element is not in right position, place the element in position, and set
                // the displaced entry as the current element.
                int l = tag[j]++;
                int i = row_ind[k];
                row_ind[k] = row_ind[l];
                col_ind[k] = col_ind[l];
                row_ind[l] = i;
                col_ind[l] = j;
                if(value) // If we are storing values in the IRowColumnDS.
                {
                    double t = x[k];
                    x[k] = x[l];
                    x[l] = t;
                }
            }

        }while(k <=nz);
    }
    catch(std::bad_alloc)
    {
        std::cerr << "Memory Exhausted " << std::endl;
        delete[] tag;
        return false;
    }
    delete[] tag;
    return true;
}

/**
 * Method max
 */
int IRowColumnDS::max(int s, int t)
{
    if(s >= t)
        return s;
    else
        return t;
}

/*
 * Purpose:         Removes duplicate entries from the column-oriented definition of the sparsity,
 *                  and compresses three class member <id:row_ind>, <id:jpntr> and <id:x> array.
 *
 * Pre-condition: 	Assumes that the sparsity pattern is stored in column-oriented definition in
 *          		<id:jpntr>, <id:row_ind> and <id:x> array by	calling computeCCS() method.
 *
 * Post-condition: 	Removes duplicate entry and reorganize <id:row_ind>, <id:jpntr> and <id:x>
 *                  array.
 *
 * Return values:   Returns number of unique nonzero items when the function is executed
 *          		successfully, otherwise returns zero.
 *
 */

// TODO: Shouldn't we also change the value array x[] here ?
int IRowColumnDS::compress()
{

    int *tag;

    try
    {
        tag = new int[M+1]; // Tag is a working array of size M+1.


        // Initialization of <id:tag> array.
        for (int i = 1; i <= M; i++)
        {
            tag[i] = 0;
        }

        int nnz = 1;
        int k;

        for (int j = 1; j <= N; ++j)
        {
            k = nnz;
            for (int jp = jpntr[j]; jp < jpntr[j + 1] ; ++jp)
            {
                int ir = row_ind[jp];
                if (tag[ir] != j)
                {
                    row_ind[nnz] = ir;
                    tag[ir] = j;
                    if(value)
                    {
                        x[nnz] = x[jp]; // TODO: Check for bugs here. Added on September 10, 2010
                    }
                    ++nnz;
                }
            }
            jpntr[j] = k;
        }
        jpntr[N + 1] = nnz;

        // Clean up.
        if(tag) delete[] tag;
        this->nnz = nnz -1 ;
        return nnz -1;
    }
    catch(std::bad_alloc)
    {
        if(tag) delete[] tag;
        return 0;
    }
}

/*
 * Purpose:         Computes Compressed Row Storage(CCS) definition of the sparse matrix. The CRS
 *          		format stores the rows of matrix A in two member arrays in IRowColumnDS
 *          		object: <id:ipntr>, <id:col_ind>. Value array <id:x> is not stored in row
 *          		oriented definition.
 *
 * Pre-condition: 	Assumes that the matrix definition is stored in CCS format in
 *          		<id:row_ind> and <id:jpntr> integer array and duplicate entries has been
 *          		removed by calling computeCCS() method and compress() method.

 *
 * Post-condition: 	Row-oriented definition of the sparse matrix is stored in the two array
 *          		<id:ipntr> and <id:col_ind>.
 *
 * Return values:   Returns true when the function is executed successfully, otherwise returns
 *          		false.
 *
 */
bool IRowColumnDS::computeCRS()
{
    /**     **********
     *
     *     subroutine setr
     *
     *     given a column-oriented definition of the sparsity pattern
     *     of an M by N matrix a, this subroutine determines a
     *     row-oriented definition of the sparsity pattern of a.
     *
     *     on input the column-oriented definition is specified by
     *     the arrays row_ind and jpntr. on output the row-oriented
     *     definition is specified by the arrays col_ind and ipntr.
     *
     *     the subroutine statement is
     *
     *       subroutine setr(M,N,row_ind,jpntr,col_ind,ipntr,tag)
     *
     *     where
     *
     *       M is a positive integer input variable set to the number
     *         of rows of a.
     *
     *       N is a positive integer input variable set to the number
     *         of columns of a.
     *
     *       row_ind is an integer input array which contains the row
     *         indices for the non-zeroes in the matrix a.
     *
     *       jpntr is an integer input array of length N + 1 which
     *         specifies the locations of the row indices in row_ind.
     *         the row indices for column j are
     *
     *               row_ind[k], k = jpntr[j],...,jpntr[j+1]-1.
     *
     *         note that jpntr[N+1]-1 is then the number of non-zero
     *         elements of the matrix a.
     *
     *       col_ind is an integer output array which contains the
     *         column indices for the non-zeroes in the matrix a.
     *
     *       ipntr is an integer output array of length M + 1 which
     *         specifies the locations of the column indices in col_ind.
     *         the column indices for row i are
     *
     *               col_ind[k], k = ipntr[i],...,ipntr[i+1]-1.
     *
     *         note that ipntr[1] is set to 1 and that ipntr[M+1]-1 is
     *         then the number of non-zero elements of the matrix a.
     *
     *       w is an integer work array of length M.
     *
     *     argonne national laboratory. minpack project. july 1983.
     *     thomas f. coleman, burton s. garbow, jorge j. more'
     *
     **/
    int *tag;
    try
    {
        tag = new int[M+1]; // Temporary working array of size N+1.
        int ir,jcol,jp;

        // Initialize the working array <id:w>
        for (int ir = 1; ir <=  M ;ir++  )
        {
            tag[ir] = 0;
        }
        // Store the number of nozero entries in w[i] for each row i.
        for (int jp = 1; jp <=  jpntr[N+1]-1 ;jp++  )
        {
            tag[row_ind[jp]] = tag[row_ind[jp]] + 1;
        }

        // Compute ipntr such that ipntr[ir] points to the beginning of column indices for row ir in
        // <id:row_ind> integer array.
        ipntr[1] = 1;
        for (ir = 1; ir <=  M ;ir++  )
        {
            ipntr[ir+1] = ipntr[ir] + tag[ir];
            tag[ir] = ipntr[ir];
        }

        // Fill col_ind traversing the row-oriented definition of the sparse Matrix A.
        for (jcol = 1; jcol <=  N ;jcol++  )
        {
            for (jp = jpntr[jcol]; jp <=  jpntr[jcol+1]-1 ;jp++  )
            {
                ir = row_ind[jp];
                col_ind[tag[ir]] = jcol;
                tag[ir] = tag[ir] + 1;
            }
        }

    }
    catch(std::bad_alloc)
    {
        if(tag) delete[] tag;
        return false;
    }
    delete[] tag;
    return true;
}


void IRowColumnDS::setIndRowEntry(int index,int entry)
{
    *(row_ind+index) = entry;
}
void IRowColumnDS::setIndColEntry(int index,int entry)
{
    *(col_ind+index) = entry;
}

int IRowColumnDS::getIndRowEntry(int index) const
{
    return *(row_ind + index);
}
int IRowColumnDS::getIndColEntry(int index) const
{
    return *(col_ind + index);
}

int IRowColumnDS::getJpntrEntry(int index) const
{
    return *(jpntr + index);
}

int IRowColumnDS::getIpntrEntry(int index) const
{
    return *(ipntr + index);
}

/*
 * Purpose:         To access the value of nonzero items.
 *
 * Pre-condition:   Values have been stored through entry() method of IRowColumnDS.
 *
 * Return values:   Value stored in x[i] location.
 *
 */
double IRowColumnDS::getX(int index) const
{
    return *(x+index);
}



int IRowColumnDS::getRowMin() const
{
    int  min = MAX_INT;
    for (int i = 2; i <= M+1; i++)
    {
        int tmp = ipntr[i] - ipntr[i-1];
        if(tmp < min)
            min = tmp;
    }
    return min;

}

int IRowColumnDS::getRowMax() const
{
    int max = 0;
    for (int i = 2; i <= M+1; i++)
    {
        int tmp = ipntr[i] - ipntr[i-1];
        if(tmp > max)
        {
            max = tmp;
            std::cout << max << "at row " << (i-1) << std::endl;
        }
    }
    return max;
}

int IRowColumnDS::getRowAvg() const
{
    return round(nnz/M);
}


int IRowColumnDS::getColMax() const
{
    int max = 0;
    for (int i = 2; i <= N+1; i++)
    {
        int tmp = jpntr[i] - jpntr[i-1];
        if(tmp > max)
        {
            max = tmp;
        }
    }
    return max;
}

int IRowColumnDS::getColMin() const
{
    int  min = MAX_INT;
    for (int i = 2; i <= N+1; i++)
    {
        int tmp = jpntr[i] - jpntr[i-1];
        if(tmp < min)
            min = tmp;
    }
    return min;
}

int IRowColumnDS::getColAvg() const
{
    return round(nnz/N);
}

int IRowColumnDS::getNNZ() const
{
    return nnz;
}

/**
 * Method min(int s, int t)
 * */
int IRowColumnDS::min(int s, int t)
{
    if(s <= t)
        return s;
    else
        return t;
}


/*
 * Purpose:         Make an entry for each nonozero item in the Sparse Matrix A.
 *
 * Pre-condition: 	None.
 *
 * Post-condition: 	An entry for nonzero location for a sparse matrix is inserted in the co-ordinate
 *          		storage defined by <id:row_ind>, <id:col_ind>, and <id:x>. For i-th invocation of
 *          		entry() method, the row and column for the nonzero item is stored in indRow[i],
 *          		and col_ind[i], where x[i] holds the nonzero item. If member variable <id:value>
 *          		is false, then the nozero items literal value is not inserted in <id:x> array.
 *
 * Return values:   void
 *
 */
void IRowColumnDS::entry(int row, int col, double value)
{
    setIndRowEntry(entry_index,row);
    setIndColEntry(entry_index,col);
    if(this->value)
        x[entry_index] = value;
    entry_index++;
}

/*
 * Purpose:         Overloaded version of method entry(row,col,value). This function is used when
 *                  the matrix object A is storing only the sparsity pattern excluding the nonzero
 *                  item itself.
 * Pre-condition: 	None.
 *
 * Post-condition: 	An entry for nonzero location for a sparse matrix is inserted in the co-ordinate
 *          		storage defined by <id:row_ind>, <id:col_ind>, without inserting value in for
 *          		<id:x> array.
 *
 * Return values:   void
 *
 */
void IRowColumnDS::entry(int row, int col)
{
    entry(row,col,0);
}

int IRowColumnDS::getM() const
{
    return M;
}
int IRowColumnDS::getN() const
{
    return N;
}

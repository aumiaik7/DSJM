#ifndef IROWCOLUMNDS_HH
#define IROWCOLUMNDS_HH

#define MAX_INT 32767

class IRowColumnDS
{
public:
    IRowColumnDS(int M,int N, int nz, bool value);
    virtual ~IRowColumnDS();

    // ========================================
    // RowColumnDS Begins
    // ========================================

    bool computeCCS();
    /**
     * Method max
     */

    int max(int s, int t);

    int compress();
    bool computeCRS();


    /**
     * Setter getter methods
     */
    void setIndRowEntry(int index,int entry);
    void setIndColEntry(int index,int entry);

    void entry(int row,int col, double value);
    void entry(int row, int col);

    int getIndRowEntry(int index) const;
    int getIndColEntry(int index) const;

    int getJpntrEntry(int index) const;
    int getIpntrEntry(int index) const;

    double getX(int index) const;



    /**
     * Statistical Methods
     */
    int getRowMax() const;
    int getRowMin() const;
    int getRowAvg() const;
    int getColMax() const;
    int getColMin() const;
    int getColAvg() const;

    int getNNZ() const;

    int getM() const;
    int getN() const;

    /**
     * Method min(int s, int t)
     * */
    int min(int s, int t);

    // RowColumnDS Ends
    // ========================================



protected:
    int M; // Row dimension
    int N; // Column Dimension
    int nz; // Reported Number of nonzeroes in the file
    int nnz; // Number of nonzeroes after the Compression.
    // In most of the cases you are going to use this one,
    // This one is going to be used after the compression
    // or the deletion of the duplicate values happen.

    int *row_ind;
    int *col_ind;
    int *jpntr;
    int *ipntr;
    bool value;

    double *x;

    int entry_index;
};

#endif

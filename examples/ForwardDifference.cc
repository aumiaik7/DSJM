#include "StructuredJacobian.hh"
#include "FD.hh"
#include <iostream>
#include <cstdlib>
#include "Eye.hh"
#include <cstring>
#include <fstream>
#include "Timer.h"
#include "CLI.h"
#include <getopt.h>

// It sets p[i] equal to the sum of c[1] through c[i-1]+1
// On output c[1..n] is overwritten with p[1..n]
void cumsum(int *p, int *c, int n)
{
    int i,nz = 1;

    if( !p || !c) return; /* Check inputs */

    for (i = 1; i <= n; i++)
    {
        p[i] = nz;
        nz += c[i];
        c[i] = p[i];
    }
    p[i] = nz;
}

Matrix* reconstruct(const Matrix* B, Matrix* matrix, int *ngrp)
{
    int maxgrp = matrix->getNumberOfColors();
    int N = matrix->getN();
    int M = matrix->getM();
    int NZ = matrix->getNNZ();
    // w:
    //    Working memroy of maxgrp size.
    //    The indexing is 1 based.
    int *w = new int[maxgrp+1];

    // Initialize the working memory w to 0
    for(int color = 1; color <= maxgrp; color++) w[color] = 0;

    // Count the number of columns for each color in w.
    for(int jcol = 1; jcol <= N; jcol++) w[ngrp[jcol]]++;

    // cpntr:
    //    Base-1 Index.
    //    Maps a color(Base-1) to its member columns. The indices of
    //    the columns can be found from cpntr[color] and
    //    cpntr[color+1].
    int *cpntr = new int[maxgrp + 2];

    // cpntr holds the indices in the colindex array to find the
    // member columns for a color.
    cumsum(cpntr,w,maxgrp);

    // colindex:
    //         Columns belonging to a color can be found using
    //         colindex[i] where cpntr[color] <= i < cpntr[color+1]
    //         Indexing is done from 1 to N. And the column indices
    //         found are also in Base-1 Format.
    int *colindex = new int[N+1];

    // We build the Colindex Array. Copying the column index to the
    // location pointed by w and incrementing the pointer to the next
    // postion.
    for(int jcol = 1; jcol <= N ; jcol++) colindex[w[ngrp[jcol]]++] = jcol;

    // inverseOfRow:
    //     In a color group, one row can belong to only one
    //     column. This array will record the column number for each
    //     row in that color group.
    int *inverseOfRow = new int[M+1];

    // So we start to reconstruct the column.
    Matrix *mT;
    mT = new Matrix(M,N,NZ,true);

    // Reconstruction takes place column by column of the Matrix B. We
    // distribute the nonzeros of the B matrix in the appropriate
    // columns.
    for(int color = 1; color <= maxgrp; color++)
    {
        // Step 1: Build the InverseofRow array first.
        // Find all the columns belonging to this color group.
        // Then for each column find all the rows and record the
        // corresponding column entry in inverseOfRow Array.
        // cout << "Color : " << color << endl;
        for(int cp = cpntr[color]; cp < cpntr[color+1]; cp++)
        {
            int jcol = colindex[cp];
            for(int jp = matrix->getJpntrEntry(jcol) ; jp < matrix->getJpntrEntry(jcol + 1); jp++)
            {
                int ir = matrix->getIndRowEntry(jp);
                inverseOfRow[ir] = jcol;
            }
        }

        // Step 2: Redistribute the nonzerof of a color to the
        // corresponding columns.
        for(int p = B->getJpntrEntry(color) ; p < B->getJpntrEntry(color+1) ; p++)
        {
            int ir = B->getIndRowEntry(p);
            double x = B->getX(p);
            int col = inverseOfRow[ir];
            mT->entry(ir,col,x);
        }
    }
    mT->computeCCS();
    mT->compress();
    mT->computeCRS();
    delete[] w;
    delete[] cpntr;
    delete[] inverseOfRow;
    delete[] colindex;
    return mT;
}

void display_usage(int argc,char *argv[])
{
    fprintf (stderr, "Usage: %s -n number -m dd/fd/cfd", argv[0]);
    exit(1);
}

int main(int argc, char *argv[])
{

    if(argc < 4)
    {
        fprintf (stderr, "Usage: %s -n number -m dd/fd/cfd", argv[0]);
        exit(1);
    }

    char optstring[] = "m:o:n:x:d:h?v";
    int opt = getopt(argc,argv, optstring);
    int n = -1;
    double d;
    char *inputfile = NULL;
    CLI::ordering_method ord = CLI::SLO;
    CLI::approximation_method method = CLI::CFD;
    bool isVerify = false;
    while(opt != -1)
    {
        switch(opt)
        {
        case 'm':
            method = CLI::parseApproximationMethod(optarg);
            break;
        case 'o':
            ord = CLI::parseMethod(optarg);
            break;
        case 'n':
            n = atoi(optarg);
            break;
        case 'x':
            if(optarg != NULL && strlen(optarg ) > 0)
            {
                if(inputfile != NULL)
                    delete[] inputfile ;
                inputfile = new char[strlen(optarg) + 1];
                strcpy(inputfile, optarg);
            }
            break;
        case 'd':
            d = atof(optarg);
            break;
        case 'v':
            isVerify = true;
            break;
        case 'h':
        case '?':
            display_usage(argc, argv);
        default:
            break;
        }
        opt = getopt(argc,argv,optstring);
    }

    if (n == -1)
        exit(1);

    // Build the Structured Jacobian Object with the $n$ parameter.
    ColPack::Timer timer;
    double time;

    timer.Start();
    StructuredJacobian structuredJacobian(n);
    Matrix *mT = NULL;

    // A Vector(point) in the R^n Space.
    vector<double> x;
    x.push_back(0);
    for(int i = 0; i < n ; i++) x.push_back(1.0);

    // ifstream f(inputfile);
    // delete[] inputfile;
    // if(!f.is_open())
    // {
    //     cerr << "Input File cannot be open " << endl;
    // }
    // for(int i = 0; i <= n; i++)
    // {
    //     double val;
    //     f >> val;
    //     x.push_back(val);
    // }
    //f.close();

    // identity matrix.

    if(method == CLI::FD)
    {
        Eye eye(n);
        mT = FD::finiteDifference(structuredJacobian,x,&eye);

    }
    else if(method == CLI::DD)
    {
        mT = structuredJacobian.derivative(x);
    }
    else if(method == CLI::CFD)
    {
        auto_ptr<Sparsity> sparsity = structuredJacobian.getSparsityPattern();
        Matrix matrix(n,n,sparsity->size()-1,false);
        for(int i = 1; i <= sparsity->size()-1; i++)
        {
            matrix.setIndRowEntry(i,(*sparsity)[i].first);
            matrix.setIndColEntry(i,(*sparsity)[i].second);
        }

        matrix.computeCCS();
        int nnz = matrix.compress();
        matrix.computeCRS();

        matrix.computedegree();
        int *list = NULL;
        int *ngrp = NULL;
        Matrix *seed = NULL;
        Matrix *m = NULL;
        try
        {
            list = new int[n+1];
            ngrp = new int[n+1];
            bool success = matrix.slo(list); // TODO: Change to the method supplied from
                                             // outside.
            int maxgrp = matrix.greedycolor(list, ngrp);

            // cout << "Coloring : " << maxgrp << endl;

            seed = matrix.getSeedMatrix(ngrp);

            m = FD::finiteDifference(structuredJacobian,x,seed);
            // cs_print(T,0);

            mT = reconstruct(m,&matrix, ngrp);
            // cs_print(A,0);
            delete[] list;
            delete[] ngrp;
            delete seed;
            delete m;
        }
        catch(std::bad_alloc)
        {
            std::cerr << "Memory Exhausted" << std::endl;
            if(list) delete[] list;
            if(ngrp) delete[] ngrp;
            if(seed) delete seed;
            if(m) delete m;
            return 1;
        }



    }
    timer.Stop();

    time = timer.GetWallTime();


    cout << "-" << endl;
    cout << "    n: " << n << endl;
    cout << "    method: " << CLI::approximation_method_to_string(method) << endl;
    cout << "    time: " << time << endl;

    delete[] inputfile;
    if(isVerify)
    {
        string s("output");
        s += CLI::approximation_method_to_string(method);
        s += ".txt" ;
        ofstream f(s.c_str());

        for(int j = 1; j <= n; j++)
        {
            for( int jp = mT->getJpntrEntry(j); jp < mT->getJpntrEntry(j+1) ; jp++)
            {
                int ir = mT->getIndRowEntry(jp);
                f << ir-1 << " " << j-1 << " " << mT->getX(jp) << endl;
            }
        }

        f.close();
    }
    if(mT)
        delete mT;
    return 0;
}

#include "StructuredJacobian.hh"
#include <iostream>
#include <utility>
#include <algorithm>
#include "Matrix.hh"

StructuredJacobian::StructuredJacobian(int n)
    :n(n),m(n)
{
}

StructuredJacobian::~StructuredJacobian()
{
}


/*
 * This function returns the Sparsity pattern
 * TODO: We have to write this sparsity pattern to a file.
 */
auto_ptr<Sparsity> StructuredJacobian::getSparsityPattern()
{
    auto_ptr<Sparsity> sparsity (new Sparsity());
    sparsity->push_back(std::make_pair(0,0));

    for (int i = 1; i <= n; i++)
    {
        if (i > 1 && (i-1) < (n-4)) sparsity->push_back(std::make_pair(i,i-1));
        if(i < (n-4))    sparsity->push_back(std::make_pair(i,i));
        if((i+1) < (n-4))    if (i < n) sparsity->push_back(std::make_pair(i,i+1));

        if (n > 4) sparsity->push_back(std::make_pair(i,n-4));
        if (n > 3) sparsity->push_back(std::make_pair(i,n-3));
        if (n > 2) sparsity->push_back(std::make_pair(i,n-2));
        if (n > 1) sparsity->push_back(std::make_pair(i,n-1));
        sparsity->push_back(std::make_pair(i,n));
    }

    return sparsity;
}

/*
 * Function Evaluation
 */
auto_ptr<vector<double> > StructuredJacobian::eval(const vector<double>& x)
{
    auto_ptr<vector<double> > y (new vector<double>());
    y->push_back(0);
    y->push_back(-2 * x[1] * x[1] + 3 * x[1] - 2 * x[2] + 3 * x[n-4] - x[n-3] - x[n-2] + 0.5 * x[n-1] - x[n] + 1);
    for(int i = 2; i < n; i++)
    {
        y->push_back(-2 * x[i] * x[i] + 3 * x[i] - x[i-1] - 2 * x[i+1] + 3 * x[n-4] - x[n-3] - x[n-2] + 0.5 * x[n-1] - x[n] + 1);
    }
    y->push_back(-2 * x[n] * x[n] + 3 * x[n] - x[n-1] + 3 * x[n-4] - x[n-3] - x[n-2] + 0.5 * x[n-1] - x[n] + 1);

    return y;
}

/*
 * Returns the Derivative or the Jacobian Matrix given a point
 * TODO: This is not wrapped into an auto_ptr ???
 */
double* StructuredJacobian::arrayderivative(const vector<double>& x)
{
    double *result = new double[n*n + 1];
    for(int i = 1; i <= n*n; i++)
    {
        *(result+i)  = 0;
    }

    for(int i =1; i <= n; i++)
    {
        *(result+ (i-1) * n + i) += -4 * x[i] + 3;

        if(i > 1)
            *(result+ (i-1) * n + i-1) += -1;

        if(i < n)
            *(result+ (i-1) * n + i+1) += -2;

        *(result+ (i-1) * n + n-4) += 3;

        *(result+ (i-1) * n + n-3) += -1;

        *(result+ (i-1) * n + n-2) += -1;

        *(result+ (i-1) * n + n-1) += 0.5;

        *(result+ (i-1) * n + n) += -1;
    }

    return result;
}

/*
 * Returns the Derivative or the Jacobian Matrix given a point
 * TODO: This is not wrapped into an auto_ptr ???
 */
Matrix* StructuredJacobian::derivative(const vector<double>& x)
{
    Matrix *matrix = new Matrix(n,n,getNzmax(n), true);

    pair<int,double> w[8];
    for(int i =1; i <= n; i++)
    {
        int size = 0;

        w[size++] = make_pair( i,-4 * x[i] + 3);

        if(i > 1)
        {
            w[size++] = make_pair(i-1, -1);
        }
        if(i < n)
        {
            w[size++] = make_pair(i+1, -2);
        }

        w[size++] = make_pair(n-4,3);
        w[size++] = make_pair(n-3,-1);
        w[size++] = make_pair(n-2,-1);
        w[size++] = make_pair(n-1,0.5);
        w[size++] = make_pair(n,-1);
        sort(w, w+size);
        for (int j = 1; j < size; j++)
        {
            if(w[j].first == w[j-1].first)
            {
                w[j].second = w[j].second + w[j-1].second;
                w[j-1].first = -1;
            }
        }

        for(int j = 0; j < size ; j++)
        {
            if(w[j].first != -1)
            {
                matrix->entry(i,w[j].first,w[j].second);
            }
        }
    }

    matrix->computeCCS();
    int nnz = matrix->compress();
    matrix->computeCRS();

    return matrix;
}

int StructuredJacobian::getNzmax(const int n) const
{
    if(n < 6) return n * n;
    else return 8 * n - 16;
}



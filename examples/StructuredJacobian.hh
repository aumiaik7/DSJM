#ifndef STRUCTUREDJACOBIAN_HH
#define STRUCTUREDJACOBIAN_HH

#include <memory>
#include <vector>
#include <utility>
#include "Matrix.hh"

using namespace std;

/*
 * Structured Jacobian Problem. This function recieves an input in
 * vector, where the indices start from 1. and return another vector.
 */

// TODO: Should We return the vector as auto_ptr or simple array ?

typedef vector<pair<int,int> > Sparsity;

class StructuredJacobian
{
public:
    StructuredJacobian(int n);
    virtual ~StructuredJacobian();

    auto_ptr<Sparsity> getSparsityPattern();

    auto_ptr<vector<double> > eval(const vector<double>& x);

    Matrix* derivative(const vector<double>& x);
    double* arrayderivative(const vector<double>& x);

    int getNzmax(const int n) const;

    int n; // Matrix Dimension.
    int m;



private:

};

#endif

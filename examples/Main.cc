#include "StructuredJacobian.hh"
#include "FD.hh"
#include <iostream>
#include <cstdlib>
#include "Eye.hh"
extern "C" {
    // Get Declaration of CSparse
#include "cs.h"
}


void display_usage(int argc,char *argv[])
{
    fprintf (stderr, "Usage: %s -n number\n", argv[0]);
    exit(1);
}

int main(int argc, char *argv[])
{

    if(argc < 2)
    {
        fprintf (stderr, "Usage: %s -n number\n", argv[0]);
        exit(1);
    }

    int opt = getopt(argc,argv, "n:h?");
    int n = -1;
    while(opt != -1)
    {
        switch(opt)
        {
        case 'n':
            n = atoi(optarg);
            break;
        case 'h':
        case '?':
            display_usage(argc, argv);
        default:
            break;
        }
        opt = getopt(argc,argv, "n:h?");
    }

    if (n == -1)
        exit(1);

    // Build the Structured Jacobian Object with the $n$ parameter.
    StructuredJacobian structuredJacobian(n);

    // Get the Sparsity Pattern from the Object.
    auto_ptr<Sparsity> sparsity = structuredJacobian.getSparsityPattern();

    std::cout << "%%MatrixMarket matrix coordinate pattern general" << std::endl;
    std::cout << n << " " << n << " " << sparsity->size()-1<< std::endl;
    // TODO: The size can be replaced by 8n-16 also for n > 6

    for(int i = 1; i <= sparsity->size()-1 ; i++)
    {
        std::cout <<  (*sparsity)[i].first << " " <<  (*sparsity)[i].second << std::endl;
    }

    // A Vector(point) in the R^n Space.
    vector<double> x;
    for(int i = 0; i <= n; i++)
        x.push_back(1);

    // Calculate y = F(x)
    auto_ptr<vector<double> > y = structuredJacobian.eval(x);

    for(int i = 1; i <= n; i++)
        std::cout << (*y)[i] << std::endl;
    std::cout << endl;

    // identity matrix.
    Eye eye(n);
    // TODO: Matrix Object do not have the capacity to store values
    //yet.
    cs* m = FD::finiteDifference(structuredJacobian,x,&eye);


    std::cout << "From Finite Difference" << std::endl;


    int *Ap = m->p, *Ai = m->i;
    double *Ax = m->x;
    for(int j = 0; j < structuredJacobian.n; j++)
    {
        for(int p = Ap[j]; p < Ap[j+1]; p++)
        {
            std::cout.width(5);
            std::cout << Ax[p] << " " ;
        }
        std::cout << endl;
    }


    cs_free(m);

    cs* T = structuredJacobian.derivative(x);

    std::cout << "From Symbolic Difference" << std::endl;
    // for(int i = 1; i <= structuredJacobian.m; i++)
    // {
    //     for(int j = 1; j <= structuredJacobian.n; j++)
    //     {
    //         std::cout.width(5);
    //         std::cout << *(mm + (i-1) * structuredJacobian.n + j) << " " ;
    //     }
    //     std::cout << endl;
    // }

    cs_free(T);

    return 0;
}

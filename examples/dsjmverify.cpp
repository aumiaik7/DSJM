#include <fstream>
#include <iostream>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <cmath>

using namespace std;

int main()
{
    char file1[] = "sortedoutputdd";
    char file2[] = "sortedoutputcfd";
    char file3[] = "sortedoutputfd";

    FILE *f1 = fopen(file1, "r");
    FILE *f2 = fopen(file2, "r");
    FILE *f3 = fopen(file3, "r");

    while(!feof(f1))
    {
        int r1,r2,r3,c1,c2,c3;
        double v1,v2,v3;
        fscanf(f1, "%d %d %lf", &r1,&c1,&v1);
        fscanf(f2, "%d %d %lf", &r2,&c2,&v2);
        fscanf(f3, "%d %d %lf", &r3,&c3,&v3);
        // printf ("%d %d %lf\n" , r1,c1, v1);
        // printf ("%d %d %lf\n" , r2,c2, v2);
        // printf ("%d %d %lf\n" , r3,c2, v3);


        if ( r1 != r2 || r1 != r3 || c1 != c2 || c1 != c3)
        {
            printf("The sorted outputs are not same");
            exit(-1);
        }

        if( fabs(v1 -v2) > 0.01 || fabs(v1 - v3) > 0.01)
        {
            printf("ERROR");
            printf ("%d %d %lf\n" , r1,c1, v1);
            printf ("%d %d %lf\n" , r2,c2, v2);
            printf ("%d %d %lf\n" , r3,c2, v3);

        }


    }
    fclose(f1);
    fclose(f2);
    fclose(f3);

}

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <fstream>


int main(int argc,char *argv[])
{

    FILE *gfp, *cfp;
    if (argc < 3)
    {
        fprintf(stderr, "Usage %s [graph-file] [coloring-file]\n",argv[0]);
        exit(1);
    }
    else
    {
        std::cout << "    vCommand: " << argv[0] <<  " " << argv[1]  <<  " " << argv[2]<< std::endl;
        if((gfp = fopen(argv[1],"r")) == NULL)
        {
            fprintf(stderr,"Graph File %s not found\n",argv[1]);
            exit(1);
        }

        if((cfp = fopen(argv[2],"r")) == NULL)
        {
            fprintf(stderr,"Coloring File %s not found\n",argv[2]);
            exit(1);
        }
    }
    /*
      std::cout << "Verfication Code for a given coloring\n";

      std::cout << "Give us two file name\n";
      std::cout << "One file will have the edges as like as dotty\n";
      std::cout << "And the other file will have the coloring\n";

    */

    /**
     * Read the number of columns or the vertices.
     **/


    int N;
    int jcol;

    fscanf(gfp,"%d",&N);

    // std::cout << "there are "<< N << " number of vertices or columns\n";

    int *coloring  = new int[N+1];
    int i,j;
    int totalColors;
    fscanf(cfp,"%d",&totalColors);
    for (int i = 1;i <= N; i++)
    {
        fscanf(cfp,"%d %d",&jcol,coloring+i);
        if( *(coloring+i) > totalColors)
            std::cerr << "We have color number " << *(coloring+i) << " for column " << i << ", which is higher than totalColors: " << totalColors << std::endl;
        // assert(totalColors >= *(coloring + i));
    }

    if(fscanf(gfp,"%d %d",&i,&j) != EOF)
        do
        {
            if(*(coloring + i) == *(coloring+j) )
            {
                std::cerr << "We have same color for connected edges or vertices\n";
                std::cerr << "So the coloring was not correct\n";

                std::cerr << "Coloring for edge " << i << " is " << *(coloring + i) << "\n";
                std::cerr << "Coloring for edge " << j << " is " << *(coloring + j) << "\n";
                exit(1);
            }
            else
            {
                // std::cout << "Coloring verified for " << i << " and "<< j << "\n";
            }
        }while(fscanf(gfp,"%d %d",&i,&j) != EOF);

    return 0;

}

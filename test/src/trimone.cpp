#include <cstdio> 
#include <iostream>
#include <cstdlib>

#define MAXLINE 100

int main(int argc, char* argv[])
{

  FILE *f; 
  char line[MAXLINE+1];
  if ( argc < 2 ) 
    { 
      // fprintf(stderr,"Usage %s [fileName]\n",argv[0]);
      // exit(1); 
      f = stdin;
    }
  else 
    {
      if (( f = fopen(argv[1],"r")) == NULL) 
	exit(1); 
    }

  fgets(line,MAXLINE,f);
  while(!feof(f))
    {
      fgets(line,MAXLINE,f);
      printf("%s",line);
    }

  fclose (f);
  return 0;
}

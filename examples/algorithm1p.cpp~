#include <fstream>
#include <algorithm> 
#include <iostream>
#include <cstdlib>

#include "gcolor.hh"

using namespace std;

// Declare variables:
int M, N, L, X;
double **d;
double **d1;
double **p;
double **r;
double **fh;
double **fy;

/*Calculate all possibilities of flaws using procedure 1*/
void possibilities()
{
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++)
			if (fh[i][j]!=0){
				for (int k=0; k<N; k++){
					for (int l=0;l<3;l++)
						if (j == (d[k][l]-1)){
							for (int k=0; k<N; k++)
								if (d[k][l]!=0 && (d[k][l]-1)!=j){
									int a = (int) d[k][l];
									p[i][a-1]=65;}
						}
				
				}			
			}
	}

	for (int i=0; i<N; i++){
		for (int j=0; j<3; j++)
			if (fy[i][j]!=0){
				for (int k=0; k<N; k++){
					for (int l=0;l<3;l++)
						if (j == (d[k][l]-1)){
							for (int k=0; k<N; k++)
								if (d[k][l]!=0){
									int b = (int) d[k][l];
									p[i][b-1]=65;}
					}
				
				}		
			}
	}	

}
 

int main()
{
// Open the file:
ifstream fin("Matrix.mtx");

// Ignore headers and comments:
while (fin.peek() == '%') fin.ignore(2048, '\n');

// Read defining parameters:
fin >> M >> N >> L;
cout <<M<<" "<<N<<" "<<L<<std::endl;
// Create your matrix:
double* matrix;			     // Creates a pointer to the array
matrix = new double[M*N];	     // Creates the array of M*N size
fill(matrix, matrix + M*N, 0.); // From <algorithm>, zeros all entries.

// Read the data
for (int l = 0; l < L; l++)
{
	int m, n;
	//double data;
	//fin >> m >> n >> data;
	fin >> m >> n ;
	//matrix[(m-1) + (n-1)*M] = data;
	matrix[(m-1) + (n-1)*M] = 1.0;
}

fin.close();

for (int m = 0; m < M; m++)
{
    for(int n = 0; n < N; n++)
        cout << matrix[m + n*M] << ' ';
    cout << endl;
}

ifstream fin1("coloring.txt");

// Ignore headers and comments:
while (fin1.peek() == '%') fin1.ignore(2048, '\n');

// Read defining parameters:
fin1 >> X;
cout <<X<<endl;
// Create your matrix:
double* matrix1;
double* matrix2;			// Creates a pointer to the array
matrix1 = new double[N];	     // Creates the array of N size
matrix2 = new double[N];
fill(matrix1, matrix1 + N, 0.); // From <algorithm>, zeros all entries.
fill(matrix2, matrix2 + N, 0.); // From <algorithm>, zeros all entries.
// Read the data
for (int l = 0; l < N; l++)
{
	int a, b;
	fin1 >> a >> b ;
	matrix1[l]=a;
	matrix2[l]=b;
}

fin.close();

d1 = new double*[N];
for (int i = 0; i < N; ++i) 
  	d1[i] = new double[X];
double i=1;
for (int i = 0; i < N; i++)
{
	int j=matrix2[i]-1;
      d1[i][j]=matrix1[i];
	//i++;
}
cout<<"Vector D1 based on CPR : "<<endl;
for (int m = 0; m < N; m++)
{
    for(int n = 0; n < X; n++)
        cout << d1[m][n] << ' ';
    cout << endl;
}
double S[N][N];
for (int m = 0; m < N; m++){
    for(int n = 0; n < N; n++){
	if (m == n || m == n-1 ||m == n+1)
		S[m][n]=1;
	else
        	S[m][n]=0;
}
}
cout<<"Vector S : "<<endl;
for (int m = 0; m < N; m++)
{
    for(int n = 0; n < N; n++)
        cout << S[m][n] << ' ';
    cout << endl;
}
//double d[N][3];
//for (int m = 0; m < N; m++)
//{
//    for(int n = 0; n < 3; n++)
//        d[m][n]=0;
//}
d = new double*[N];
for (int i = 0; i < N; ++i) 
  	d[i] = new double[3];
 i=1;
int n=0;
for (int m = 0; m < N; m++)
{
        d[m][n]=i;
	i++;
	n++;
	if (n>2)
		n=0;
}

cout<<"Vector D based on CPR : "<<endl;
for (int m = 0; m < N; m++)
{
    for(int n = 0; n < 3; n++)
        cout << d[m][n] << ' ';
    cout << endl;
}
	p = new double*[N];
	for (int i = 0; i < N; ++i) 
  		p[i] = new double[N];
	
	r = new double*[N];
	for (int i = 0; i < N; ++i) 
  		r[i] = new double[N];
	
	fh = new double*[N];
	for (int i = 0; i < N; ++i) 
  		fh[i] = new double[N];
	
	fy = new double*[N];
	for (int i = 0; i < N; ++i) 
  		fy[i] = new double[N];

	
        double sp1[N][N];
        int y[N][3];
        //double s[N][N];	


	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			sp1[i][j]=0;
		
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++){
		if (matrix[i + j*M]==1 && S[i][j]==0)
			sp1[i][j]=59;
		else if (matrix[i + j*M]==1 && S[i][j]==1)
			sp1[i][j]=1;
		}

	cout<<"The Sparse matrix S*: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<sp1[i][j]<<" ";
	}
	cout<<endl;
	}	
/* calculating the Y vector*/

	for(int i=0;i<N;i++){
            for(int j=0;j<3;j++){
		y[i][j]=0;
                for(int k=0;k<N;k++){
                   y[i][j]=y[i][j]+(sp1[i][k]*d[k][j]);
		}
            }
        }

	cout<<"Vector Y (H*d = y) : "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<3; j++)			
			cout<<y[i][j]<<" ";
		cout<<endl;
	}
	cout<<endl;

	int o;
	for (int i=0; i<N; i++)
		for (int j=0; j<3; j++){
			if(y[i][j]!=0 && y[i][j]<N+1){	
				int a= y[i][j];
				r[i][a-1]=1;}
			if(y[i][j]!=0 && y[i][j]>N+1){
				o= y[i][j]%59;
				if (o==0)
					o=y[i][j]/59;
				if (o>N)
					o=o/X;
				
				r[i][o-1]=-1;
			}
		}
	
	cout<<"The Sparse matrix H: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<r[i][j]<<" ";
	}
	cout<<endl;
	}
	/* calculating the flaws*/

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++){
			if (i == j){
				if (matrix[i + j*M]-r[i][j] > 0)
					fh[i][j]=S[i][j];
			}
			else
			{
				if (matrix[i + j*M]-r[i][j] > 0){
					fh[i][j]=S[i][j];
					fh[j][i]=S[j][i];
					}
			}
			//cout<<fh[i][j]<<" ";
		}
	
	
	cout<<"The F^h matrix: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<fh[i][j]<<" ";
	}
	cout<<endl;
	}
	
	for (int i=0; i<N; i++)
		for (int j=0; j<3; j++)
			for (int k=0; k<N; k++){
				if (S[i][k]==0 && (y[i][j]%59==0)){
					fy[i][j]=y[i][j];
				}
				else if (y[i][j]%59!=0) 
					fy[i][j]=0;	
			}				
	

	cout<<"The F^y matrix:: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<3; j++)			
			cout<<fy[i][j]<<" ";
		cout<<endl;
	}
	
	/*calculating the possibilities*/

	possibilities();
	cout<<"The P matrix: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<p[i][j]<<" ";
	}
	cout<<endl;
	}
	/* eliminating non symmetric elements from possibilities*/
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++){
			if (j>=i){
			if (p[i][j]!=p[j][i])		
				p[i][j]=p[j][i]=0;
		}
	}
		
	cout<<"The P^T matrix: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<p[i][j]<<" ";
	}
	cout<<endl;
	}

	/*S=So+p^T*/
	double si[N][N];
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++){		
				si[i][j]=S[i][j]+p[i][j];
	}	

	cout<<"The S matrix: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<si[i][j]<<" ";
	}
	cout<<endl;
	}
	/*##################*/
	double sp2[N][N];
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			sp2[i][j]=0;
		
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++){
		if (matrix[i + j*M]==1 && si[i][j]==65)
			sp2[i][j]=1;
		else if (matrix[i + j*M]==1 && si[i][j]==1)
			sp2[i][j]=1;
		else if (si[i][j]==65)
			sp2[i][j]=65;
		
	}
	cout<<"The Sp# matrix: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<sp2[i][j]<<" ";
	}
	cout<<endl;
	}

	/* calculating new H */
	
	double y1[N][X];
	for(int i=0;i<N;i++){
            for(int j=0;j<X;j++){
		y1[i][j]=0;
                for(int k=0;k<N;k++){
                   y1[i][j]=y1[i][j]+(sp2[i][k]*d1[k][j]);
		}
            }
        }

	cout<<"Vector Y1 (H*d = y) : "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<X; j++)			
			cout<<y1[i][j]<<" ";
		cout<<endl;
	}
	cout<<endl;

	double r1[N][N];
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++)		
			r1[i][j]=0;
	}
	int o1,o2;
	for (int i=0; i<N; i++)
		for (int j=0; j<X; j++){
			if(y1[i][j]!=0 && y1[i][j]<N+1)	{
				o2= (int) y1[i][j];		
				r1[i][o2-1]=1;}
			else if(y1[i][j]>N+1){
				o1= (int) y1[i][j] % 65;
				
				r1[i][o1-1]=2;	
					
			}	
			
		}
	
	cout<<"The New H matrix: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<r1[i][j]<<" ";
	}
	cout<<endl;
	}

	double sz[N][N];
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++)		
			sz[i][j]=0;
	}
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++){
			if (j>=i){
				int z = min(r1[i][j],r1[j][i]);
				if (z < 1)
					sz[i][j]=sz[j][i]=1;
			}
		}

	cout<<"The Sz matrix: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<sz[i][j]<<" ";
	}
	cout<<endl;
	}
	
	/*Final S*/
	
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++){	
				if (sz[i][j]==1){
					r1[i][j]=0;
					si[i][j]=0;}
				if (si[i][j]==65)
					si[i][j]=1;
				
	}
	cout<<endl;
	cout<<"The Final H matrix: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<r1[i][j]<<" ";
	}
	cout<<endl;
	}	
	cout<<endl;

	cout<<"The Final S matrix: "<<endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){		
			cout<<si[i][j]<<" ";
	}
	cout<<endl;
	}


	delete []  d1 ;
	delete []  d ;
	delete []  p ;
	delete []  r ;
  	delete []  fh ;
 	delete []  fy ;

	return 0;
}

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include<vector>
#include <cstdlib>

using namespace std;


vector<string> split(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

int main(int argc, char *argv[])
{

	string line;
	ifstream fin;
	ofstream fout;
	vector<string> tokens;
	if(argc > 1)
	{
		string fileName = argv[1];	
		if(fileName.find(".col") != string::npos)
		{
			fin.open(argv[1]);
			tokens = split(fileName,'.');
			fout.open((tokens[0]+".mtx").c_str());
		}
		else
		{
			cout<<"DIMACS format file only"<<endl;
			exit(0);
		}
		
	}
	else
	{
		cout<<"PLease provide the DIMACS file path"<<endl;
	}

	char first;
	do
	{
		getline(fin,line);
		first = line.at(0);
	}while(toupper(first)=='C');
	
	// cout<<tokens[2]<<endl;
	fout<<"%%MatrixMarket matrix coordinate pattern general"<<endl;

	tokens = split(line,' ');
	vector<string> tok = split(tokens[3],'\r'); // not needed in linux
	//fout<<tokens[3]<<" "<<tokens[2]<<" "<<atoi(tokens[3].c_str())*2<<endl;
	int M = atoi(tok[0].c_str());
	// int M = tokens[3];
	int N =  atoi(tokens[2].c_str());
	int NNZ = atoi(tokens[3].c_str())*2;
	fout<<M<<" "<<N<<" "<<NNZ<<endl;
	for(int i=1; i<= M; i++)
	{
		getline(fin,line);
		tokens = split(line,' ');
		tok = split(tokens[2],'\r'); // not needed in linux
		fout<<i<<" "<<tokens[1]<<endl;
		fout<<i<<" "<<tok[0]<<endl;
	
	}
	return 0;
}
#ifndef _BUMPY_PREPARE_H_
#define _BUMPY_PREPARE_H_

#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdio.h>

using namespace std;

void stringSplit(string & s_in, int N, double * output, char delim = ' ')
{
	stringstream ss(s_in);
	int count = 0;
	string token;
	while (getline(ss, token, delim)) 
	{
        output[count] = atof(token.data());
		count += 1;
    }
	return;
}

void initializaSavePath(string & filename, char ** argv)
{
    string space = "_";
	string fileext = ".out";
    filename.append(argv[1]);
	filename.append(space);
	filename.append(argv[2]);
	filename.append(space);
	filename.append(argv[3]);
	filename.append(space);
	filename.append(argv[4]);
    filename.append(space);
	filename.append(argv[5]);
    filename.append(fileext);
}

void initializaLoadPath(string & loadname, char ** argv)
{
    double mu = atof(argv[2]);
    double phi_target = atof(argv[4]);
    string space = "_";
	string fileext = ".out";
    string loadfileext = ".txt";
    loadname.append(argv[1]);
	loadname.append(space);
    ostringstream muS;
    muS << fixed << setprecision(4) << mu;
	loadname.append(muS.str());
	loadname.append(space);
	loadname.append(argv[3]);
	loadname.append(space);
	ostringstream phiS;
    phiS << fixed << setprecision(4) << phi_target;
	loadname.append(phiS.str());
    loadname.append(space);
	loadname.append(argv[5]);
    loadname.append(loadfileext);
}

void loadInitialConfiguration(string & loadname, int & n_l, int & n_s, int & Ns
, int & Nl, double & Lx, double & Ly, double & dt, double * x, double * y, double * th, double * Dn
, double * R_eff, double * m, double * I, double * r_shape_l, double * th_shape_l, double * r_shape_s, double * th_shape_s)
{
	int N;
	string line;
	ifstream myfile (loadname);
	if (myfile.is_open())
	{
		getline (myfile,line);
		N = atoi(line.data());
		getline (myfile,line);
		n_l = atoi(line.data());
		getline (myfile,line);
		n_s = atoi(line.data());
		getline (myfile,line);
		Nl = atoi(line.data());
		getline (myfile,line);
		Ns = atoi(line.data());
		getline (myfile,line);
		Lx = atof(line.data());
		getline (myfile,line);
		Ly = atof(line.data());
		getline (myfile,line);
		stringSplit(line, N, x);
		getline (myfile,line);
		stringSplit(line, N, y);
		getline (myfile,line);
		stringSplit(line, N, th);
		getline (myfile,line);
		stringSplit(line, N, Dn);
		getline (myfile,line);
		stringSplit(line, N, R_eff);
		getline (myfile,line);
		stringSplit(line, N, m);
		getline (myfile,line);
		stringSplit(line, N, I);
		getline (myfile,line);
		stringSplit(line, n_l, r_shape_s);
		getline (myfile,line);
		stringSplit(line, n_l, th_shape_s);
		getline (myfile,line);
		stringSplit(line, n_l, r_shape_l);
		getline (myfile,line);
		stringSplit(line, n_l, th_shape_l);
		getline (myfile,line);
		dt = atof(line.data());
		myfile.close();
	}

	
}

#endif
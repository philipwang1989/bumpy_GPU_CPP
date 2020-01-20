#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstring>
#include "bumpy_generate.hpp"
#include "md_print.hpp"
#include "md.hpp"

using namespace std;

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592;

int main(int argc, char **argv)
{
    // TO COMPILE (MAC or LINUX?)
    // clang++ -std=c++14 -stdlib=libc++ Bumpyjam_v0.cpp 
    // TO COMPILE (WINDOWS w/ VS_community 2017+)
    // cl -EHsc -std:c++14 Bumpyjam_v0.cpp 
    // TO COMPILE (WINDOWS w/ clang - not perfect, only works for -m32)
    // clang++ -std=c++14 -m32 -stdlib=libc++ Bumpyjam_v0.cpp

    // TO RUN (in the pseudo-sense...)
    // ./compiled_program Nc n /mu phi_target seed
    // you can just do the following:
    // ./compiled_program 6 16 1 1 1

    //Save output file
    string space = "_";
	string fileext = ".out";

    string debugfile = "";
    string debugfileext = ".xy";

    // string filename = "E:/Dropbox/Yale/C++/Bumpy/output_"; // Windows
    string filename = "/Users/philipwang/Dropbox/Yale/C++/Bumpy/output_"; // Mac
    filename.append(argv[1]);
	filename.append(space);
	filename.append(argv[2]);
	filename.append(space);
	filename.append(argv[3]);
	filename.append(space);
	filename.append(argv[4]);
    filename.append(space);
	filename.append(argv[5]);
    debugfile.append(filename);
    filename.append(fileext);
    debugfile.append(debugfileext);
    char *filechar = new char[filename.length() + 1];
	strcpy(filechar, filename.c_str());
	FILE *out = fopen(filechar, "w+");
    printf("File char is: %s\n", filechar);

    char *debugfilechar = new char[debugfile.length() + 1];
    strcpy(debugfilechar, debugfile.c_str());
    FILE *debug = fopen(debugfilechar, "w+");

    int Nc;// = 6;
    int n;// = 16; // n-mers
    double mu;// = 1.0;
    double phi_target;// = 1.00;
    int seed;// = 1;

    Nc = atoi(argv[1]);
    n = atoi(argv[2]);
    mu = atof(argv[3]);
    phi_target = atof(argv[4]);
    seed = atoi(argv[5]);

    double temp, U, dphi, dG, Utol, tol, gam, phitot;
    int count, count_max, C;
    int N = n * Nc;

    mt19937 generator(seed);

    // create particles
    int Nsmall = Nc / 2;
    int Nbig = Nc / 2;

    // create basic nmers
    double * xval = new double[n];
    double * yval = new double[n];
    double * lengths = new double[n];
    space_equally_circle(n, xval, yval, lengths);
    // printBasicNmers(n, xval, yval, lengths);

    double l = average(n, lengths);
    double D = (l/2) * sqrt(1 + mu * mu) / mu;
    for (int i=0;i<n;i++)
    {
        xval[i] = xval[i] / D;
        yval[i] = yval[i] / D;
        lengths[i] = lengths[i] / D;
    }
    double rad = 1.0 / D;
    D = 1.0;
    double E0 = 1.0;
    double K = 2.0 * E0 * D * D;
    double G = 1.4;

    // set spatial quantities
    double Lx = round(10.0 * G * sqrt(Nc) * (D + rad));
    double Ly = round(10.0 * G * sqrt(Nc) * (D + rad));

    int long Nt = 1e6;
    double N_per_coll = 100.0;
    double dt = 2.0 * pi * sqrt(1.0 / K) / N_per_coll;

    // set intermediate MD arrays
    double * Dn0 = new double[Nc];
    double * R_eff0 = new double[Nc];
    double * m0 = new double[Nc];
    double * I0 = new double[Nc]; // moment of inertia of disks along center = 1/2 * m * r * r    
    double * x0 = new double[Nc];
    double * y0 = new double[Nc];
    double * th0 = new double[Nc];
    double * vx0 = new double[Nc];
    double * vy0 = new double[Nc];
    double * w0 = new double[Nc];
    double * ax_old0 = new double[Nc];
    double * ay_old0 = new double[Nc];
    double * alph_old0 = new double[Nc];
    // build nmers
    double * x_shape0 = new double [N];
    double * y_shape0 = new double [N];
    double * r_shape0 = new double [N];
    double * th_shape0 = new double [N];

    // set MD arrays
    double * Dn = new double[Nc];
    double * R_eff = new double[Nc];
    double * m = new double[Nc];
    double * I = new double[Nc]; // moment of inertia of disks along center = 1/2 * m * r * r    
    double * x = new double[Nc];
    double * y = new double[Nc];
    double * th = new double[Nc];
    double * vx = new double[Nc];
    double * vy = new double[Nc];
    double * w = new double[Nc];
    double * ax_old = new double[Nc];
    double * ay_old = new double[Nc];
    double * alph_old = new double[Nc];

    // local MD arrays
    double * ax = new double[Nc];
    double * ay = new double[Nc];
    double * alph = new double[Nc];

    // build nmers
    double * x_shape = new double [N];
    double * y_shape = new double [N];
    double * r_shape = new double [N];
    double * th_shape = new double [N];
    // set nmer local variables
    long int MCpoints = 1e5;
    double * xlist = new double[n];
    double * ylist = new double[n];
    // double * Rlist = new double[n];
    double * x_mc = new double[MCpoints];
    double * y_mc = new double[MCpoints];
    int * in_mc = new int[MCpoints];

    // stress tensor
    double * stress = new double [4];
    double P, Ptol;

    for(int i = 0; i < Nc / 2; i++) R_eff[i] = D + rad;
	for(int i = Nc / 2; i < Nc; i++) R_eff[i] = G * (D + rad);
    for(int i = 0; i < Nc / 2; i++) Dn[i] = D; // small
	for(int i = Nc / 2; i < Nc; i++) Dn[i] = G; // big

    uniform_real_distribution<double> distribution(0.0, Lx); // Generate number between 0, Lx    
    // seed random initial particle positions
	/*
    do
	{
		for(int i = 0; i < Nc; i++) x[i] = distribution(generator);
        for(int i = 0; i < Nc; i++) y[i] = distribution(generator);
	}while(anytouch(Nc, x, y, R_eff, Lx, Ly));//re-seed until no particles touch
    */

    for(int nn=0;nn<Nc;++nn){
        // cout << nn <<endl;
        x[nn] = distribution(generator);
        y[nn] = distribution(generator);
        for(int mm=0;mm<nn;++mm){
            double dx,dy,Dnm,dnm2;
            dx = x[mm] - x[nn];
            dx = dx - round(dx/Lx)*Lx;
            dy = y[mm] - y[nn];
            dy = dy - round(dy/Ly)*Ly;
			//sum of 2 radii
			Dnm = R_eff[nn] + R_eff[mm];
            dnm2 = dx*dx + dy*dy;
            if(dnm2<Dnm*Dnm)
            {
                // cout << nn << "," <<mm<<endl;
                --nn;
                break;
            }
        }
    }

    printf("Particles ready!\n");

    nmersBuild(n, Nc, G, xval, yval, x_shape, y_shape, r_shape, th_shape);
    nmersProperty(n, Nc, MCpoints, xlist, ylist, x_mc, y_mc, in_mc, m, I, x_shape, y_shape, G/2.0, D/2.0);
    // printNmers(n, Nc, x_shape, y_shape, r_shape, th_shape);
    // printf("%e\n",-INFINITY);
    // printParticles(Nc, m, I);
    temp = 0;
    for(int i = 0; i < Nc; i++) temp += m[i];
    phitot = temp / (Lx * Ly);
    // printf("Lx=%e,Lx*Ly=%e,phi_ini=%e.\n",Lx,Lx*Ly,phitot);
    
    for (int i=0;i<Nc;i++)
    {
        temp = -INFINITY;
        for (int j=0;j<n;j++) if (r_shape[i*n+j] > temp) temp = r_shape[i*n+j]; // get max
        R_eff[i] = temp + Dn[i] / 2;
    }

    // update dt
    // double m_mean = average(Nc, m);
    // dt = 2.0 * pi * sqrt(m_mean / K) / N_per_coll;

    // First print out a file to plot the unjam states

    double * Fx = new double[Nc];
    double * Fy = new double[Nc];
    double * T = new double[Nc];
    int * Cn = new int[Nc];

    count_max = 1e6;
    gam = 0;
    Utol = 1e-14;
    Ptol = 1e-7;
    tol = 1e-7;
    dphi = 1e-3;
    count = 0;
    C = 0;
    U = 0;
    P = 0;

    double Kmax = 0;

    // Print initial states
    fprintf(debug, "%d\n", Nc);
    fprintf(debug, "%d\n", n);
    fprintf(debug, "%d\n", N);
    fprintf(debug, "%3.4f\n", mu);
    fprintf(debug, "%03.16f\n", Lx);
    fprintf(debug, "%03.16f\n", Ly);
    fprintf(debug, "%03.20f\n", gam);
    for (int j = 0; j < Nc; j++)
	{
		fprintf(debug, "%03.20f", Dn[j]); // this is Dn!!!
		fprintf(debug, ", %03.20f", x[j]);
		fprintf(debug, ", %03.20f", y[j]);
        fprintf(debug, ", %03.20f", th[j]);
        for (int k=0;k<n;k++) fprintf(debug, ", %03.20f", r_shape[j*n+k]);
        for (int k=0;k<n;k++) fprintf(debug, ", %03.20f", th_shape[j*n+k]);
        fprintf(debug, ", %d", Cn[j]);
		fprintf(debug, "\n");
	}

    fclose(debug);

    // Compression jamming
    while (P < Ptol || P > (1.0 + tol) * Ptol)
    // need to fix bug here
    // while (U < Utol || U > 2.0 * Utol)
    {
        // if (U < Utol)
        if (P < Ptol)
        {
            dphi = fabs(dphi);
            // copy all inputs
            copy(Dn, Dn + Nc, Dn0);
            copy(R_eff, R_eff + Nc, R_eff0);
            copy(m, m + Nc, m0);
            copy(I, I + Nc, I0);
            copy(x, x + Nc, x0);
            copy(y, y + Nc, y0);
            copy(th, th + Nc, th0);
            copy(vx, vx + Nc, vx0);
            copy(vy, vy + Nc, vy0);
            copy(w, w + Nc, w0);
            copy(ax_old, ax_old + Nc, ax_old0);
            copy(ay_old, ay_old + Nc, ay_old0);
            copy(alph_old, alph_old + Nc, alph_old0);
            copy(x_shape, x_shape + N, x_shape0);
            copy(y_shape, y_shape + N, y_shape0);
            copy(r_shape, r_shape + N, r_shape0);
            copy(th_shape, th_shape + N, th_shape0);
        }
        else if (P > (1.0 + tol) * Ptol)
        // else if (U > 2 * Utol)
        // else if (U > 2 * Utol && C >= (2 * Nc - 1))
        {
            dphi = -fabs(dphi) / 2;
            // copy all inputs
            copy(Dn0, Dn0 + Nc, Dn);
            copy(R_eff0, R_eff0 + Nc, R_eff);
            copy(m0, m0 + Nc, m);
            copy(I0, I0 + Nc, I);
            copy(x0, x0 + Nc, x);
            copy(y0, y0 + Nc, y);
            copy(th0, th0 + Nc, th);
            copy(vx0, vx0 + Nc, vx);
            copy(vy0, vy0 + Nc, vy);
            copy(w0, w0 + Nc, w);
            copy(ax_old0, ax_old0 + Nc, ax_old);
            copy(ay_old0, ay_old0 + Nc, ay_old);
            copy(alph_old0, alph_old0 + Nc, alph_old);
            copy(x_shape0, x_shape0 + N, x_shape);
            copy(y_shape0, y_shape0 + N, y_shape);
            copy(r_shape0, r_shape0 + N, r_shape);
            copy(th_shape0, th_shape0 + N, th_shape);
        }
        U = bumpy_2D_comp_VV(Nc, N, n, gam, Dn, r_shape, th_shape, m, I, R_eff, x, y, th, Fx, Fy, T, Cn, vx, vy, w, ax, ay, alph, ax_old, ay_old, alph_old, dt, Nt, K, Lx, Ly, Utol, dphi);
        P = bumpy_2D_stress(Nc, N, n, x, y, th, r_shape, th_shape, stress, K, R_eff, Dn, Lx, Ly, gam);
        Kmax = getEk(Nc, m, I, vx, vy, w);
        temp = 0;
        count += 1;
        // get C
        C = bumpy_2D_shrjam_getC(Nc, N, n, x, y, th, r_shape, th_shape, K, R_eff, Dn, Lx, Ly, gam);
        for(int i = 0; i < Nc; i++) temp += m[i];
        phitot = temp / (Lx * Ly);
        if (count % 100 == 0) printf("Step %d, phi=%1.7f, dphi=%e, C=%d, U/K/N=%e, P=%e\n", count, phitot, dphi, C, U, P);
        if (count > count_max) break;
        // m_mean = average(Nc, m);
        // dt = 2.0 * pi * sqrt(m_mean / K) / N_per_coll;
    }

    printf("U/K/N=%e, Kmax=%e, P=%e\n", U, Kmax, P);

    // print final compression jammed state
    fprintf(out, "%d\n", Nc);
    fprintf(out, "%d\n", n);
    fprintf(out, "%d\n", N);
    fprintf(out, "%3.4f\n", mu);
    fprintf(out, "%03.16f\n", Lx);
    fprintf(out, "%03.16f\n", Ly);
    fprintf(out, "%03.20f\n", gam);
    for (int j = 0; j < Nc; j++)
	{
		// fprintf(out, "%3d, ", j);
		fprintf(out, "%03.20f", Dn[j]); // this is Dn!!!
		fprintf(out, ", %03.20f", x[j]);
		fprintf(out, ", %03.20f", y[j]);
        fprintf(out, ", %03.20f", th[j]);
        for (int k=0;k<n;k++) fprintf(out, ", %03.20f", r_shape[j*n+k]);
        for (int k=0;k<n;k++) fprintf(out, ", %03.20f", th_shape[j*n+k]);
        fprintf(out, ", %d", Cn[j]);
		fprintf(out, "\n");
	}
    printf("\n");
    
    fclose(out);

    delete[] xval, yval, lengths;
    delete[] Dn0, R_eff0, m0, I0, x0, y0, th0, vx0, vy0, w0, ax_old0, ay_old0, alph_old0;
    delete[] x_shape0, y_shape0, r_shape0, th_shape0;
    delete[] ax, ay, alph;
    delete[] Dn, R_eff, m, I, x, y, th, vx, vy, w, ax_old, ay_old, alph_old;
    delete[] x_shape, y_shape, r_shape, th_shape;
    
    delete[] xlist, ylist, x_mc, y_mc, in_mc;

    return 0;
}
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstring>
#include "bumpy_generate.hpp"

using namespace std;

// void space_equally_circle(int n, double * xval, double * yval, double * lengths)
void space_equally_circle(int n, double xval[], double yval[], double lengths[])
{
    // This function finds the equally spaced point on a unit circle of radius 1
    double circumference = 2.0 * M_PI;
    double dtheta = circumference / n;
    for (int i=0;i<n;i++)
    {
        xval[i] = cos(dtheta * i);
        yval[i] = sin(dtheta * i);
    }
    double diff = sqrt((xval[1] - xval[0]) * (xval[1] - xval[0]) + (yval[0] - yval[1]) * (yval[0] - yval[1]));
    for (int i=0;i<n;i++)
    {
        lengths[i] = diff;
    }
}

void nmersBuild(int n, int Nc, double G, double * xval, double * yval, double * x_shape, double * y_shape, double * r_shape, double * th_shape)
{
    // fill the shape arrays
    for (int j=0;j<Nc;j++)
    {
        for (int i=0;i<n;i++)
        {
            x_shape[j*n+i] = xval[i];
            y_shape[j*n+i] = yval[i];
        }
    }
    // scale with particle size
    for (int j=Nc/2;j<Nc;j++)
    {
        for (int i=0;i<n;i++)
        {
            x_shape[j*n+i] = x_shape[j*n+i] * G;
            y_shape[j*n+i] = y_shape[j*n+i] * G;
        }
    }
    // fill r and th array
    for (int j=0;j<Nc;j++)
    {
        for (int i=0;i<n;i++)
        {
            r_shape[j*n+i] = sqrt(x_shape[j*n+i] * x_shape[j*n+i] + y_shape[j*n+i] * y_shape[j*n+i]);
            th_shape[j*n+i] = atan2(y_shape[j*n+i],x_shape[j*n+i]);
        }
    }
    // recalculate x,y array
    for (int j=0;j<Nc;j++)
    {
        for (int i=0;i<n;i++)
        {
            x_shape[j*n+i] = r_shape[j*n+i] * cos(th_shape[j*n+i]);
            y_shape[j*n+i] = r_shape[j*n+i] * sin(th_shape[j*n+i]);
        }
    }

}

int inpolygon(int nvert, double * vertx, double * verty, double testx, double testy)
{
    int i, j, c = 0;
    for (i = 0, j = nvert-1; i < nvert; j = i++) 
    {
        if ( ((verty[i]>testy) != (verty[j]>testy)) &&
        (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
        c = !c;
    }
    return c;
}

void inpart(int n ,long int MCpoints, int * in, double * x, double * y, double * xp, double * yp, double R)
{
    double R2 = R * R;
    bool tempbool;
    // inpolygon
    for (long int i=0;i<MCpoints;i++)
    {
        in[i] = inpolygon(n,xp,yp,x[i],y[i]);
    }
    for (int j=0;j<n;j++)
    {
        for (long int i=0;i<MCpoints;i++)
        {
            tempbool = ((x[i] - xp[j]) * (x[i] - xp[j]) + (y[i] - yp[j]) * (y[i] - yp[j])) < R2;
            in[i] = in[i] | tempbool;
        }
    }
}

void nmersProperty(int n, int Nc, int long MCpoints, double * xlist, double * ylist, double * x, double * y, int * in, double * m, double * I, double * x_shape, double * y_shape, double Rl, double Rs)
{
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);

    // long int MCpoints = 1e5;
    // double * xlist = new double[n];
    // double * ylist = new double[n];
    // double * Rlist = new double[n];
    // double * x = new double[MCpoints];
    // double * y = new double[MCpoints];
    // int * in = new int[MCpoints];
    double Lx_mc, Ly_mc, xmin, ymin, Atemp, Itemp, temp, count, sum, xcm, ycm, xsum, ysum, temp_max, temp_min;

    // first deal with small particles
    for (int i=0;i<n;i++)
    {
        xlist[i] = x_shape[i];
        ylist[i] = y_shape[i];
    }
    temp_max = -INFINITY;
    temp_min = INFINITY;
    for (int i=0;i<n;i++)
    {
        if ((xlist[i] + Rs) > temp_max) temp_max = (xlist[i] + Rs);
        if ((xlist[i] - Rs) < temp_min) temp_min = (xlist[i] - Rs);
    }
    Lx_mc = temp_max - temp_min;
    xmin = temp_min;
    temp_max = -INFINITY;
    temp_min = INFINITY;
    for (int i=0;i<n;i++)
    {
        if ((ylist[i] + Rs) > temp_max) temp_max = (ylist[i] + Rs);
        if ((ylist[i] - Rs) < temp_min) temp_min = (ylist[i] - Rs);
    }
    Ly_mc = temp_max - temp_min;
    ymin = temp_min;
    for (long int i=0;i<MCpoints;i++)
    {
        x[i] = Lx_mc * distribution(generator) + xmin;
        y[i] = Ly_mc * distribution(generator) + ymin;
    }
    inpart(n, MCpoints, in, x, y, xlist, ylist, Rs);
    // Area
    sum = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) sum += 1.0;
    Atemp = (sum / double(MCpoints)) * Lx_mc * Ly_mc;
    // CM
    xsum = 0;
    ysum = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) xsum += x[i];
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) ysum += y[i];
    xcm = xsum / sum;
    ycm = ysum / sum;
    // I
    temp = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) temp += ((x[i] - xcm)*(x[i] - xcm) + (y[i] - ycm)*(y[i] - ycm));
    Itemp = Atemp * temp / sum;
    for (int i=0;i<Nc/2;i++)
    {
        m[i] = Atemp;
        I[i] = Itemp;
    }
    // done with small particles

    // then deal with large particles
    for (int i=0;i<n;i++)
    {
        xlist[i] = x_shape[(Nc/2)*n+i];
        ylist[i] = y_shape[(Nc/2)*n+i];
    }
    temp_max = -INFINITY;
    temp_min = INFINITY;
    for (int i=0;i<n;i++)
    {
        if ((xlist[i] + Rl) > temp_max) temp_max = (xlist[i] + Rl);
        if ((xlist[i] - Rl) < temp_min) temp_min = (xlist[i] - Rl);
    }
    Lx_mc = temp_max - temp_min;
    xmin = temp_min;
    temp_max = -INFINITY;
    temp_min = INFINITY;
    for (int i=0;i<n;i++)
    {
        if ((ylist[i] + Rl) > temp_max) temp_max = (ylist[i] + Rl);
        if ((ylist[i] - Rl) < temp_min) temp_min = (ylist[i] - Rl);
    }
    Ly_mc = temp_max - temp_min;
    ymin = temp_min;
    for (long int i=0;i<MCpoints;i++)
    {
        x[i] = Lx_mc * distribution(generator) + xmin;
        y[i] = Ly_mc * distribution(generator) + ymin;
    }
    inpart(n, MCpoints, in, x, y, xlist, ylist, Rl);
    // Area
    sum = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) sum += 1.0;
    Atemp = (sum / double(MCpoints)) * Lx_mc * Ly_mc;
    // CM
    xsum = 0;
    ysum = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) xsum += x[i];
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) ysum += y[i];
    xcm = xsum / sum;
    ycm = ysum / sum;
    // I
    temp = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) temp += ((x[i] - xcm)*(x[i] - xcm) + (y[i] - ycm)*(y[i] - ycm));
    Itemp = Atemp * temp / sum;
    for (int i=Nc/2;i<Nc;i++)
    {
        m[i] = Atemp;
        I[i] = Itemp;
    }

    // delete[] xlist, ylist, x, y, in;
}

void bumpy_2D_create_particles(double * Dn, double rad, int Nc, int n_org, double mu)
{
    
}
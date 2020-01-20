#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstring>
#include "md_print.hpp"

using namespace std;

void printBasicNmers(int n, double * xval, double * yval, double * lengths)
{
    printf("Basic nmers:\n");
    printf("i, xval, yval, lengths\n");
    for (int i=0;i<n;i++)
    {
        printf("%d, %1.6f, %1.6f, %1.6f\n",i,xval[i],yval[i],lengths[i]);
    }
}

void printNmers(int n, int Nc, double * x_shape, double * y_shape, double * r_shape, double * th_shape)
{
    printf("x_shape, y_shape, r_shape, th_shape\n");
    for (int j=0;j<Nc;j++)
    {
        for (int i=0;i<n;i++)
        {
            printf("%1.6f, %1.6f, %1.6f, %1.6f\n",x_shape[j*n+i],y_shape[j*n+i],r_shape[j*n+i],th_shape[j*n+i]);
        }
        printf("\n");
    }

}

void printParticles(int Nc, double * m, double * I)
{
    printf("m, I\n");
    for (int j=0;j<Nc;j++)
    {
        printf("%1.6f, %1.6f\n",m[j],I[j]);
    }
}



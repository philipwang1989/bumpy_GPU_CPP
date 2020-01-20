#ifndef _MD_PRINT_H_
#define _MD_PRINT_H_

#include "md_print.cpp"

void printBasicNmers(int n, double * xval, double * yval, double * lengths);
void printNmers(int n, int Nc, double * x_shape, double * y_shape, double * r_shape, double * th_shape);
void printParticles(int Nc, double * m, double * I);

#endif
#ifndef _MDFUNC_H_
#define _MDFUNC_H_

#include "md.cpp"

bool anytouch(int N, double * x, double * y, double * R_eff, double Lx, double Ly);
double average(int N, double * numList);
double arrminf(int N, double * arr);
double arrmaxf(int N, double * arr);
double arrsumf(int N, double * arr);
void CMVelocityZeroing(int N, double * m, double * vx);
void VV_pos_integration(int N, double * x, double * vx, double * ax_old, double dt, double half_dt2);
void VV_vel_integration(int N, double * vx, double * ax, double * ax_old, double dt);
void getAcceleration(int N, double * Fx, double * ax, double * m);
void copyAcceleration(int N, double * ax, double * ax_old);
void correctionFIRE(int Nc, double * vx, double * Fx, double a);
void NonContactZeroing(int N, double * vx, int * Cn);
double getEk(int N, double * m, double * I, double * vx, double * vy, double * w);
double bumpy_2D_shrjam_getU(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape, double K, double * R_eff, double * Dn, double Lx, double Ly, double gam);
int bumpy_2D_shrjam_getC(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape, double K, double * R_eff, double * Dn, double Lx, double Ly, double gam);
double bumpy_2D_stress(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape, double * stress, double K, double * R_eff, double * Dn, double Lx, double Ly, double gam);
double bumpy_2D_Force(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape, double * Fx, double * Fy, double * T, int * Cn, double K, double * R_eff, double * Dn, double Lx, double Ly, double gam);
double bumpy_2D_comp_VV(int Nc, int N, int n, double gam, double * Dn, double * r_shape, double * th_shape, double * m, double * I, double * R_eff, double * x, double * y, double * th, double * Fx, double * Fy, double * T, int * Cn, double * vx, double * vy, double * w, double * ax, double * ay, double * alph, double * ax_old, double * ay_old, double * alph_old, double dt, long int Nt, double K, double Lx, double Ly, double Utol, double dphi);

#endif
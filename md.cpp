#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstring>
#include "md.hpp"

using namespace std;

// check if any particles are touching
bool anytouch(int N, double * x, double * y, double * R_eff, double Lx, double Ly)
{
    double dx, dy, Dnm, dnm2;    
	for(int nn = 0; nn < N-1; nn++)
	{
		for(int mm = nn + 1; mm < N; mm++)
		{
			//vector pointing from particle i to j
            dx = x[mm] - x[nn];
            dx = dx - round(dx/Lx)*Lx;
            dy = y[mm] - y[nn];
            dy = dy - round(dy/Ly)*Ly;
			//sum of 2 radii
			Dnm = R_eff[nn] + R_eff[mm];
            dnm2 = dx*dx + dy*dy;
			//they touch if the distance between them is less than Dnm square
			if(dnm2 < Dnm * Dnm) return true;
		}
	}
	return false;
}

double average(int N, double * numList)
{
    double average;
    double temp = 0.0;
    for (int i=1;i<N;i++)
    {
        temp += numList[i];
    }
    average = temp / N;
    return average;
}

double arrminf(int N, double * arr)
{
    double temp_min = INFINITY;
    for (int i=0;i<N;i++) if (arr[i] < temp_min) temp_min = arr[i];
    return temp_min;
}

double arrmaxf(int N, double * arr)
{
    double temp_max = -INFINITY;
    for (int i=0;i<N;i++) if (arr[i] > temp_max) temp_max = arr[i];
    return temp_max;
}

double arrsumf(int N, double * arr)
{
    double temp_sum = 0;
    for (int i=0;i<N;i++) temp_sum += arr[i];
    return temp_sum;
}
// Set CM velocity to 0
void CMVelocityZeroing(int N, double * m, double * vx)
{
    double M = 0;
    double P = 0;
    for (int i = 0; i < N; i++)
    {
        M += m[i];
        P += (m[i] * vx[i]);
    }
    double CM_v = P / M;
    for (int i = 0; i < N; i++)
    {
        vx[i] = vx[i] - CM_v;
    }
}

// Velocity Verlet position integration
void VV_pos_integration(int N, double * x, double * vx, double * ax_old, double dt, double half_dt2)
{
    for (int i = 0; i < N; i++)
    {
        x[i] = x[i] + vx[i] * dt + ax_old[i] * half_dt2;
    }
}

// Velocity Verlet velocity integration
void VV_vel_integration(int N, double * vx, double * ax, double * ax_old, double dt)
{
    for (int i = 0; i < N; i++)
    {
        vx[i] = vx[i] + (ax[i] + ax_old[i]) * dt / 2.0;
    }
}

void getAcceleration(int N, double * Fx, double * ax, double * m)
{
    for (int i = 0; i < N; i++) ax[i] = Fx[i] / m[i];
}

void copyAcceleration(int N, double * ax, double * ax_old)
{
    for (int i = 0; i < N; i++) ax_old[i] = ax[i];
}

void correctionFIRE(int Nc, double * vx, double * Fx, double a)
{
    double normF, normv;
    normF = 0;
    normv = 0;
    for (int i=0;i<Nc;i++)
    {
        normF += Fx[i] * Fx[i];
        normv += vx[i] * vx[i];
    }
    normF = sqrt(normF);
    normv = sqrt(normv);
    for (int i=0;i<Nc;i++)
    {
        vx[i] = (1.0 - a) * vx[i] + a * (Fx[i]/normF) * normv;
    }
}

// Set noncontacting particles to zero velocity
void NonContactZeroing(int N, double * vx, int * Cn)
{
    for (int i = 0; i < N; i++)
    {
        if (Cn[i] == 0)
        {
            vx[i] = 0.0;
        }
    }
}

// Calculate max Ek
double getEk(int N, double * m, double * I, double * vx, double * vy, double * w)
{
    double * K = new double[N];
    for (int i = 0; i < N; i++) K[i] = 0.0;
    for (int i = 0; i < N; i++)
    {
        K[i] += (m[i] * vx[i] * vx[i] + m[i] * vy[i] * vy[i] + I[i] * w[i] * w[i])/2;
    }
    double Kmax = 0.0;
    for (int i = 0; i < N; i++) if (K[i] > Kmax) Kmax = K[i];

    return Kmax;
}

double bumpy_2D_shrjam_getU(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape, double K, double * R_eff, double * Dn, double Lx, double Ly, double gam)
{
    double U, dx, dy, dxx, dyy, im, Dnm, dnm, DDnm, ddnm, rymmm, rynnn, rxmmm, rxnnn;
    U = 0;
    int nn, mm, nnn, mmm;
    for (nn=0;nn<Nc;nn++)
    {
        for (mm=nn+1;mm<Nc;mm++)
        {
            dy = y[mm] - y[nn];
            im = round(dy/Ly);
            dy = dy - im * Ly;
            Dnm = R_eff[nn] + R_eff[mm];
            if (abs(dy) < Dnm)
            {
                dx = x[mm] - x[nn];
                dx = dx - round(dx/Lx-im*gam)*Lx - im*gam*Lx;
                dnm = sqrt(dx * dx + dy * dy);
                if (dnm < Dnm)
                {
                    for (nnn=0;nnn<n;nnn++)
                    {
                        for (mmm=0;mmm<n;mmm++)
                        {
                            rymmm = r_shape[mm*n+mmm] * sin(th[mm]+th_shape[mm*n+mmm]);
                            rynnn = r_shape[nn*n+nnn] * sin(th[nn]+th_shape[nn*n+nnn]);
                            dyy = (y[mm] + rymmm) - (y[nn]+rynnn);
                            im = round(dyy/Ly);
                            dyy = dyy-im*Ly;
                            DDnm = (Dn[nn]+Dn[mm])/2;
                            if (abs(dyy) < Dnm)
                            {
                                rxmmm = r_shape[mm*n+mmm] * cos(th[mm]+th_shape[mm*n+mmm]);
                                rxnnn = r_shape[nn*n+nnn] * cos(th[nn]+th_shape[nn*n+nnn]);
                                dxx = (x[mm] + rxmmm) - (x[nn]+rxnnn);
                                dxx = dxx-round(dxx/Lx-im*gam)*Lx-im*gam*Lx;
                                ddnm=sqrt(dxx * dxx + dyy * dyy);
                                if (ddnm < DDnm)
                                {
                                    U += K * ((DDnm - ddnm) * (DDnm - ddnm)/2);
                                }
                            }

                        }
                    }
                }

            }
        }
    }
    U = U / (K * Nc);
    return U;
}

int bumpy_2D_shrjam_getC(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape, double K, double * R_eff, double * Dn, double Lx, double Ly, double gam)
{
    double dx, dy, dxx, dyy, im, Dnm, dnm, DDnm, ddnm, rymmm, rynnn, rxmmm, rxnnn;
    int C = 0;
    int nn, mm, nnn, mmm;
    for (nn=0;nn<Nc;nn++)
    {
        for (mm=nn+1;mm<Nc;mm++)
        {
            dy = y[mm] - y[nn];
            im = round(dy/Ly);
            dy = dy - im * Ly;
            Dnm = R_eff[nn] + R_eff[mm];
            if (abs(dy) < Dnm)
            {
                dx = x[mm] - x[nn];
                dx = dx - round(dx/Lx-im*gam)*Lx - im*gam*Lx;
                dnm = sqrt(dx * dx + dy * dy);
                if (dnm < Dnm)
                {
                    for (nnn=0;nnn<n;nnn++)
                    {
                        for (mmm=0;mmm<n;mmm++)
                        {
                            rymmm = r_shape[mm*n+mmm] * sin(th[mm]+th_shape[mm*n+mmm]);
                            rynnn = r_shape[nn*n+nnn] * sin(th[nn]+th_shape[nn*n+nnn]);
                            dyy = (y[mm] + rymmm) - (y[nn]+rynnn);
                            im = round(dyy/Ly);
                            dyy = dyy-im*Ly;
                            DDnm = (Dn[nn]+Dn[mm])/2;
                            if (abs(dyy) < Dnm)
                            {
                                rxmmm = r_shape[mm*n+mmm] * cos(th[mm]+th_shape[mm*n+mmm]);
                                rxnnn = r_shape[nn*n+nnn] * cos(th[nn]+th_shape[nn*n+nnn]);
                                dxx = (x[mm] + rxmmm) - (x[nn]+rxnnn);
                                dxx = dxx-round(dxx/Lx-im*gam)*Lx-im*gam*Lx;
                                ddnm=sqrt(dxx * dxx + dyy * dyy);
                                if (ddnm < DDnm)
                                {
                                    C += 1;
                                }
                            }

                        }
                    }
                }

            }
        }
    }
    return C;
}

double bumpy_2D_stress(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape, double * stress, double K, double * R_eff, double * Dn, double Lx, double Ly, double gam)
{
    double P, dx, dy, dxx, dyy, im, Dnm, dnm, DDnm, ddnm, rymmm, rynnn, rxmmm, rxnnn, F;
    int nn, mm, nnn, mmm;

    for (int i=0;i<4;i++) stress[i] = 0.0;
    // Overlap detection
    for (nn=0;nn<Nc;nn++)
    {
        for (mm=nn+1;mm<Nc;mm++)
        {
            dy = y[mm] - y[nn];
            im = round(dy/Ly);
            dy = dy - im * Ly;
            Dnm = R_eff[nn] + R_eff[mm];
            if (abs(dy) < Dnm)
            {
                dx = x[mm] - x[nn];
                dx = dx - round(dx/Lx-im*gam)*Lx - im*gam*Lx;
                dnm = sqrt(dx * dx + dy * dy);
                if (dnm < Dnm)
                {
                    for (nnn=0;nnn<n;nnn++)
                    {
                        for (mmm=0;mmm<n;mmm++)
                        {
                            rymmm = r_shape[mm*n+mmm] * sin(th[mm]+th_shape[mm*n+mmm]);
                            rynnn = r_shape[nn*n+nnn] * sin(th[nn]+th_shape[nn*n+nnn]);
                            dyy = (y[mm] + rymmm) - (y[nn]+rynnn);
                            im = round(dyy/Ly);
                            dyy = dyy-im*Ly;
                            DDnm = (Dn[nn]+Dn[mm])/2;
                            if (abs(dyy) < Dnm)
                            {
                                rxmmm = r_shape[mm*n+mmm] * cos(th[mm]+th_shape[mm*n+mmm]);
                                rxnnn = r_shape[nn*n+nnn] * cos(th[nn]+th_shape[nn*n+nnn]);
                                dxx = (x[mm] + rxmmm) - (x[nn]+rxnnn);
                                dxx = dxx-round(dxx/Lx-im*gam)*Lx-im*gam*Lx;
                                ddnm = sqrt(dxx * dxx + dyy * dyy);
                                if (ddnm < DDnm)
                                {
                                    F = -K*(DDnm/ddnm-1);
                                    stress[0] -= F * dxx * dx;
                                    stress[1] -= 0.5 * F * (dxx*dy + dyy*dx);
                                    stress[2] -= 0.5 * F * (dyy*dx + dxx*dy);
                                    stress[3] -= F * dyy * dy;
                                }
                            }

                        }
                    }
                }

            }
        }
    }
    for (int i=0;i<4;i++) stress[i] /= (Lx * Ly);
    P = (stress[0] + stress[3]) / 2.0;
    return P;
}

double bumpy_2D_Force(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape,
double * Fx, double * Fy, double * T, int * Cn, double K, double * R_eff, double * Dn, double Lx, double Ly, double gam)
{
    double U, dx, dy, dxx, dyy, im, Dnm, dnm, DDnm, ddnm, rymmm, rynnn, rxmmm, rxnnn, F;
    U = 0;
    int nn, mm, nnn, mmm;
    // Zeroing Fx, Fy and T
    for (int i = 0; i < Nc; i++)
    {        
        Fx[i] = 0.0;
        Fy[i] = 0.0;
        T[i] = 0.0;
        Cn[i] = 0;
    }
    // Overlap detection
    for (nn=0;nn<Nc;nn++)
    {
        for (mm=nn+1;mm<Nc;mm++)
        {
            dy = y[mm] - y[nn];
            im = round(dy/Ly);
            dy = dy - im * Ly;
            Dnm = R_eff[nn] + R_eff[mm];
            if (abs(dy) < Dnm)
            {
                dx = x[mm] - x[nn];
                dx = dx - round(dx/Lx-im*gam)*Lx - im*gam*Lx;
                dnm = sqrt(dx * dx + dy * dy);
                if (dnm < Dnm)
                {
                    for (nnn=0;nnn<n;nnn++)
                    {
                        for (mmm=0;mmm<n;mmm++)
                        {
                            rymmm = r_shape[mm*n+mmm] * sin(th[mm]+th_shape[mm*n+mmm]);
                            rynnn = r_shape[nn*n+nnn] * sin(th[nn]+th_shape[nn*n+nnn]);
                            dyy = (y[mm] + rymmm) - (y[nn]+rynnn);
                            im = round(dyy/Ly);
                            dyy = dyy-im*Ly;
                            DDnm = (Dn[nn]+Dn[mm])/2;
                            if (abs(dyy) < Dnm)
                            {
                                rxmmm = r_shape[mm*n+mmm] * cos(th[mm]+th_shape[mm*n+mmm]);
                                rxnnn = r_shape[nn*n+nnn] * cos(th[nn]+th_shape[nn*n+nnn]);
                                dxx = (x[mm] + rxmmm) - (x[nn]+rxnnn);
                                dxx = dxx-round(dxx/Lx-im*gam)*Lx-im*gam*Lx;
                                ddnm = sqrt(dxx * dxx + dyy * dyy);
                                if (ddnm < DDnm)
                                {
                                    Cn[nn] += 1;
                                    Cn[mm] += 1;
                                    F = -K*(DDnm/ddnm-1);
                                    Fx[nn] += F * dxx;
                                    Fx[mm] -= F * dxx;
                                    Fy[nn] += F * dyy;
                                    Fy[mm] -= F * dyy;
                                    T[nn] += F * (rxnnn * dyy - rynnn * dxx);
                                    T[mm] -= F * (rxmmm * dyy - rymmm * dxx);
                                    U += K * ((DDnm - ddnm) * (DDnm - ddnm)/2);
                                }
                            }

                        }
                    }
                }

            }
        }
    }
    U = U / (K * Nc);
    return U;
}

double bumpy_2D_comp_VV(int Nc, int N, int n, double gam, double * Dn, double * r_shape, double * th_shape, double * m, double * I, double * R_eff, double * x, double * y, 
double * th, double * Fx, double * Fy, double * T, int * Cn, double * vx, double * vy, double * w, double * ax, double * ay, double * alph, double * ax_old, double * ay_old, double * alph_old, double dt, long int Nt, double K, double Lx, double Ly, double Utol, double dphi)
{
    // double * ax = new double[Nc];
    // double * ay = new double[Nc];
    // double * alph = new double[Nc];

    double phitot, rsc, U, im, dt2, half_dt2, Ek_trans, Ek_rot, Ek_tot;
    int i, j;
    long int nt;
    // zeroing initial acceleration
    for (i=0;i<Nc;i++)
    {
        ax[i] = 0.0;
        ay[i] = 0.0;
        alph[i] = 0.0;
    }
    dt2 = dt * dt;
    half_dt2 = dt2 / 2.0;
    phitot = arrsumf(Nc, m) / (Lx * Ly);
    rsc = sqrt(1.0 + dphi / phitot);
    // Grow
    for (i=0;i<Nc;i++)
    {
        Dn[i] = Dn[i] * rsc;
        R_eff[i] = R_eff[i] * rsc;
        m[i] = m[i] * rsc * rsc;
        I[i] = I[i] * rsc * rsc * rsc * rsc;
    }
    for (i=0;i<N;i++) r_shape[i] = r_shape[i] * rsc;
    phitot = arrsumf(Nc, m) / (Lx * Ly);

    U = bumpy_2D_Force(Nc,N,n,x,y,th,r_shape,th_shape,Fx,Fy,T,Cn,K,R_eff,Dn,Lx,Ly,gam);

    // FIRE COEFF. 2.0
    int NPP, NPN, ND, NPNM;
    double finc, fdec, astart, a, at, fa, dtmax, dtmin, P;
    bool fire = true;
    bool initialdelay = true;
    // bool fire = false;
    finc = 1.1;
    fdec = 0.5;
    astart = 0.25;
    a = astart;
    fa = 0.99;
    dtmax = 10.0 * dt;
    dtmin = 0.02 * dt;
    NPP = 0;
    NPN = 0;
    ND = 20;
    NPNM = 2000;

    for(int i=0;i<Nc;++i){
        vx[i] = 0;
        vy[i] = 0;
        w[i] = 0;
    }

    double B = 1.0;
    double Kmax;

    nt = 0;
    while (U > 0 && nt < Nt)
    {
        // FIRE 2.0
        P = 0.0;
        if (fire)
        {
            for (int i=0;i<Nc;i++)
            {
                P += (vx[i] * Fx[i] + vy[i] * Fy[i] + w[i] * T[i]);
            }
            if (P <= 0)
            {
                NPP = 0;
                NPN++;
                if(NPN>NPNM) break;
                if( !(initialdelay && nt<ND) )
                {
                    if(dt*fdec>=dtmin) dt=dt*fdec;
                    a = astart;
                }
                for(int i=0;i<Nc;++i)
                {
                    x[i] -= 0.5*dt*vx[i];
                    y[i] -= 0.5*dt*vy[i];
                    th[i] -= 0.5*dt*w[i];
                    vx[i] = 0;
                    vy[i] = 0;
                    w[i] = 0;
                }
            }
            else if (P > 0)
            {
                NPP++;
                NPN = 0;
                at = a;
                if(NPP>ND)
                {
                    if (dt * finc < dtmax) dt = dt * finc;
                    else dt = dtmax;
                    a = a * fa;
                }
            }
            correctionFIRE(Nc, vx, Fx, at);
            correctionFIRE(Nc, vy, Fy, at);
            correctionFIRE(Nc, w, T, at);
        }

        // //Verlet integration
        for(int i=0;i<Nc;++i)
        {
            vx[i] += 0.5*dt*Fx[i]/m[i];
            vy[i] += 0.5*dt*Fy[i]/m[i];
            w[i] += 0.5*dt*T[i]/I[i];
        }

        correctionFIRE(Nc, vx, Fx, at);
        correctionFIRE(Nc, vy, Fy, at);
        correctionFIRE(Nc, w, T, at);

        for(int i=0;i<Nc;++i)
        {
            x[i] += dt*vx[i];
            y[i] += dt*vy[i];
            th[i] += dt*w[i];
        }

        U = bumpy_2D_Force(Nc,N,n,x,y,th,r_shape,th_shape,Fx,Fy,T,Cn,K,R_eff,Dn,Lx,Ly,gam);

        for(int i=0;i<Nc;++i)
        {
            vx[i] += 0.5*dt*Fx[i]/m[i];
            vy[i] += 0.5*dt*Fy[i]/m[i];
            w[i] += 0.5*dt*T[i]/I[i];
        }

        // MOD
        for (int i=0;i<Nc;i++)
        {
            im = floor(y[i] / Ly);
            x[i] = fmod(x[i]-im*gam*Lx,Lx);
            if (x[i] < 0) x[i] += Lx;
            y[i] = fmod(y[i],Ly);
            if (y[i] < 0) y[i] += Ly;
            th[i] = fmod(th[i],2.0*M_PI);
            if (th[i] < 0) th[i] += 2.0*M_PI;
        }
        // zeroing
        CMVelocityZeroing(Nc, m, vx);
        CMVelocityZeroing(Nc, m, vy);

        Kmax = getEk(Nc, m, I, vx, vy, w);
        
        if (nt > 0 && Kmax < 1e-28) 
        {
            // printf("Break by Ek=%e at %ld.\n",Kmax,nt);
            // print("\n");
            break;
        }

        nt += 1;
    }

    // printf("U=%e.\n",U);
    // delete[] ax, ay, alph;
    return U;
}

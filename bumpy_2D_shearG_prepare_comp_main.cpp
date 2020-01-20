#include <iostream>
#include <cstring>
#include <random>
#include "bumpy_2D_prepare.hpp"
#include "debugUtility.hpp"

using namespace std;

// srun --pty -p gpu -c 2 -t 2:00:00 --gres=gpu:v100:1 bash
// module restore cuda10 
// nvcc -x cu nested_for_loop.cpp -o test.exe -gencode arch=compute_70,code=sm_70

/*
__global__ void grow_cluster_gpu(double * Dn, double * R_eff, double * rsc)
{
    int i = blockIdx.x;
    Dn[i] *= *rsc;
    R_eff[i] *= *rsc;
}

__global__ void correct_mass_gpu(int N, double * m, double * rsc)
{
    int i = blockIdx.x;
    if (i < 2*N) m[i] *= (*rsc * *rsc);
    else m[i] *= (*rsc * *rsc * *rsc * *rsc);
}

__global__ void grow_bumpy_gpu(double * r_shape, double * rsc)
{
    int i = blockIdx.x;
    r_shape[i] *= *rsc;
}

__global__ void updatePos(double * pos, double * vel, double * acc_old, double * dt)
{
    int i = blockIdx.x;
    pos[i] += (vel[i] * *dt + acc_old[i] * (*dt * *dt) * 0.5);
}

__global__ void zeroArr(double * arr)
{
    int i = blockIdx.x;
    arr[i] = 0.0;
}

__global__ void reduce(int s, double * arr)
{
	int i = blockIdx.x;
	if(i % (2*s) == 0 && i+s < gridDim.x) arr[i] += arr[i + s];
}

void sum_gpu(int N, double * arr)
{
    for(unsigned int s=1; s < N; s *= 2)
    {
        reduce<<<N,1>>>(s, arr);
    }
}

__global__ void bumpy_2D_getForce(int N, int n, double * pos, double * R_eff, double * Dn, double * r_shape, double * th_shape, double * force, double * U)
{
    // bumpy_2D_getForce<<<dim3(N,N,n*n),1>>>(N, n, pos, R_eff, Dn, r_shape, th_shape, force, U);
	int nn = blockIdx.x; // 0 to N-1
	int mm = blockIdx.y; // 0 to N-1
	int tid = blockIdx.z; // 0 to n*n-1
	int nnn = tid % n;
	int mmm = tid / n;
	// printf("%d ", nnn);
	// printf("%d ", mmm);
	if (mm > nn)
	{
		double dy = pos[2*mm+1] - pos[2*nn+1] - round(pos[2*mm+1] - pos[2*nn+1]);
		double dx = pos[2*mm] - pos[2*nn] - round(pos[2*mm] - pos[2*nn]);
		double Dnm = R_eff[nn] + R_eff[mm];
		if (dx * dx + dy * dy < Dnm * Dnm)
		{
			double rymmm = r_shape[mmm] * sin(pos[2*mm+2] + th_shape[mmm]);
			double rynnn = r_shape[nnn] * sin(pos[2*nn+2] + th_shape[nnn]);
			double dyy = (pos[2*mm+1]+rymmm) - (pos[2*nn+1]+rynnn);
			double DDnm = (Dn[nn] + Dn[mm]) / 2.0;
			double rxmmm = r_shape[mmm] * cos(pos[2*mm+2] + th_shape[mmm]);
			double rxnnn = r_shape[nnn] * cos(pos[2*nn+2] + th_shape[nnn]);
			double dxx=(pos[2*mm]+rxmmm) - (pos[2*nn]+rxnnn);
			if (dxx * dxx + dyy * dyy < DDnm * DDnm)
			{
				double ddnm2 = dxx * dxx + dyy * dyy;
                double ddnm = sqrt(ddnm2);
				double DDnm2 = DDnm * DDnm;
                double F = -(DDnm / ddnm - 1.0);
                atomicAdd(&force[2*nn], (double)(F * dxx)); // Fx, nn
                atomicAdd(&force[2*nn+1], (double)(F * dyy)); // Fy, nn
                atomicAdd(&force[2*nn+2], (double)(F * (rxnnn * dyy - rynnn * dxx))); // T, nn
                atomicAdd(&force[2*mm], (double)(-F * dxx)); // Fx, mm
                atomicAdd(&force[2*mm+1], (double)(-F * dyy)); // Fy, mm
                atomicAdd(&force[2*mm+2], (double)(-F * (rxmmm * dyy - rymmm * dxx))); // T, mm
                atomicAdd(&U[nn], (double)(0.5 * F * F));
			}
		}
	}
}

void bumpy_2D_force(int N, int n, double * pos, double * R_eff, double * Dn, double * r_shape, double * th_shape, double * force, double * U)
{
    // step 1: zero all associated arrays: force, U
    zeroArr<<<3*N,1>>>(force);
    zeroArr<<<N,1>>>(U);

    // step 2: contact detection
    bumpy_2D_getForce<<<dim3(N,N,n*n),1>>>(N, n, pos, R_eff, Dn, r_shape, th_shape, force, U);

    // step 3: add up all the U values
    sum_gpu(N, U);
}

void bumpy_2D_comp_VV(int Nc, int n, double * d_pos, double * d_R_eff, double * d_Dn, double * d_r_shape, double * d_th_shape, double * d_U, , double h_dt, double dphi)
{
    int nt = 0;
    double * d_dt;
    cudaMalloc((void**) & d_dt, sizeof(double));
    cudaMemcpy(d_dt, h_dt, sizeof(double), cudaMemcpyHostToDevice);

    // step 1: grow using dphi: Dn, r_shape, m, I
    double * d_phitot;
    cudaMalloc((void**) & d_phitot, Nc*sizeof(double));
    cudaMemcpy(d_phitot, d_m, Nc*sizeof(double), cudaMemcpyHostToDevice);
    sum_gpu(Nc, d_phitot);
    double h_phitot;
    cudaMemcpy(h_phitot, d_phitot, sizeof(double), cudaMemcpyDeviceToHost);
    double h_rsc = sqrt(1.0 + dphi / h_phitot);
    double * d_rsc;
    cudaMalloc((void**) &d_rsc, sizeof(double));
    cudaMemcpy(d_rsc, &h_rsc, sizeof(double), cudaMemcpyHostToDevice);
    grow_cluster_gpu<<<Nc,1>>>(d_Dn, d_m, d_R_eff, d_I, d_rsc);
    grow_bumpy_gpu<<<Nc*n,1>>>(d_r_shape, d_rsc);

    // step 2: get initial U
    double h_U;
    bumpy_2D_force(Nc, n, d_pos, d_R_eff, d_Dn, d_r_shape, d_th_shape, d_force, d_U);
    cudaMemcpy(h_U, d_U, sizeof(double), cudaMemcpyDeviceToHost);

    // step 3: prepare local MD arrays, device level
    double * d_vel;
    double * d_acc_old;
    double * d_acc;
    cudaMalloc((void**) & d_vel, 3*Nc*sizeof(double));
    cudaMalloc((void**) & d_acc_old, 3*Nc*sizeof(double));
    cudaMalloc((void**) & d_acc, 3*Nc*sizeof(double));

    // step 4: minimization
    while (h_U > 0 && nt < Nt)
    {
        // step 1: set CM velocity to 0

        // step 2: update position
        updatePos<<<3*Nc,1>>>(d_pos, d_vel, d_acc_old, d_dt);

        // step 3: get force
        bumpy_2D_force(Nc, n, d_pos, d_R_eff, d_Dn, d_r_shape, d_th_shape, d_force, d_U);

        // step 4: get accel
        updateAcc<<<3*Nc,1>>>(d_acc, d_force, d_mass);

    }

    // step ?: rattler detection and stress calculation
}

void compression_jamming()
{
    // step 1: find initial P

    // step 2: binary search to get to P target

    // step 3: collect and return results
}

*/

int main(int argc, char ** argv)
{
    //Save output file
    string debugfile = "";
    string debugfileext = ".xy";

    // string filename = "E:/Dropbox/Yale/C++/Bumpy/output_"; // Windows
    string filename = "/Users/philipwang/Dropbox/Yale/C++/Bumpy/output_"; // Mac
    // string filename = "/gpfs/loomis/home.grace/ohern/pw374/project/2D_bumpy_CPP_FIRE_N_6_shearG_scan/OUTPUT_TEST/output_"; // Linux
    initializaSavePath(filename, argv);
    debugfile.append(filename);
    debugfile.append(debugfileext);
    char *filechar = new char[filename.length() + 1];
	strcpy(filechar, filename.c_str());
	FILE *out = fopen(filechar, "w+");
    printf("Save char is: %s\n", filename.c_str());
    string loadname = "/Users/philipwang/Dropbox/Yale/Matlab_Code/2D/bumpy_jamming/shear_modulus/input_"; // Mac
    // string loadname = "/gpfs/loomis/home.grace/ohern/pw374/project/2D_bumpy_CPP_FIRE_N_6_shearG_scan/INIT/input_"; // Linux

    initializaLoadPath(loadname, argv);
    printf("Load char is: %s\n", loadname.c_str());

    int Nc;// = 6;
    double mu;// = 1.0;
    int n;// = 15; // n-mers
    double phi_target;// = 1.00;
    int seed;// = 10;
    double G = 1.4;

    Nc = atoi(argv[1]);
    mu = atof(argv[2]);
    n = atoi(argv[3]);
    phi_target = atof(argv[4]);
    seed = atoi(argv[5]);
    double N_temp = Nc * n * G;
    N_temp = N_temp + 0.5 - (N_temp<0);
    int N = int(N_temp);

    mt19937 generator(seed);

    // initialize temp arrays
    int n_l = n * G;
    int n_s = n;
    double * r_shape_l = new double[n_l];
	double * th_shape_l = new double[n_l];
	double * r_shape_s = new double[n_l];
	double * th_shape_s = new double[n_l];

    // initialize particle arrays
    double * h_Dn = new double[Nc];
    double * h_R_eff = new double[Nc];
    double * h_m = new double[Nc];
    double * h_I = new double[Nc]; // moment of inertia of disks along center = 1/2 * m * r * r    
    double * x = new double[Nc];
    double * y = new double[Nc];
    double * th = new double[Nc];
    double * h_r_shape = new double[N];
    double * h_th_shape = new double[N];

    // initialize measurement
    double h_stress[4] = {0}; // Sxx, Sxy, Syx, Syy
    double h_P = 0;
    
    int Ns, Nl;
    double Lx, Ly, dt;

    loadInitialConfiguration(loadname, n_l, n_s, Ns, Nl, Lx, Ly, dt, x, y, th
    , h_Dn, h_R_eff, h_m, h_I, r_shape_l, th_shape_l, r_shape_s, th_shape_s);
    // finally, assign values to r_shape, th_shape
		for (int i=0; i<Nc/2; i++)
		{
			for (int j=0; j<n_l; j++)
			{
				h_r_shape[i*n_l+j] = r_shape_s[j];
				h_th_shape[i*n_l+j] = th_shape_s[j];
			}
		}
        printf("\n");
		for (int i=Nc/2; i<Nc; i++)
		{
			for (int j=0; j<n_l; j++)
			{
				h_r_shape[i*n_l+j] = r_shape_l[j]; // memory overflow?
				h_th_shape[i*n_l+j] = th_shape_l[j]; // memory overlfow? reported by valgrind
			}
		}
    // printArray(Nc, x);
    // printArray(N, h_r_shape);
    // printArray(N, h_th_shape);
    
    delete[] r_shape_l, th_shape_l, r_shape_s, th_shape_s;

    // initialize host arrays
    double * h_pos = new double[3*Nc];
    for (int i=0; i<Nc; i++)
    {
        h_pos[3*i] = x[i];
        h_pos[3*i+1] = y[i];
        h_pos[3*i+2] = th[i];
    }
    // initialize device arrays
    double * d_pos;
    double * d_Dn;
    double * d_R_eff;
    double * d_m;
    double * d_I;
    double * d_r_shape;
    double * d_th_shape;
    // F is arranged as [Fx[0],Fy[0],T[0],...,Fx[Nc-1],Fy[Nc-1],T[Nc-1]]
    double * d_F; // forces per particle, in x, y, th
    int * d_Cn; // contacts per particle
    double * d_U; // U per particle
    double * d_stress;

    // initialize device variables (not array)

    cout << "Good up to here!" << endl;
    // allocate memory to GPU
    cudaMalloc((void**) & d_pos, 3*Nc*sizeof(double));
    cudaMalloc((void**) & d_Dn, Nc*sizeof(double));
    cudaMalloc((void**) & d_R_eff, Nc*sizeof(double));
    cudaMalloc((void**) & d_m, Nc*sizeof(double));
    cudaMalloc((void**) & d_I, Nc*sizeof(double));
    cudaMalloc((void**) & d_r_shape, N*sizeof(double));
    cudaMalloc((void**) & d_th_shape, N*sizeof(double));
    cudaMalloc((void**) & d_F, 3*Nc*sizeof(double));
    cudaMalloc((void**) & d_Cn, Nc*sizeof(int));
    cudaMalloc((void**) & d_U, Nc*sizeof(double));
    // cudaMalloc((void**) & d_stress, 4*sizeof(double));
	// cout << endl;
    // copy initial data to GPU
    cudaMemcpy(d_pos, h_pos, 3*Nc*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Dn, h_Dn, Nc*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_R_eff, h_R_eff, Nc*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_m, h_m, Nc*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_I, h_I, Nc*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r_shape, h_r_shape, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_th_shape, h_th_shape, N*sizeof(double), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_stress, h_stress, 4*sizeof(double), cudaMemcpyHostToDevice);

    
    delete[] h_Dn, h_R_eff, h_m, h_I, x, y, th, h_r_shape, h_th_shape;
    fclose(out);
//	cout << endl;
    return 0;
}

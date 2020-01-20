#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

// srun --pty -p gpu -c 2 -t 2:00:00 --gres=gpu:v100:1 bash
// module restore cuda10 
// nvcc -x cu nested_for_loop.cpp -o test.exe -gencode arch=compute_70,code=sm_70

/*

__global__ void calc(int N, int n, double * pos, double * R_eff, double * Dn, double * r_shape, double * th_shape, double * overlap)
{
	int mm = blockIdx.x;
	int nnn = blockIdx.y;
	int mmm = blockIdx.z;
	for (int nn=0; nn<N; nn++)
	{
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
			double ddnm2 = dxx * dxx + dyy * dyy;
			double DDnm2 = DDnm * DDnm;
			if (dxx * dxx + dyy * dyy < DDnm * DDnm)
			{
				overlap[nn*(N*n*n)+mm*(n*n)+nnn*(n)+mmm] = ddnm2 / DDnm2;
			}
		}
	}
	}
}

*/

__global__ void calc2(int N, int n, double * pos, double * R_eff, double * Dn, double * r_shape, double * th_shape, double * overlap)
{
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
				double DDnm2 = DDnm * DDnm;
				// overlap[nn*(N*n*n)+mm*(n*n)+nnn*(n)+mmm] = ddnm2 / DDnm2;
				// atomicAdd(&del[2*i], (double)(overlap*rij[0]));
				atomicAdd(&overlap[nn], (double) (ddnm2 / DDnm2));
				atomicAdd(&overlap[mm], (double) (ddnm2 / DDnm2));
			}
		}
	}
}

int main(int argc, char ** argv)
{
	//  nvcc -x cu -arch=compute_30 nested_for_loop.cpp -o test.exe
	// ./test.exe 30
    const int N = 3;
    const int n = 21;
	double h_R_eff[] = {0.0557285110124916,0.0557285110124916,0.0557285110124916};
    double h_pos [] = {0.2461858838516268, 0.5, 4.0847976733740090, 0.5744337289871293, 0.1, 6.1226191813248754, 0.0820619612838756,0.7, 3.2923723973920778};//positions of disks, x0, y0, th0..., xN-1, yN-1, thN-1
    double h_th_shape[] = {0.0000000000000000,0.2991993003418849,0.5983986006837699,0.8975979010256550,1.1967972013675399,1.4959965017094248,1.7951958020513099,2.0943951023931948,2.3935944027350797,2.6927937030769646,2.9919930034188496,-2.9919930034188522,-2.6927937030769673,-2.3935944027350819,-2.0943951023931970,-1.7951958020513115,-1.4959965017094261,-1.1967972013675408,-0.8975979010256554,-0.5983986006837702,-0.2991993003418847};//positions of disks
    double h_r_shape[] = {0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613,0.0677762469981613};//positions of disks
	double h_Dn [n] = {0};
	for (int i=0; i<n; i++)
	{
		h_Dn[i] = 0.0142857142857143;
	}
	// rescale
	double scale = atof(argv[1]);
	for (int i=0; i<n; i++)
	{
		h_Dn[i] = h_Dn[i] * scale;
		h_r_shape[i] = h_r_shape[i] * scale;
	}
	for (int i=0; i<N; i++)
	{
		h_R_eff[i] = h_R_eff[i] * scale;
	}
	double overlap [N] = {0};
	double dx, dy, Dnm, dxx, dyy, DDnm, ddnm2, DDnm2;
	double rymmm, rynnn, rxmmm, rxnnn;
	for (int nn=0; nn<N; nn++)
	{
		for (int mm=nn+1; mm<N; mm++)
		{
			Dnm = h_R_eff[nn] + h_R_eff[mm];
			dy = h_pos[2*mm+1] - h_pos[2*nn+1] - round(h_pos[2*mm+1] - h_pos[2*nn+1]);
			if (abs(dy) < Dnm)
			{
				dx = h_pos[2*mm] - h_pos[2*nn] - round(h_pos[2*mm] - h_pos[2*nn]);
				if (dx * dx + dy * dy < Dnm * Dnm)
				{
					for (int nnn=0; nnn<n; nnn++)
					{
						for (int mmm=0; mmm<n; mmm++)
						{
							rymmm = h_r_shape[mmm] * sin(h_pos[2*mm+2] + h_th_shape[mmm]);
							rynnn = h_r_shape[nnn] * sin(h_pos[2*nn+2] + h_th_shape[nnn]);
							dyy=(h_pos[2*mm+1]+rymmm)-(h_pos[2*nn+1]+rynnn);
							DDnm=(h_Dn[nn]+h_Dn[mm])/2.0;
							if (abs(dyy)<Dnm)
							{
								rxmmm = h_r_shape[mmm] * cos(h_pos[2*mm+2] + h_th_shape[mmm]);
								rxnnn = h_r_shape[nnn] * cos(h_pos[2*nn+2] + h_th_shape[nnn]);
								dxx=(h_pos[2*mm]+rxmmm)-(h_pos[2*nn]+rxnnn);
								ddnm2 = dxx * dxx + dyy * dyy;
								DDnm2 = DDnm * DDnm;
								if (dxx * dxx + dyy * dyy < DDnm * DDnm)
								{
									overlap[nn] += ddnm2 / DDnm2;
									overlap[mm] += ddnm2 / DDnm2;
								}
							}
						}
					}

				}
			}
			

		}
	}
	printf("CPU result:\n");
	for (int i=0; i<N; i++)
	{
		//if (overlap[i] > 1e-14) 
		printf("%2.16f ", overlap[i]);

	}
	 
	printf("\n");
	printf("GPU result:\n");

    double * d_pos;
	double * d_R_eff;
	double * d_th_shape;
	double * d_r_shape;
	double * d_Dn;
	double * d_overlap;

	double h_overlap [N] = {0};

	cudaMalloc((void**) & d_pos, 3*N*sizeof(double));
	cudaMalloc((void**) & d_R_eff, N*sizeof(double));
	cudaMalloc((void**) & d_th_shape, n*sizeof(double));
	cudaMalloc((void**) & d_r_shape, n*sizeof(double));
	cudaMalloc((void**) & d_Dn, n*sizeof(double));
	cudaMalloc((void**) & d_overlap, N*sizeof(double));

	cudaMemcpy(d_pos, h_pos, 3*N*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_R_eff, h_R_eff, N*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_th_shape, h_th_shape, n*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r_shape, h_r_shape, n*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Dn, h_Dn, n*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_overlap, h_overlap, N*sizeof(double), cudaMemcpyHostToDevice);

	calc2<<<dim3(N,N,n*n),1>>>(N, n, d_pos, d_R_eff, d_Dn, d_r_shape, d_th_shape, d_overlap);

	cudaMemcpy(h_overlap, d_overlap, N*sizeof(double), cudaMemcpyDeviceToHost);

	for (int i=0; i<N; i++)
	{
		// if (overlap[i] > 1e-14) 
		printf("%2.16f ", h_overlap[i]);

	}
	printf("\n");
	return 0;
}

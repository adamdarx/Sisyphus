#include "defs.h"
#include "init.h"
#include "mks.h"
#include "advance.h"
#include "io.h"

int main()
{
	MKS::init_metricFunc(metricFunc);
	prim = (double (*)[N2M][N3M][NPRIM])(&(a_prim[NG][NG][NG][0]));
	prim_half = (double (*)[N2M][N3M][NPRIM])(&(a_prim_half[NG][NG][NG][0]));
	slope = (double (*)[N2M][N3M][NPRIM])(&(a_slope[NG][NG][NG][0]));
	F1 = (double (*)[N2M][N3M][NPRIM])(&(a_F1[NG][NG][NG][0]));
	F2 = (double (*)[N2M][N3M][NPRIM])(&(a_F2[NG][NG][NG][0]));
	F3 = (double (*)[N2M][N3M][NPRIM])(&(a_F3[NG][NG][NG][0]));
	// Initialize the grid and primitive variables
	init_torus();
	fix_prim(prim);
	bound_prim(prim);
	write_prim(fopen("data0.bin", "wb"));
	double sum = 0;
	for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) sum += prim[i][j][k][RHO];
	printf("Total mass: %lf\n", sum);
	const double tf = 2000;
	int counter = 0;
	double ndt1 = 0, ndt2 = 0, ndt3 = 0;
	while (t < tf) {
		step_ch(ndt1, ndt2, ndt3); // Advance the simulation by one time step)
		char name[32];
		sprintf(name, "data%d.bin", ++counter);
		write_prim(fopen(name, "wb"));
	}
	return 0;
}
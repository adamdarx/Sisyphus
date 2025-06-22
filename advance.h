#pragma once
#include "boundary.h"
#include "defs.h"
#include "fixup.h"
#include "flux.h"
#include "utoprim_2d.h"

void printMatrix(double mat[4][4]) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			printf("%lf ", mat[i][j]);
		}
		printf("\n");
	}
}

double advance(
    double pi[][N2M][N3M][NPRIM],
    double pb[][N2M][N3M][NPRIM],
    double Dt,
    double pf[][N2M][N3M][NPRIM],
    double& ndt1,
    double& ndt2,
    double& ndt3
)
{
    double ndt, U[NPRIM], dU[NPRIM];
    int flag = 0;
    double flag_sum = 0.;
    double gcov[4][4], gcon[4][4], gdet, conn[4][4][4];

    for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) for (int m = 0; m < NPRIM; m++) pf[i][j][k][m] = pi[i][j][k][m];        /* needed for Utoprim */

    ndt1 = fluxcalc(pb, F1, 1);
    ndt2 = fluxcalc(pb, F2, 2);
    ndt3 = fluxcalc(pb, F3, 3);

    // fix_flux(F1, F2, F3);

    flux_ct(F1, F2, F3);

    for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) {
        get_metric(i, j, k, gcov, gcon, gdet);
        get_conn(i, j, k, conn);

        primtosrc(gcov, gcon, gdet, conn, pb[i][j][k], dU);
        primtoconv(gcov, gcon, gdet, pi[i][j][k], U);
		//printf("prim[%d][%d][%d] = %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", i, j, k,
		//	pb[i][j][k][RHO], pb[i][j][k][UU], pb[i][j][k][U1], pb[i][j][k][U2], pb[i][j][k][U3],
		//	pb[i][j][k][B1], pb[i][j][k][B2], pb[i][j][k][B3], pb[i][j][k][BSQ]);
        for (int m = 0; m < NPRIM; m++) {
            U[m] += Dt * (
#if( N1 != 1 )
                - (F1[i + 1][j][k][m] - F1[i][j][k][m]) / dx1
#endif
#if( N2 != 1 )
                - (F2[i][j + 1][k][m] - F2[i][j][k][m]) / dx2
#endif
#if( N3 != 1 )
                - (F3[i][j][k + 1][m] - F3[i][j][k][m]) / dx3
#endif
                + dU[m]
                );
        }
        flag = Utoprim_2d(U, gcov, gcon, gdet, pf[i][j][k]);
        flag_sum += flag;
        pflag[i][j][k] = flag == 0;
    }
    printf("flag_average: %lf\n", flag_sum / (N1 * N2 * N3));
    return(0);
}

void step_ch(double& ndt1, double& ndt2, double& ndt3)
{
    double ndt;
    printf("advancing.\n");
    advance(prim, prim, 0.5 * dt, prim_half, ndt1, ndt2, ndt3);   /* time step primitive variables to the half step */

    printf("fixing.\n");
    fix_prim(prim_half);         /* Set updated densities to floor, set limit for gamma */
    printf("bounding.\n");
    bound_prim(prim_half);    /* Set boundary conditions for primitive variables, flag bad ghost zones */
    printf("u2prim fixing.\n");
    fixup_utoprim(prim_half);  /* Fix the failure points using interpolation and updated ghost zone values */

    //ZLOOP eVol(p,p,prim_half,0.5*dt,i,j,k,0,1); /* Account for floor heat */
    printf("bounding.\n");
    bound_prim(prim_half);    /* Reset boundary conditions with fixed up points */



    /* Repeat and rinse for the full time (aka corrector) step:  */
    printf("advancing.\n");
    advance(prim, prim_half, dt, prim, ndt1, ndt2, ndt3);

    printf("fixing.\n");
    fix_prim(prim);
    bound_prim(prim);
    fixup_utoprim(prim);

    bound_prim(prim);


    /* increment time */
    ndt = 1. / (1. / (ndt1) + 1. / (ndt2) + 1. / (ndt3));
    if (ndt > SAFE * dt) ndt = SAFE * dt;
    dt = ndt;
    if (t + dt > tf) dt = tf - t;
    t += dt;
	printf("dt: %lf, t: %lf\n", dt, t);
}


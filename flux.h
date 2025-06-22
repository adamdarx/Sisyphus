#pragma once
#include "defs.h"
#include "interp.h"
#include "metric.h"
#include "phys.h"

/** algorithmic choices **/

/* use local lax-friedrichs or HLL flux:  these are relative weights on each numerical flux */

#define HLLF  (0.0)
#define LAXF  (1.0)

/** end algorithmic choices **/

inline double max(double a, double b) {
	return (a > b) ? a : b;
}

inline double min(double a, double b) {
	return (a < b) ? a : b;
}

/***********************************************************************************************/
/***********************************************************************************************
 fluxcalc():
 ---------
 -- sets the numerical fluxes, avaluated at the cell boundaries using the slope limiter
 slope_lim();

 -- only has HLL and Lax-Friedrichs  approximate Riemann solvers implemented;

 ***********************************************************************************************/

double fluxcalc(
    double pr[][N2M][N3M][NPRIM],
    double F[][N2M][N3M][NPRIM],
    int dir
)
{
    double p_l[NPRIM], p_r[NPRIM], F_l[NPRIM], F_r[NPRIM], U_l[NPRIM], U_r[NPRIM];
    double cmax_l, cmax_r, cmin_l, cmin_r, cmax, cmin, ndt, dtij;
    double ctop, s_l, s_r;
    double cmax_lw, cmax_rw, cmin_lw, cmin_rw, cmaxw, cminw;
    double ctopw, dummy;
    double tmp;
    double bsq;
	for (int i = -1; i <= N1; i++) for (int j = -1; j <= N2; j++) for (int k = -1; k <= N3; k++) for (int m = 0; m < NPRIM; m++) {
        slope[i][j][k][m] = 0.;
        F[i][j][k][m] = 0.;
	}

    int idel = dir == 1 ? 1 : 0;
    int jdel = dir == 2 ? 1 : 0;
    int kdel = dir == 3 ? 1 : 0;
    int ihalfdel = dir == 1 ? 0.5 : 0;
    int jhalfdel = dir == 2 ? 0.5 : 0;
    int khalfdel = dir == 3 ? 0.5 : 0;

    /* then evaluate slopes on active grid plus 1-cell-wide layer of ghost cells */
    for (int i = -1; i <= N1; i++) for (int j = -1; j <= N2; j++) for (int k = -1; k <= N3; k++) for (int m = 0; m < NPRIM; m++) {
          slope[i][j][k][m] = slope_lim(
                                  pr[i - idel][j - jdel][k - kdel][m],
                                  pr[i][j][k][m],
                                  pr[i + idel][j + jdel][k + kdel][m]
                                  );
    }

    ndt = 1.e9;

    for (int i = -1 * !idel; i <= N1; i++) for (int j = -1 * !jdel; j <= N2; j++) for (int k = -1 * !kdel; k <= N3; k++) {
        if (dir == 2 && (j == 0 || j == N2)) {
            //FIXFLUXMARK: fix_flux-like zero out of fluxes in x2-dir
            for (int m = 0; m < NPRIM; m++) F[i][j][k][m] = 0.;
        }
        else {
            for (int m = 0; m < NPRIM; m++) {
                p_l[m] = pr[i - idel][j - jdel][k - kdel][m] + 0.5 * slope[i - idel][j - jdel][k - kdel][m];
                p_r[m] = pr[i][j][k][m] - 0.5 * slope[i][j][k][m];
            }

            double gcov_l[4][4], gcon_l[4][4], gdet_l;
            double gcov_r[4][4], gcon_r[4][4], gdet_r;

            get_metric(i - ihalfdel, j - jhalfdel, k - khalfdel, gcov_l, gcon_l, gdet_l);
            get_metric(i + ihalfdel, j + jhalfdel, k + khalfdel, gcov_r, gcon_r, gdet_r);

            primtoflux(gcov_l, gcon_l, gdet_l, dir, p_l, F_l);
            primtoflux(gcov_r, gcon_r, gdet_r, dir, p_r, F_r);

            primtoflux(gcov_l, gcon_l, gdet_l, dir, p_l, U_l);
            primtoflux(gcov_r, gcon_r, gdet_r, dir, p_r, U_r);

            velocity(gcov_l, gcon_l, gdet_l, dir, p_l, cmax_l, cmin_l);
            velocity(gcov_r, gcon_r, gdet_r, dir, p_r, cmax_r, cmin_r);

            cmax = fabs(max(max(0., cmax_l), cmax_r));
            cmin = fabs(max(max(0., -cmin_l), -cmin_r));
            ctop = max(cmax, cmin);

            for (int m = 0; m < NPRIM; m++) {
                F[i][j][k][m] =
                    HLLF * (
                        (cmax * F_l[m] + cmin * F_r[m]
                    - cmax * cmin * (U_r[m] - U_l[m])) /
                        (cmax + cmin + SMALL)
                        ) +
                    LAXF * (
                        0.5 * (F_l[m] + F_r[m]
                        - ctop * (U_r[m] - U_l[m]))
                        );
            }

            switch (dir)
            {
            case 1: dtij = cour * dx1 / ctop; break;
            case 2: dtij = cour * dx2 / ctop; break;
            case 3: dtij = cour * dx3 / ctop; break;
            default:break;
            }
            if (dtij < ndt) ndt = dtij;

        }
    }
    return(ndt);

}

/***********************************************************************************************/
/***********************************************************************************************
 flux_ct():
 ---------
 -- performs the flux-averaging used to preserve the del.B = 0 constraint (see Toth 2000);

 ***********************************************************************************************/
void flux_ct(double F1[][N2M][N3M][NPRIM], double F2[][N2M][N3M][NPRIM], double F3[][N2M][N3M][NPRIM])
{
    /* calculate EMFs */
#define DOE1 (N2>1 && N3>1)
#define DOE2 (N1>1 && N3>1)
#define DOE3 (N1>1 && N2>1)

    for(int i = 0; i <= N1; i++) for (int j = 0; j <= N2; j++) for (int k = 0; k <= N3; k++) {
        //E_i = e_{ijk} F_j B_k
#if(DOE1)
    //E1
        emf[1][i][j][k] = 0.25 * ((F2[i][j][k][B3] + F2[i][j][k - 1][B3])
            - (F3[i][j][k][B2] + F3[i][j - 1][k][B2])
            );
#endif

#if(DOE2)
        //E2
        emf[2][i][j][k] = 0.25 * ((F3[i][j][k][B1] + F3[i - 1][j][k][B1])
            - (F1[i][j][k][B3] + F1[i][j][k - 1][B3])
            );
#endif

#if(DOE3)
        //E3
        emf[3][i][j][k] = 0.25 * ((F1[i][j][k][B2] + F1[i][j - 1][k][B2])
            - (F2[i][j][k][B1] + F2[i - 1][j][k][B1])
            );
#endif
    }
    /* rewrite EMFs as fluxes, after Toth */
    for (int i = 0; i <= N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) {
#if(N1>1)
        F1[i][j][k][B1] = 0.;
#endif
#if(DOE3)
        F1[i][j][k][B2] = 0.5 * (emf[3][i][j][k] + emf[3][i][j + 1][k]);
#endif
#if(DOE2)
        F1[i][j][k][B3] = -0.5 * (emf[2][i][j][k] + emf[2][i][j][k + 1]);
#endif
    }

    for (int i = 0; i < N1; i++) for (int j = 0; j <= N2; j++) for (int k = 0; k < N3; k++) {
#if(DOE3)
        F2[i][j][k][B1] = -0.5 * (emf[3][i][j][k] + emf[3][i + 1][j][k]);
#endif
#if(DOE1)
        F2[i][j][k][B3] = 0.5 * (emf[1][i][j][k] + emf[1][i][j][k + 1]);
#endif
#if(N2>1)
        F2[i][j][k][B2] = 0.;
#endif
    }

    for (int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k <= N3; k++) {
#if(DOE2)
        F3[i][j][k][B1] = 0.5 * (emf[2][i][j][k] + emf[2][i + 1][j][k]);
#endif
#if(DOE1)
        F3[i][j][k][B2] = -0.5 * (emf[1][i][j][k] + emf[1][i][j + 1][k]);
#endif
#if(N3>1)
        F3[i][j][k][B3] = 0.;
#endif
    }
}


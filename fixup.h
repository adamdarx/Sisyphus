#pragma once
#include "defs.h"

//compute bsq
double _bsq_cal(int i, int j, int k) {
    double b[4], bsq;  // b^mu and bsq
    b[0] = 0.;
    for (int m = 1; m < 4; m++) {
        for (int n = 0; n < 4; n++) {
            b[0] += gdd_mks[i][j][k][m][n] * prim[i][j][k][B1 + m - 1] * prim[i][j][k][U0 + n - 1];
        }
    }
    for (int m = 1; m < 4; m++) {
        b[m] = (prim[i][j][k][B1 + m - 1] + b[0] * prim[i][j][k][U1 + m - 1]) / (SMALL + prim[i][j][k][U0]);
    }
    bsq = 0.;
    for (int m = 0; m < 4; m++) {
        for (int n = 0; n < 4; n++) {
            bsq += gdd_mks[i][j][k][m][n] * b[m] * b[n];
        }
    }
    bsq = fabs(bsq);
    return bsq;
}

//fix primitive variable by adding density rho
void fix_prim(double prim[N1M][N2M][N3M][NPRIM]) {
    double r, rho_floor, ug_floor, bsq, sigma;
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                r = BL_coord1[i][j][k];
                rho_floor = RHOMIN * pow(r, -3. / 2.);
                ug_floor = UUMIN * pow(r, -3. / 2. * gam);
                if (prim[i][j][k][RHO] < rho_floor) prim[i][j][k][RHO] = rho_floor;
                if (prim[i][j][k][UU] < ug_floor) prim[i][j][k][UU] = ug_floor;
                bsq = _bsq_cal(i, j, k);
                sigma = bsq / prim[i][j][k][RHO];
                if (sigma > SIGMAMAX) prim[i][j][k][RHO] = bsq / SIGMAMAX;
            }
        }
    }
}

/**************************************************************************************
 INTERPOLATION STENCILS:
 ------------------------
   -- let the stencils be characterized by the following numbering convention:

           1 2 3
           8 x 4      where x is the point at which we are interpolating
           7 6 5
*******************************************************************************************/

//!!!ATCH: leave the fixup stencils unchanged, i.e., never couple cells at different k for fixups (for now)

/* 12345678 */
#define AVG8(pr,i,j,k,m)  \
        (0.125*(pr[i-1][j+1][k][m]+pr[i][j+1][k][m]+pr[i+1][j+1][k][m]+pr[i+1][j][k][m]+pr[i+1][j-1][k][m]+pr[i][j-1][k][m]+pr[i-1][j-1][k][m]+pr[i-1][j][k][m]))

/* 2468  */
#define AVG4_1(pr,i,j,k,m) (0.25*(pr[i][j+1][k][m]+pr[i][j-1][k][m]+pr[i-1][j][k][m]+pr[i+1][j][k][m]))

/* 1357  */
#define AVG4_2(pr,i,j,k,m) (0.25*(pr[i+1][j+1][k][m]+pr[i+1][j-1][k][m]+pr[i-1][j+1][k][m]+pr[i-1][j-1][k][m]))

/* + shaped,  Linear interpolation in X1 or X2 directions using only neighbors in these direction */
/* 48  */
#define AVG2_X1(pr,i,j,k,m) (0.5*(pr[i-1][j  ][k][m]+pr[i+1][j  ][k][m]))
/* 26  */
#define AVG2_X2(pr,i,j,k,m) (0.5*(pr[i  ][j-1][k][m]+pr[i  ][j+1][k][m]))

/* x shaped,  Linear interpolation diagonally along both X1 and X2 directions "corner" neighbors */
/* 37  */
#define AVG2_1_X1X2(pr,i,j,k,m) (0.5*(pr[i-1][j-1][k][m]+pr[i+1][j+1][k][m]))
/* 15  */
#define AVG2_2_X1X2(pr,i,j,k,m) (0.5*(pr[i-1][j+1][k][m]+pr[i+1][j-1][k][m]))

void fixup_utoprim(double prim[N1M][N2M][N3M][NPRIM])
{
    int i, j, k, m;
    static int pf[9];

    /* Fix the interior points first */
    for(int i = 0; i < N1; i++) for (int j = 0; j < N2; j++) for (int k = 0; k < N3; k++) {
        if (pflag[i][j][k] == 0) {
            pf[1] = pflag[i - 1][j + 1][k];   pf[2] = pflag[i][j + 1][k];  pf[3] = pflag[i + 1][j + 1][k];
            pf[8] = pflag[i - 1][j][k];                                    pf[4] = pflag[i + 1][j][k];
            pf[7] = pflag[i - 1][j - 1][k];   pf[6] = pflag[i][j - 1][k];  pf[5] = pflag[i + 1][j - 1][k];

            /* Now the pf's  are true if they represent good points */

      //      if(      pf[1]&&pf[2]&&pf[3]&&pf[4]&&pf[5]&&pf[6]&&pf[7]&&pf[8] ){ FLOOP prim[i][j][k] = AVG8(            prim,i,j,k)                   ; }
      //      else if(        pf[2]&&       pf[4]&&       pf[6]&&       pf[8] ){ FLOOP prim[i][j][k] = AVG4_1(          prim,i,j,k)                   ; }
      //      else if( pf[1]&&       pf[3]&&       pf[5]&&       pf[7]        ){ FLOOP prim[i][j][k] = AVG4_2(          prim,i,j,k)                   ; }
      //      else if(               pf[3]&&pf[4]&&              pf[7]&&pf[8] ){ FLOOP prim[i][j][k] = 0.5*(AVG2_1_X1X2(prim,i,j,k)+AVG2_X1(prim,i,j,k)); }
      //      else if(        pf[2]&&pf[3]&&              pf[6]&&pf[7]        ){ FLOOP prim[i][j][k] = 0.5*(AVG2_1_X1X2(prim,i,j,k)+AVG2_X2(prim,i,j,k)); }
      //      else if( pf[1]&&              pf[4]&&pf[5]&&              pf[8] ){ FLOOP prim[i][j][k] = 0.5*(AVG2_2_X1X2(prim,i,j,k)+AVG2_X1(prim,i,j,k)); }
      //      else if( pf[1]&&pf[2]&&              pf[5]&&pf[6]               ){ FLOOP prim[i][j][k] = 0.5*(AVG2_2_X1X2(prim,i,j,k)+AVG2_X2(prim,i,j,k)); }
      //      else if(               pf[3]&&                     pf[7]        ){ FLOOP prim[i][j][k] = AVG2_1_X1X2(     prim,i,j,k)                   ; }
      //      else if( pf[1]&&                     pf[5]                      ){ FLOOP prim[i][j][k] = AVG2_2_X1X2(     prim,i,j,k)                   ; }
      //      else if(        pf[2]&&                     pf[6]               ){ FLOOP prim[i][j][k] = AVG2_X2(         prim,i,j,k)                   ; }
      //      else if(                      pf[4]&&                     pf[8] ){ FLOOP prim[i][j][k] = AVG2_X1(         prim,i,j,k)                   ; }

      // Old way:
            if (pf[2] && pf[4] && pf[6] && pf[8]) {
                for (int m = 0; m < NPRIM; m++) prim[i][j][k][m] = AVG4_1(prim, i, j, k, m);
            }
            else if (pf[1] && pf[3] && pf[5] && pf[7]) {
                for (int m = 0; m < NPRIM; m++) prim[i][j][k][m] = AVG4_2(prim, i, j, k, m);
            }
            else {
                for (m = RHO; m <= UU; m++) {
                    prim[i][j][k][m] = 0.5 * (AVG4_1(prim, i, j, k, m) + AVG4_2(prim, i, j, k, m));
                }
                prim[i][j][k][U1] = prim[i][j][k][U2] = prim[i][j][k][U3] = 0.;
            }

            pflag[i][j][k] = 1;                  /* The cell has been fixed so we can use it for interpolation elsewhere */
            fix_prim(prim);                      /* Floor and limit gamma the interpolated value */
        }
    }

    return;
}

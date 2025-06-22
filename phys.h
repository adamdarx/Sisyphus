#pragma once
/*
In this file we define the functions that are used to transform
primitive variables to conserved variables. For example, U, F, S and so on.
*/

#include <cmath>
#include "la.h"

/*
所有的变换必须基于单一格点
*/

void get_u_and_b(double gcov[4][4], double gcon[4][4], double gdet,
    double pr[NPRIM],
    double ucon[4], double ucov[4], double bcon[4], double bcov[4])
{
    double alpha, bsq, Bv, gamma[3][3], Gamma, Gvsq;
    alpha = 1.0 / sqrt(-gcon[0][0]);
    for (int row = 1; row < 4; ++row) {
        for (int col = 1; col < 4; ++col) {
            gamma[row - 1][col - 1] = gcov[row][col];
        }
    }
    double Gv[3] = { pr[U1], pr[U2], pr[U3] };
    Gvsq = square(Gv, gamma);
    Gamma = sqrt(1 + Gvsq);

    // calculate u and b
    // first ucon
    ucon[0] = Gamma / alpha;
    ucon[1] = pr[U1] - ucon[0] * gcon[0][1] * alpha * alpha;
    ucon[2] = pr[U2] - ucon[0] * gcon[0][2] * alpha * alpha;
    ucon[3] = pr[U3] - ucon[0] * gcon[0][3] * alpha * alpha;
    // then ucov
    for (int row = 0; row < 4; ++row)
    {
        ucov[row] = 0;
        for (int col = 0; col < 4; ++col)
            ucov[row] += gcov[row][col] * ucon[col];
    }
    // now bcon
    bcon[0] = pr[B1] * ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3];
    bcon[1] = (pr[B1] + bcon[0] * ucon[1]) / ucon[0];
    bcon[2] = (pr[B2] + bcon[0] * ucon[2]) / ucon[0];
    bcon[3] = (pr[B3] + bcon[0] * ucon[3]) / ucon[0];
    // finally bcov
    for (int row = 0; row < 4; ++row)
    {
        bcov[row] = 0;
        for (int col = 0; col < 4; ++col)
            bcov[row] += gcov[row][col] * bcon[col];
    }
    return;
}



void get_mhd(double gcov[4][4], double gcon[4][4], double gdet,
    int dir, double pr[NPRIM], double mhd[4])
{
    const double gam = (5. / 3.);
    double ucon[4], ucov[4], bcon[4], bcov[4];
    get_u_and_b(gcov, gcon, gdet, pr, ucon, ucov, bcon, bcov);

    double rho = pr[RHO];    // density
    double u = pr[UU];      // internal energy
    double P = (gam - 1.) * u;
    double w = P + rho + u;
    double bsq = dot(bcon, bcov);
    double eta = w + bsq;
    double ptot = P + 0.5 * bsq;
    for (int index = 0; index < 4; ++index)
        mhd[index] = eta * ucon[dir] * ucov[index] + ptot * (dir == index ? 1 : 0) - bcon[dir] * bcov[index];
}



void primtoflux(double gcov[4][4], double gcon[4][4], double gdet,
    int dir, double pr[NPRIM], double flux[NPRIM])
{
    const double gam = (5. / 3.);
    double ucon[4], ucov[4], bcon[4], bcov[4];
    double mhd[4];
    get_u_and_b(gcov, gcon, gdet, pr, ucon, ucov, bcon, bcov);
    get_mhd(gcov, gcon, gdet, dir, pr, mhd);
    for (int m = 0; m < NPRIM; m++) flux[m] = 0.;
    flux[RHO] = pr[RHO] * ucon[dir];

    flux[UU] = mhd[0] + flux[0];
    flux[U1] = mhd[1];
    flux[U2] = mhd[2];
    flux[U3] = mhd[3];

    flux[B1] = bcon[1] * ucon[dir] - bcon[dir] * ucon[1];
    flux[B2] = bcon[2] * ucon[dir] - bcon[dir] * ucon[2];
    flux[B3] = bcon[3] * ucon[dir] - bcon[dir] * ucon[3];
    for (int index = 0; index < NPRIM; ++index)
        flux[index] *= sqrt(fabs(gdet));
    return;
}



void primtoconv(double gcov[4][4], double gcon[4][4], double gdet,
    double pr[NPRIM], double conv[NPRIM])
{
    primtoflux(gcov, gcon, gdet, 0, pr, conv);
    return;
}



void primtosrc(
    double gcov[4][4], double gcon[4][4], double gdet, double conn[4][4][4],
    double pr[NPRIM],
    double src[NPRIM])
{
    const double gam = (5. / 3.);
    double mhds[4][4];
    for (int dir = 0; dir < 4; ++dir) {
        get_mhd(gcov, gcon, gdet,
            dir, pr,
            mhds[dir]);
    }

    for (int index = 0; index < NPRIM; ++index)
        src[index] = 0.;

    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col) {
            src[UU] += mhds[row][col] * conn[col][0][row];
            src[U1] += mhds[row][col] * conn[col][1][row];
            src[U2] += mhds[row][col] * conn[col][2][row];
            src[U3] += mhds[row][col] * conn[col][3][row];
        }

    for (int index = 0; index < NPRIM; ++index)
        src[index] *= sqrt(fabs(gdet));
}



void velocity(
    double gcov[4][4], double gcon[4][4], double gdet,
    int dir,
    double pr[NPRIM],
    double& vmax, double& vmin)
{
    double discr, vp, vm, bsq, EE, EF, va2, cs2, cms2, rho, u;
    double Acov[4], Bcov[4], Acon[4], Bcon[4];
    double Asq, Bsq, Au, Bu, AB, Au2, Bu2, AuBu, A, B, C;
    double gam = (5. / 3.);

    for (int index = 0; index < 4; ++index) Acov[index] = 0.;
    Acov[dir] = 1.;
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            Acon[row] += gcon[row][col] * Acov[col];

    for (int index = 0; index < 4; ++index) Bcov[index] = 0.;
    Bcov[0] = 1.;
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            Bcon[row] += gcon[row][col] * Bcov[col];

    double ucon[4], ucov[4], bcon[4], bcov[4];
    get_u_and_b(gcov, gcon, gdet,
        pr,
        ucon, ucov, bcon, bcov);

    bsq = dot(bcon, bcov);
    rho = pr[RHO];
    u = pr[UU];
    EF = rho + gam * u;
    EE = bsq + EF;
    va2 = bsq / EE;
    cs2 = gam * (gam - 1.) * u / EF;
    cms2 = cs2 + va2 - cs2 * va2;

    if (cms2 < 0.) {
        cms2 = 1e-10;
    }
    if (cms2 > 1.) {
        cms2 = 1.;
    }

    Asq = dot(Acon, Acov);
    Bsq = dot(Bcon, Bcov);
    Au = dot(Acov, ucon);
    Bu = dot(Bcov, ucon);
    AB = dot(Acon, Bcov);
    Au2 = Au * Au;
    Bu2 = Bu * Bu;
    AuBu = Au * Bu;

    A = Bu2 - (Bsq + Bu2) * cms2;
    B = 2. * (AuBu - (AB + AuBu) * cms2);
    C = Au2 - (Asq + Au2) * cms2;

    discr = B * B - 4. * A * C;
    if ((discr < 0.0) && (discr > -1.e-10)) discr = 0.0;
    else if (discr < -1.e-10) discr = 0.;

    discr = sqrt(discr);
    vp = -(-B + discr) / (2. * A);
    vm = -(-B - discr) / (2. * A);

    if (vp > vm) {
        vmax = vp;
        vmin = vm;
    }
    else {
        vmax = vm;
        vmin = vp;
    }
}

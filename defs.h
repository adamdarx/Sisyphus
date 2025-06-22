#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
constexpr auto a = (0.9375);
constexpr auto h = (0.);
constexpr auto SMALL = (1.e-5);
constexpr auto theta = 0.3;
constexpr auto NDIM = (4);
constexpr auto N1 = (64);
constexpr auto N2 = (64);
constexpr auto N3 = (8);
constexpr auto NG = (2);
constexpr auto N1M = N1 + 2 * NG;
constexpr auto N2M = N2 + 2 * NG;
constexpr auto N3M = N3 + 2 * NG;
constexpr auto PI = (3.14159265358979323846);
constexpr auto X1min = (0.19325057145871735);
constexpr auto X1max = (7.824046010856292);
constexpr auto X2min = (1e-2);
constexpr auto X2max = (PI - 1e-2);
constexpr auto X3min = (1e-2);
constexpr auto X3max = (2. * PI - 1e-2);
constexpr auto R0 = (0.);
constexpr auto cour = 0.9;
constexpr auto LEFT = 0;
constexpr auto RIGHT = 1;
constexpr auto NEG = 0;
constexpr auto POS = 1;
constexpr auto RHOMIN = (1.e-6);
constexpr auto UUMIN = (1.e-8);
constexpr auto SIGMAMAX = (50.);

//FM_torus disk parameter
constexpr auto rin = (6.);
constexpr auto rmax = (12.);
constexpr auto beta = (100.);
constexpr auto gam = (5. / 3.);
constexpr auto kappa = (1.e-3);

unsigned short max_iter = 10;		// maximum of iteration
double tol = 1e-10;					// tolerance of root devation
auto physicalTimeMax = 1000;		// maximum of physical time
double frame = 0.1;
//MKS grid
double Xgrid1[N1][N2][N3];
double Xgrid2[N1][N2][N3];
double Xgrid3[N1][N2][N3];
//MKS grid spacing
double dx1 = (X1max - X1min) / (N1);
double dx2 = (X2max - X2min) / (N2);
double dx3 = (X3max - X3min) / (N3);
//KS grid
double KS_coord1[N1][N2][N3];
double KS_coord2[N1][N2][N3];
double KS_coord3[N1][N2][N3];
//BL grid
double BL_coord1[N1][N2][N3];
double BL_coord2[N1][N2][N3];
double BL_coord3[N1][N2][N3];

//metric at grid point
double gdd_bl[N1][N2][N3][NDIM][NDIM];
double guu_bl[N1][N2][N3][NDIM][NDIM];
double gdet_bl[N1][N2][N3];              /*sqrt(-g_bl)*/
double gdd_ks[N1][N2][N3][NDIM][NDIM];
double guu_ks[N1][N2][N3][NDIM][NDIM];
double gdet_ks[N1][N2][N3];              /*sqrt(-g_ks)*/
double gdd_mks[N1][N2][N3][NDIM][NDIM];
double guu_mks[N1][N2][N3][NDIM][NDIM];
double gdet_mks[N1][N2][N3];             /*sqrt(-g_mks)*/

//Jacobian matrix at grid point
double J_bl2ks[N1][N2][N3][NDIM][NDIM];
double J_ks2bl[N1][N2][N3][NDIM][NDIM];
double J_ks2mks[N1][N2][N3][NDIM][NDIM];
double J_mks2ks[N1][N2][N3][NDIM][NDIM];

//primitive variables
constexpr auto NPRIM = (10);
constexpr auto NPR = (9);
constexpr auto RHO = (0);
constexpr auto UU = (1);
constexpr auto U0 = (2);
constexpr auto U1 = (3);
constexpr auto U2 = (4);
constexpr auto U3 = (5);
constexpr auto B1 = (6);
constexpr auto B2 = (7);
constexpr auto B3 = (8);
constexpr auto BSQ = (9);
double A[N1 + 1][N2 + 1][N3 + 1];
double a_prim[N1M][N2M][N3M][NPRIM];
double a_prim_half[N1M][N2M][N3M][NPRIM];
double a_slope[N1M][N2M][N3M][NPRIM];
double a_F1[N1M][N2M][N3M][NPRIM];
double a_F2[N1M][N2M][N3M][NPRIM];
double a_F3[N1M][N2M][N3M][NPRIM];
double (*prim)[N2M][N3M][NPRIM];
double (*prim_half)[N2M][N3M][NPRIM];
double (*slope)[N2M][N3M][NPRIM];
double (*F1)[N2M][N3M][NPRIM];
double (*F2)[N2M][N3M][NPRIM];
double (*F3)[N2M][N3M][NPRIM];
double emf[NDIM][N1 + 1][N2 + 1][N3 + 1];
int pflag[N1][N2][N3]; // pflag[i][j][k] = 0 means the cell is good for utoprim, otherwise it is bad
typedef double(*MetricComponent)(double, double, double, double);
MetricComponent metricFunc[4][4];
double t = 1e-5, dt = 1e-5; // physical time and time step
const double tf = 2000;
const double SAFE = (1.3);
enum limiter_type { MC, VANL, MINM };
static enum limiter_type lim = MC; // default slope limiter is Woodward's Monotonized Central (MC) limiter

// functions

//advance.h
double advance(double[][N2M][N3M][NPRIM], double[][N2M][N3M][NPRIM], double, double[][N2M][N3M][NPRIM], double&, double&, double&);
void step_ch(double&, double&, double&);
//boundary.h
void bound_prim(double [N1M][N2M][N3M][NPRIM]);
//fixup.h
void fix_prim(double[N1M][N2M][N3M][NPRIM]);
void fixup_utoprim(double[N1M][N2M][N3M][NPRIM]);
//flux.h
double fluxcalc(double[][N2M][N3M][NPRIM], double[][N2M][N3M][NPRIM], int);
void flux_ct(double[][N2M][N3M][NPRIM], double[][N2M][N3M][NPRIM], double[][N2M][N3M][NPRIM]);
//init.h
void init_torus();
//interp.h
double slope_lim(double, double, double);
//io.h
int write_prim(FILE*);
//la.h
double dot(double a[4], double b[4]);
double dot(double rowVec[3], double colVec[3], double metric[3][3]);
double square(double vec[3], double metric[3][3]);
double contract(double A[3][3], double B[3][3]);
double determinant(double matrix[4][4]);
int inverseMatrix(double[4][4], double[4][4], int);
//metric.h
void get_metric(int i, int j, int k, double gcov[4][4], double gcon[4][4], double &gdet);
void get_conn(int i, int j, int k, double conn[4][4][4]);
//mks.h
/////////////////////////////
/* Check mks.h for details */
/////////////////////////////
//phys.h
void get_u_and_b(double gcov[4][4], double gcon[4][4], double gdet,
    double pr[NPRIM],
    double ucon[4], double ucov[4], double bcon[4], double bcov[4]);
void get_mhd(double gcov[4][4], double gcon[4][4], double gdet,
    int dir, double pr[NPRIM], double mhd[4]);
void primtoflux(double gcov[4][4], double gcon[4][4], double gdet,
    int dir, double pr[NPRIM], double flux[NPRIM]);
void primtoconv(double gcov[4][4], double gcon[4][4], double gdet,
    double pr[NPRIM], double conv[NPRIM]);
void primtosrc(
    double gcov[4][4], double gcon[4][4], double gdet, double conn[4][4][4],
    double pr[NPRIM],
    double src[NPRIM]);
void velocity(
    double gcov[4][4], double gcon[4][4], double gdet,
    int dir,
    double pr[NPRIM],
    double& vmax, double& vmin);
//ranc.h
double ranc(int iseed);
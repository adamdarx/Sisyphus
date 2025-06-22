#pragma once
#include "u2p_defs.h"


void primtoU_g( double prim[], double gcov[][4], double gcon[][4], double gdet, 
		double U[] );
void ucon_calc_g(double prim[],double gcov[][4],double gcon[][4],double ucon[]);
void raise_g(double vcov[], double gcon[][4], double vcon[]);
void lower_g(double vcon[], double gcov[][4], double vcov[]);
void ncov_calc(double gcon[][4],double ncov[]) ;
void bcon_calc_g(double prim[],double ucon[],double ucov[],double ncov[],double bcon[]); 
double pressure_rho0_u(double rho0, double u);
double pressure_rho0_w(double rho0, double w);


/**********************************************************************/
/******************************************************************
   primtoU_g(): 

       -- calculates the conserved variables from the primitive variables 
            and the metric;
       -- assumes that the conserved and primitive variables are defined ala HARM:

              /  rho u^t           \
         U =  |  T^t_\mu + rho u^t |  sqrt(-det(g_{\mu\nu}))
              \   B^i              /

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /

**************************************************************************/
void primtoU_g(
	       double prim[8],       /* primitive variables */
	       double gcov[4][4],    /* covariant (index dn) form of metric */
	       double gcon[4][4],    /* contravariant (index up) form of metric */
	       double gdet,                /* det(g_{\mu \nu}) */
	       double U[8]           /* matrix of derivatives */
	       ) {
  int i,j ;
  double rho0 ;
  static double ucon[4],ucov[4],bcon[4],bcov[4],ncov[4] ;
  double gamma,n_dot_b,bsq,u,p,w, alpha ;

	
  /* Calculate auxiliary quantities: */
  alpha = 1.0/sqrt(-gcon[0][0]);

  ucon_calc_g(prim,gcov,gcon,ucon) ;
  lower_g(ucon,gcov,ucov) ;
  ncov_calc(gcon,ncov) ;

  gamma = -ncov[0]*ucon[0] ;

  bcon_calc_g(prim,ucon,ucov,ncov,bcon) ;
  lower_g(bcon,gcov,bcov) ;

  n_dot_b = 0. ;
  for(i=0;i<4;i++) n_dot_b += ncov[i]*bcon[i] ;
  bsq = 0. ;
  for(i=0;i<4;i++) bsq += bcov[i]*bcon[i] ;

  rho0 = prim[RHO] ;
  u = prim[UU] ;
  p = pressure_rho0_u(rho0,u) ;
  w = rho0 + u + p ;

  // Now set the conserved variables themselves, using HARM's definition:
  U[RHO] = ucon[0]*rho0 ;

  for( i = 0; i < 4; i++) {
    U[QCOV0+i] = gamma*(w + bsq)*ucov[i] 
      - (p + bsq/2.)*ncov[i] 
      + n_dot_b*bcov[i] ;

    U[QCOV0+i] /= alpha;
  }

  U[QCOV0] = U[QCOV0] + U[RHO];
  U[BCON1] = prim[BCON1] ;
  U[BCON2] = prim[BCON2] ;
  U[BCON3] = prim[BCON3] ;

  for(i = 0; i < 8; i++ ) {
    U[i] *= sqrt(fabs(gdet));
  }

  return ;
}

/**********************************************************************/
/******************************************************************

  ucon_calc_g(): 
    
       -- calculates the contravariant (up) components of the four-velocity
          given the primitive variables, of which the velocity is 
          \tilde{u}^i = \gamma v^j  where v^j is the velocity of the 
          flow w.r.t a normal observer to the coordinates;

       -- also requires the metric and inverse metric;

       -- assumes:

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /

******************************************************************/
void ucon_calc_g(double prim[8],double gcov[4][4],double gcon[4][4],
		 double ucon[4])
{
  double u_tilde_con[4] ;
  double u_tilde_sq ;
  double gamma,lapse ;
  int i,j ;
	
  u_tilde_con[0] = 0. ;
  u_tilde_con[1] = prim[UTCON1] ;
  u_tilde_con[2] = prim[UTCON2] ;
  u_tilde_con[3] = prim[UTCON3] ;

  u_tilde_sq = 0. ;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      u_tilde_sq += gcov[i][j]*u_tilde_con[i]*u_tilde_con[j] ;
  u_tilde_sq = fabs(u_tilde_sq) ;

  gamma = sqrt(1. + u_tilde_sq) ;

  lapse = sqrt(-1./gcon[0][0]) ;

  for(i=0;i<4;i++) ucon[i] = u_tilde_con[i] - lapse*gamma*gcon[0][i] ;

  return ;
}

/**********************************************************************/
/******************************************************************
   
    raise_g():
 
         -- calculates the contravariant form of a covariant tensor, 
            using the inverse of the metric;

******************************************************************/
void raise_g(double vcov[4], double gcon[4][4], double vcon[4])
{
  int i,j;

  for(i=0;i<4;i++) {
    vcon[i] = 0. ;
    for(j=0;j<4;j++) 
      vcon[i] += gcon[i][j]*vcov[j] ;
  }

  return ;
}

/**********************************************************************/
/******************************************************************

     lower_g():
  
          -- calculates the ocvariant form of a contravariant tensor 
             using the metric;

******************************************************************/
void lower_g(double vcon[4], double gcov[4][4], double vcov[4])
{
  int i,j;

  for(i=0;i<4;i++) {
    vcov[i] = 0. ;
    for(j=0;j<4;j++) 
      vcov[i] += gcov[i][j]*vcon[j] ;
  }

  return ;
}

/**********************************************************************/
/******************************************************************

     ncov_calc(): 

         -- calculates the covariant form of the normal vector to our 
            spacelike hypersurfaces ala the ADM formalism.

         -- requires the inverse metric;

******************************************************************/
void ncov_calc(double gcon[4][4],double ncov[4]) 
{
  double lapse ;
  int i;

  lapse = sqrt(-1./gcon[0][0]) ;

  ncov[0] = -lapse ;
  for( i = 1; i < 4; i++) { 
    ncov[i] = 0. ;
  }

  return ;
}

/**********************************************************************/
/******************************************************************

    bcon_calc_g(): 
  
        -- using the primitive variables, contra-/co-variant 4-vel., 
           and covariant normal vector, calculate the contravariant 
           form of the magnetic 4-vector b^\mu (the small "b" in HARM);
       -- assumes:

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /


******************************************************************/
void bcon_calc_g(double prim[8],double ucon[4],double ucov[4],
		 double ncov[4],double bcon[4]) 
{
  static double Bcon[4] ;
  double u_dot_B ;
  double gamma ;
  int i ;

  // Bcon = \mathcal{B}^\mu  of the paper:
  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = -ncov[0] * prim[BCON1+i-1] ;

  u_dot_B = 0. ;
  for(i=0;i<4;i++) u_dot_B += ucov[i]*Bcon[i] ;

  gamma = -ucon[0]*ncov[0] ;
  for(i=0;i<4;i++) bcon[i] = (Bcon[i] + ucon[i]*u_dot_B)/gamma ;
}


/**********************************************************************/
/******************************************************************

    gamma_calc_g(): 
  
        -- using the primitive variables, contra-/co-variant 4-vel., 
           and covariant normal vector, calculate the contravariant 
           form of the magnetic 4-vector b^\mu (the small "b" in HARM);

       -- assumes:

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /

******************************************************************/
int gamma_calc_g(double *pr, double gcov[4][4], double *gamma)
{
  double utsq ;
  int j,k;

  utsq =    gcov[1][1]*pr[UTCON1]*pr[UTCON1]
    + gcov[2][2]*pr[UTCON2]*pr[UTCON2]
    + gcov[3][3]*pr[UTCON3]*pr[UTCON3]
    + 2.*(  gcov[1][2]*pr[UTCON1]*pr[UTCON2]
	    + gcov[1][3]*pr[UTCON1]*pr[UTCON3]
	    + gcov[2][3]*pr[UTCON2]*pr[UTCON3]) ;

  if(utsq<0.0){
    if(fabs(utsq)>1.E-10){ // then assume not just machine precision
      return (1);
    }
    else utsq=fabs(utsq); // set floor
  }

  *gamma = sqrt(1. + utsq) ;

  return(0) ;
}


/**************************************************/
/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/

/* 
pressure as a function of rho0 and u 
this is used by primtoU and Utoprim_?D
*/
double pressure_rho0_u(double rho0, double u)
{
  return((GAMMA - 1.) * u) ;
}


  
/* 
pressure as a function of rho0 and w = rho0 + u + p 
this is used by primtoU and Utoprim_1D
*/
double pressure_rho0_w(double rho0, double w)
{
  return((GAMMA - 1.)*(w - rho0)/GAMMA) ;
}

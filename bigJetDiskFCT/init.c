#include "pluto.h"

real vler(real dbl,real dbr);
real mmod(real dbl,real dbr);

void compute_pdisk(real p0d, real b0d, real mb, real hitst, int nptst, real *ztst, real *ptst, real *vtst);
real zforce(real z, real pnow, real p0d, real b0d, real mb);
real rforce(real z, real pnow, real p0d, real b0d, real mb);

/* ************************************************************** */
void Init (double *us, double x1, double x2, double x3)
/* 
 *
 * 
 * NAME 
 * 
 *   INIT
 *
 *
 * PURPOSE
 *
 *   Set inital condition. 
 *
 *
 * ARGUMENTS 
 *
 *   us (OUT)           a pointer to a array_1D of primitive 
 *                      (when *prim_var = 1) or conservative
 *                      (when *prim_var = 0) variables;
 *
 *   x1, x2, x3 (IN)    the three coordinates; the meaning is different
 *                      depending on the geometry:
 *
 *                       x1      x2      x3
 *                       -------------------                
 *                        x       y       z     in cartesian geometry
 *                        r       z       -     in cylindrical geometry 
 *                        r      phi      z     in polar geometry 
 *                        r     theta    phi    in spherical geometry
 *                  
 *   i , j , k  (IN)    the integer indexes of the cell containing the 
 *                      point x1, x2, x3.
 *                
 *   *prim_var (OUT)    an integer flag. When *prim_var = 1 initial 
 *                      conditions are assigned in terms of 
 *                      primitive values. When *prim_var = 0 initial
 *                      conditions are given in terms of conservative 
 *                      value.
 *
 *
 * Variable names are accessed as us[nv], where
 *
 *   nv = RHO is density
 *   nv = PRS is pressure
 *  
 *   Vector components are labelled always as VX1,VX2,VX3, but alternative
 *   labels may be used:
 *
 *   Cartesian       VX1     VX2      VX3
 *   Cylindrical    iVR    iVX3     iVPHI
 *   Polar          iVR    iVPHI   iVX3
 *   Spherical      iVR    iVTH    iVPHI
 *
 * 
 *
 **************************************************************** */
{

  int ind, iloc;
  real eps, d0d, p0d, d0a, p0a, b0d, mb, alpha;
  real ddisk, pdisk, vdisk, datmo, patmo, scrh;
  real rad, subk; 
  real aa, bb, xloc, dzloc, func, vunc;

  static int   nptst, indx;
  static real  hitst;
  static real *ptst, *vtst, *ztst;

  eps = 0.1;
  d0d = 1.;
  p0d = eps*eps;
  d0a = 1.e-4;
  p0a = d0a*(g_gamma-1.)/g_gamma;
  b0d = sqrt(2.*p0d*0.3);
  b0d = 0;
  mb = 0.935;
  alpha = 0.9;
  subk = 1. - eps*eps*g_gamma/(g_gamma-1.);
 
  nptst = 10000;
  if (ptst == NULL) {
   hitst = 1.1*sqrt(1.-subk*subk)/subk;
   ptst = array_1D(nptst);
   vtst = array_1D(nptst);
   ztst = array_1D(nptst);
   compute_pdisk(p0d,b0d,mb,hitst,nptst,ztst,ptst,vtst);
    for (ind=0 ; ind<nptst ; ind++) {
     if (ptst[ind] > 0.0) indx = ind; 
    }
  }
 
  rad = sqrt(x1*x1+x2*x2);
  
  datmo = d0a*pow(rad,-3./2.);
  patmo = p0a*pow(rad,-5./2.);

  xloc = fabs(x2/x1);
  xloc = dmin(xloc,ztst[nptst-1]);
 dzloc = hitst/(real)nptst;
  iloc = (int)(xloc/dzloc);

  func = 0.0;
  vunc = 0.0;

  if (iloc < nptst-1) {
   if ((xloc > ztst[iloc+1]) || (xloc < ztst[iloc]))  printf("%d %f %f %f \n",iloc,xloc,ztst[iloc],ztst[iloc+1]);
   aa = (ptst[iloc+1]-ptst[iloc])/(ztst[iloc+1]-ztst[iloc]);
   bb = (ztst[iloc+1]*ptst[iloc]-ztst[iloc]*ptst[iloc+1])/(ztst[iloc+1]-ztst[iloc]); 
   func = aa*x2*p0d/pow(x1,7./2.)+bb*p0d/pow(x1,5./2.);
   func = dmax(func,0.0);
   aa = (vtst[iloc+1]-vtst[iloc])/(ztst[iloc+1]-ztst[iloc]);
   bb = (ztst[iloc+1]*vtst[iloc]-ztst[iloc]*vtst[iloc+1])/(ztst[iloc+1]-ztst[iloc]);
   vunc = aa*x2/pow(x1,3./2.)+bb/pow(x1,1./2.);
   vunc = dmax(vunc,0.0);
  }

/* Simplified: only hydrostatic equilibrium */ 
  scrh = (1./rad-subk/x1)*(g_gamma-1.)/g_gamma/p0d;
  ddisk = d0d*pow(dmax(scrh,0.),3./2.);
  pdisk = p0d*pow(dmax(scrh,0.),5./2.);
  vdisk = sqrt(subk/x1);

/* Complicated... */
  pdisk = func;
  ddisk = pow(pdisk/p0d,3./5.);
  vdisk = vunc;

  /* Kluzniak Kita */
  real r=x1;
  real z=x2;
  real h=sqrt(5.0)*0.1*x1;
real rm1=x1-20./256.;
real hm1=0.1*(rm1);

 // real alpha=1.0;
  real lambda=(11./5.)/(1+ 64.*alpha*alpha/25.);
	real omega0=pow(r,-1.5);
  real omega1=0;
  real dlnhdlnr= (log(h)-log(hm1) ) /log(r)-log(rm1);
  dlnhdlnr= 1;
  real omega2=omega0*pow(h/r,2.)*(-3./4. + 0.5*dlnhdlnr+2./15.*alpha*alpha*lambda*(1-6*z*z/h/h));

  real cs0=sqrt((h*h-z*z)/3.0/pow(r,3.0));
  real cs=cs0;
  real u0=0;
  real u1=-alpha*(h*h/pow(r,5./2.))*(-lambda*(1-z*z/h/h) -lambda*(32*alpha*alpha/15.) + 2*dlnhdlnr);


   vdisk +=  pow(h/r,2.)* 2./15.*alpha*alpha*lambda*(1-6*z*z/h/h)/sqrt(r);
   real vxdisk= u1;
   real vydisk=  z*vxdisk/r;


  if ( pdisk > patmo)
   {
      us[RHO] = ddisk;
      us[PRS] = pdisk;
      us[VX1] = vxdisk;
      us[VX2] = vydisk;
      us[VX3] = vdisk;
   }
   else
   {
      us[RHO] = datmo;
      us[RHO] = dmax(datmo, 7e-8 );
      us[PRS] = patmo;
      us[PRS] = dmax(patmo, 2e-10 );
      us[VX3] = 0.0;
   }

/*  if (fabs(pdisk-patmo)/patmo < 1.) 
   {
      us[RHO] = 0.5*(ddisk+datmo);
      us[PRS] = 0.5*(pdisk+patmo);
      us[VX3] = 0.5*vdisk;
   } */
      //us[VX1] = 0.;
      //us[VX2] = 0.;

      scrh = (2.*x2*x2+3./4.*mb*mb*x1*x1)/pow(x2*x2+mb*mb*x1*x1,13./8.);
      us[BX2] = scrh*4./3.*b0d*pow(mb,5./4.);
      scrh = 5./4.*x1*x2/pow(x2*x2+mb*mb*x1*x1,13./8.);
      us[BX1] = scrh*4./3.*b0d*pow(mb,5./4.);

		b0d=1e-2;
   #ifdef STAGGERED_MHD
      us[AX] = 0.0;
      us[AY] = 0.0;
		real sm=9.0; /*smoothing radius */
		real mb2=mb*mb;
      //scrh = x1/pow(x2*x2+mb*mb*x1*x1+sm*sm,5./8.);
      //us[AZ] = scrh*4./3.*b0d*pow(mb,5./4.);
		real  alpha_psi = 0.5;
		real   beta =  (2.-alpha_psi)/2.;
		  real const1;
			real rint=6.0;
		   const1 =  pow(mb,beta)/alpha_psi* b0d*pow(rint,2-alpha_psi);

		b0d=0.01;
      us[AZ] = b0d* pow(r,alpha_psi+2*beta-1)*pow(mb2*r*r +z*z, -beta);
		/*
		if (x1<1){
			  us[AZ] = 4./3.*b0d*pow(mb,5./4.);
		}
		*/
   #endif

      us[BX3] = 0.0;
 
}
/* **************************************************************** */
#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 * NAME 
 *
 *   GRAVITY
 *
 *
 * PURPOSE
 *
 *   Define gravitational acceleration array_1D as 
 *   function of the three spatial coordinates
 *
 * 
 * ARGUMENTS
 *
 *   w       (IN)      a array_1D containing cell-centered primitive
 *                     values ofthe state array_1D;
 *
 *   x1,x2,x3 (IN)     the three coordinates;
 *  
 *   g        (OUT)    a array_1D containing the three components
 *                     of the gravitational acceleration.
 * 
 **************************************************************** */
{

  real rad, scrh;

  rad = sqrt(x1*x1+x2*x2);
  scrh = rad*rad*rad;

  g[IDIR] = -x1/scrh;
  g[JDIR] = -x2/scrh;
  g[KDIR] = 0.0;

}
#endif
#if BACKGROUND_FIELD == YES
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/* 
 *
 *
 * NAME
 * 
 *   BACKGROUND_FIELD
 *
 *
 * PURPOSE
 *
 *   Define the component of a static, curl-free background 
 *   magnetic field.
 *
 *
 * ARGUMENTS
 *
 *   x1, x2, x3  (IN)    coordinates
 *
 *   B0         (OUT)    array_1D component of the background field.
 *
 *
 **************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif


/* ************************************************************** */
//void USERDEF_BOUNDARY (real ***uu[], int  side,  int grid_type, 
//                       int i0, int i1, int j0, int j1,
//                       int k0, int k1, Grid *grid) 
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 * 
 * NAME
 *   
 *   USERDEF_BOUNDARY
 *
 *
 * PURPOSE
 *
 *   Define user-defined boundary conditions.
 *
 *
 * ARGUMENTS
 * 
 *   uu (IN/OUT)       a three-dimensional array_1D uu[nv][k][j][i], 
 *                     containing the solution and boundary values 
 *                     to be assigned.
 *
 *   side (IN)         tells on which side boundary conditions need 
 *                     to be assigned. side can assume the following 
 *                     pre-definite values: X1_BEG, X1_END,
 *                                          X2_BEG, X2_END, 
 *                                          X3_BEG, X3_END
 *
 *                     When 
 *
 *                      side = Xn_BEG  -->  b.c. are assigned at the 
 *                                          beginning of the xn
 *                                          direction.
 *                      side = Xn_BEG  -->  b.c. are assigned at the 
 *                                          beginning of the xn
 *                                          direction.
 *
 *   grid_type (IN)    Tells whether boundary conditions have
 *                     to be assigned on cell-centered quantities
 *                     (grid_type == CELL_CENTER) or 
 *                     face centered (grid_type == FACE_CENTER),
 *                     as it is the case for the magnetic field
 *                     components in the Flux-CT algoritm.          
 *
 *   i0, i1   (IN)     first and last index spanning the 
 *                     x1 direction. When side = X1_BEG, or 
 *                     side = X1_END, i0 and i1 are the leftmost and
 *                     rightmost points in the boundary. 
 *                     Otherwise i0 and i1 are the first and last points
 *                     inside the domain;
 *
 *   j0, j1   (IN)    first and last index spanning the x2 direction;
 *
 *   k0, k1   (IN)    first and last index spanning the x3 direction;
 *
 *   grid    (IN)     a pointer to an array of grid structures.
 *                     
 *
 *
 **************************************************************** */
{
  int   i, j, k, nv;
  real  x1, x2, x3;
  real  dul, dur, du, slpb, slpv;
  real slpbz;
  real slpd;
  real slppr;
  real  scrh1, scrh2, scrh3, scrh4;



  if (side == X1_BEG){

    if (grid_type == CELL_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}

    }else if (grid_type == FACE_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}

      i = IBEG - 2;
      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k]; 
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j]; 
 
      }}

    }

  } else if (side == X1_END) {

    if (grid_type == CELL_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  


   dur = uu[BX2][k][j+1][i0-1]-uu[BX2][k][ j ][i0-1];
   dul = uu[BX2][k][ j ][i0-1]-uu[BX2][k][j-1][i0-1];
   if (j == j0) dul = 0.0;
   du = vler(dul,dur);
   scrh1 = -du-uu[BX1][k][j][i0-1]/grid[IDIR].x[i0-1]*grid[IDIR].dx[i0-1];
   dur = uu[BX1][k][j][i0-2]-uu[BX1][k][j][i0-3];
   dul = uu[BX1][k][j][i0-1]-uu[BX1][k][j][i0-2];
   scrh2 = vler(dul,dur);
   du = vler(scrh1,scrh2);

   scrh1 = grid[IDIR].x[i0-1]/grid[IDIR].x[i0-2];
   scrh2 = grid[IDIR].x[i0-2]/grid[IDIR].x[i0-3];

   slpv = 0.0;
  if ( (uu[VX3][k][j][i0-1] > 0.0) && (uu[VX3][k][j][i0-2] > 0.0)  && (uu[VX3][k][j][i0-3] > 0.0) ) {
   scrh3 = uu[VX3][k][j][i0-1]/uu[VX3][k][j][i0-2];
   scrh4 = uu[VX3][k][j][i0-2]/uu[VX3][k][j][i0-3];
   dul = log(scrh3)/log(scrh1);
   dur = log(scrh4)/log(scrh2);
   slpv = mmod(dul,dur);
   slpv = dmin(slpv,0.0);
  }



   slpbz = 0.0;
  if ( (uu[BX2][k][j][i0-1] > 0.0) && (uu[BX2][k][j][i0-2] > 0.0)  && (uu[BX2][k][j][i0-3] > 0.0) ) {
   scrh3 = uu[BX2][k][j][i0-1]/uu[BX2][k][j][i0-2];
   scrh4 = uu[BX2][k][j][i0-2]/uu[BX2][k][j][i0-3];
   dul = log(scrh3)/log(scrh1);
   dur = log(scrh4)/log(scrh2);
   slpbz = mmod(dul,dur);
   slpbz = dmin(slpbz,0.0);
  }


  slpb = 0.0;
  if ( (uu[BX3][k][j][i0-1] < 0.0) && (uu[BX3][k][j][i0-2] < 0.0)  && (uu[BX3][k][j][i0-3] < 0.0) ) {
   scrh3 = uu[BX3][k][j][i0-1]/uu[BX3][k][j][i0-2];
   scrh4 = uu[BX3][k][j][i0-2]/uu[BX3][k][j][i0-3];

   dul = log(scrh3)/log(scrh1);
   dur = log(scrh4)/log(scrh2);
   slpb = mmod(dul,dur);
   slpb = dmin(slpb,0.0);
  }



   slppr = 0.0;
   scrh3 = uu[PRS][k][j][i0-1]/uu[PRS][k][j][i0-2];
   scrh4 = uu[PRS][k][j][i0-2]/uu[PRS][k][j][i0-3];
   dul = log(scrh3)/log(scrh1);
   dur = log(scrh4)/log(scrh2);
   slppr = mmod(dul,dur);
   slppr = dmin(slppr,0.0);

   slpd = 0.0;
   scrh3 = uu[RHO][k][j][i0-1]/uu[RHO][k][j][i0-2];
   scrh4 = uu[RHO][k][j][i0-2]/uu[RHO][k][j][i0-3];
   dul = log(scrh3)/log(scrh1);
   dur = log(scrh4)/log(scrh2);
   slpd = mmod(dul,dur);
   slpd = dmin(slpd,0.0);


if (uu[RHO][k][j][i0-1] < 1.e-5 ) {
        slpd = 0.0;
        slppr = 0.0;
      }

      slpd = dmax(slpd,3./5.*slppr);
      slppr = 5./3.*slpd;

   for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  


       scrh1 = x1/grid[IDIR].x[i0-1];
       uu[VX1][k][j][i] = uu[VX1][k][j][i0-1];
       uu[VX2][k][j][i] = uu[VX2][k][j][i0-1];
       if (x2 < 0.1*x1) {
         uu[VX1][k][j][i] = dmin(uu[VX1][k][j][i],0.0);
         uu[VX2][k][j][i] = dmin(uu[VX2][k][j][i],0.0);
       }
       uu[VX3][k][j][i] = dmax(uu[VX3][k][j][i0-1],0.0)*pow(scrh1,slpv);


      real keplerian=sqrt(1/x1);
       uu[VX3][k][j][i] = dmin( uu[VX3][k][j][i], keplerian);
      // uu[VX3][k][j][i] = dmax( uu[VX3][k][j][i], 0.99*keplerian);

       uu[BX1][k][j][i] = uu[BX1][k][j][i0-1]+du*(i-i0+1);
       uu[BX2][k][j][i] = dmax(uu[BX2][k][j][i0-1],0.0)*pow(scrh1,slpbz);
       uu[BX3][k][j][i] = dmin(uu[BX3][k][j][i0-1],0.0)*pow(scrh1,slpb);
       //uu[PRS][k][j][i] = uu[PRS][k][j][i0-1]*pow(scrh1,slppr);
       //uu[RHO][k][j][i] = uu[RHO][k][j][i0-1]*pow(scrh1,slpd);
       uu[PRS][k][j][i] = uu[PRS][k][j][i0-1];//*pow(scrh1,slppr);
       uu[RHO][k][j][i] = uu[RHO][k][j][i0-1];//pow(scrh1,slpd);

         
      }}}

    }else if (grid_type == FACE_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}

      i = IEND + 1;
      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k]; 
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j]; 
 
      }}

    }

  } else if (side == X2_BEG){

    if (grid_type == CELL_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}

    }else if (grid_type == FACE_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}

      j = JBEG - 2;
      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];

      }}

    }

  } else if (side == X2_END) {

    if (grid_type == CELL_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}

    }else if (grid_type == FACE_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}

      j = JEND + 1;
      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
 
      }}

    }

  } else if (side == X3_BEG){

    if (grid_type == CELL_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}

    }else if (grid_type == FACE_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}
 
      k = KBEG - 2;
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];

      }}

    }

  } else if (side == X3_END) {

    if (grid_type == CELL_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}

    }else if (grid_type == FACE_CENTER){

      for (k = k0; k <= k1; k++) { x3 = grid[KDIR].x[k];  
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];  
         
      }}}

      k = KEND + 1;
      for (j = j0; j <= j1; j++) { x2 = grid[JDIR].x[j];  
      for (i = i0; i <= i1; i++) { x1 = grid[IDIR].x[i];

      }}

    }
  }

}


#if ADD_INTERNAL_BOUNDARY == YES
/* ************************************************************** */

void INTERNAL_BOUNDARY(real ***uu[], real ***uus[],
                       Grid *grid, 
                       int i0, int i1, int j0, int j1, int k0, int k1)

/* 
 *
 *
 * NAME
 *
 *   INTERNAL_BOUNDARY
 *
 *
 * PURPOSE
 *
 *   Allow the user to control 
 *
 *
 * ARGUMENTS
 *
 *   uu      (IN/OUT)    three-dimensional array_1D containing the solution data;
 *
 *   uu_old  (IN/OUT)    old, kept for backward compatibility;
 *
 *   grid    (IN)        pointer to grid structures;
 * 
 *   i0, j0, k0 (IN)     indexes of the lower-coordinate point inside
 *                       the domain; 
 *
 *   i1, j1, k1, (IN)    indexes of the upper-coordinate point inside
 *                       the domain.
 *   
 *  
 *
 *
 **************************************************************** */
{
  int  i, j, k;
  real x1, x2, x3;
  real btotp, btot, scrh;
  real csound, calfv, caltot, cfast, clim, tp;

  real dur, dul, du;
  real dud, dup, dubx, duby, dubz;
  real duvx;
  real duvy;
  real duvz;

 
  int nptsx, nptsy;

  clim = 100.;
  clim = 30.;
  // In order not to have a decrease in timestep try making clim smaller ...
  //clim = 5.;
  nptsx = NPTSX;
  nptsy = NPTSY;
	


  for (k = k0; k <= k1; k++)
    {
      x3 = grid[KDIR].x[k];
      for (j = j0; j <= j1; j++)
   {
     x2 = grid[JDIR].x[j];
     for (i = i0; i <= i1; i++)
       {
         x1 = grid[IDIR].x[i];
         if (  uu[PRS][k][j][i] <  0.01*pow(uu[RHO][k][j][i], g_gamma)    )
         {
         uu[PRS][k][j][i] =  0.01*pow(uu[RHO][k][j][i], g_gamma) ;
         };

       }
   }
    }

	
  if ( (  grid[IDIR].is_gbeg == 1 ) && (grid[JDIR].is_gbeg == 1  ) )  {
  


  k = 0;
  j = j0+nptsy;
  for (i = i0; i <= i0+2-1; i++) {
   csound = 5./3.*uu[PRS][k][j][i]/uu[RHO][k][j][i];
   btotp = uu[BX1][k][j][i]*uu[BX1][k][j][i]+uu[BX2][k][j][i]*uu[BX2][k][j][i];
   btot  = btotp+uu[BX3][k][j][i]*uu[BX3][k][j][i];
   calfv = btotp/uu[RHO][k][j][i];
   caltot = btot/uu[RHO][k][j][i]; 

   cfast = 0.5*(csound+caltot+sqrt((csound+caltot)*(csound+caltot)-4.*csound*calfv)); 

   tp = uu[PRS][k][j][i]/uu[RHO][k][j][i]; 
   if (cfast > clim) {
	uu[RHO][k][j][i] = uu[RHO][k][j][i]*cfast/clim;
	uu[RHO][k][j][i]=dmax (1e-9,uu[RHO][k][j][i]);
	}

   uu[RHO][k][j][i]=1e-5;
    uu[VX2][k][j][i] = 1.01;
   uu[PRS][k][j][i] = tp*uu[RHO][k][j][i];
  }
   
  for (k = k0; k <= k1; k++)      { x3 = grid[KDIR].x[k];  
  for (j = j0; j <= j0+nptsy-1; j++) { x2 = grid[JDIR].x[j];  


   dur = uu[RHO][k][j][i0+nptsx+2]-uu[RHO][k][j][i0+nptsx+1];
   dul = uu[RHO][k][j][i0+nptsx+1]-uu[RHO][k][j][i0+nptsx+0];
   dud = vler(dul,dur);
   dud = dmin(dud, 0.0);

   dur = uu[PRS][k][j][i0+nptsx+2]-uu[PRS][k][j][i0+nptsx+1];
   dul = uu[PRS][k][j][i0+nptsx+1]-uu[PRS][k][j][i0+nptsx+0];
   dup = vler(dul,dur);
	/*
	 * Limit on dup
	 * */
	real duptest=0;
	real vphi= uu[VX3][k][j][i0+nptsx];
	x1 = grid[IDIR].x[i0+nptsx];
	duptest= uu[RHO][k][j][i0+nptsx]*(-x1/pow(x1*x1+x2*x2, 1.5) + vphi*vphi/x1)*grid[IDIR].dx[i0+nptsx];
	dup = dmax( dup,duptest);
   dup = dmin(dup, 0.0);

   dur = uu[VX1][k][j][i0+nptsx+2]-uu[VX1][k][j][i0+nptsx+1];
   dul = uu[VX1][k][j][i0+nptsx+1]-uu[VX1][k][j][i0+nptsx+0];
   duvx = vler(dul,dur);
   duvx = dmin(duvx, 0.0);

   dur = uu[VX2][k][j][i0+nptsx+2]-uu[VX2][k][j][i0+nptsx+1];
   dul = uu[VX2][k][j][i0+nptsx+1]-uu[VX2][k][j][i0+nptsx+0];
   duvy = vler(dul,dur);
   duvy = dmin(duvy, 0.0);

   real eps=0.1;
   real subk= 1. - eps*eps*g_gamma/(g_gamma-1.);
   dur = uu[VX3][k][j][i0+nptsx+2]-uu[VX3][k][j][i0+nptsx+1];
   dul = uu[VX3][k][j][i0+nptsx+1]-uu[VX3][k][j][i0+nptsx+0];
   duvz = vler(dul,dur);
   duvz = dmin(duvz, 0.0);

  for (i = i0; i <= i0+nptsx-1; i++) { x1 = grid[IDIR].x[i];  

/*
    real scrh7=sqrt(  uu[PRS][k][j][i0+nptsx]/ uu[RHO][k][j][i0+nptsx] * x1*x1*x1) ;
    if ( x2 < scrh7  )
    {
    uu[VX1][k][j][i] = dmin( uu[VX1][k][j][i0+nptsx], 0.0);
    uu[VX2][k][j][i] = dmin( uu[VX2][k][j][i0+nptsx], 0.0);
    }
*/
    scrh = (real)dmin(i0+nptsx-i,3);

	 uu[VX1][k][j][i] =  uu[VX1][k][j][i0+nptsx]-duvx*scrh;
	 uu[VX2][k][j][i] =  uu[VX2][k][j][i0+nptsx]-duvy*scrh;
    uu[VX1][k][j][i] = dmin( uu[VX1][k][j][i], 0.0);
	 
	 real eps=0.1;
	 real subk= 1. - eps*eps*g_gamma/(g_gamma-1.);
	 uu[VX3][k][j][i] =  uu[VX3][k][j][i0+nptsx]-duvz*scrh;
	 uu[VX3][k][j][i] =  dmin( uu[VX3][k][j][i], sqrt(subk/x1));

	 uu[RHO][k][j][i] =  uu[RHO][k][j][i0+nptsx]-dud*scrh;
	 uu[PRS][k][j][i] =  uu[PRS][k][j][i0+nptsx]-dup*scrh;

	 /*
	 real rad= sqrt(x1*x1+x2*x2);
	 real eps=0.1;
	 real subk= 1. - eps*eps*g_gamma/(g_gamma-1.);
	 real p0d= eps*eps;
	 real d0d= 1.0;
    real scrh = (1./rad-subk/x1)*(g_gamma-1.)/g_gamma/p0d;
  	 real ddisk = d0d*pow(dmax(scrh,0.),3./2.);
  	 real pdisk = p0d*pow(dmax(scrh,0.),5./2.);

  real d0a = 1.e-4;
  real p0a = d0a*(g_gamma-1.)/g_gamma;
  	real	datmo = d0a*pow(rad,-3./2.);
  	real 	patmo = p0a*pow(rad,-5./2.);
	 if (patmo > pdisk)
	 {
	 uu[RHO][k][j][i] =  datmo;
	 uu[PRS][k][j][i] =  patmo;
	 }
	 else{
	 uu[RHO][k][j][i] =  ddisk;
	 uu[PRS][k][j][i] =  pdisk;
	 //uu[VX3][k][j][i] =  sqrt(subk/x1);
	 }
	 //zzzzzzzzzzzzzzzzzzzz
	 */

    //uu[PRS][k][j][i] =  uu[PRS][k][j][i0+nptsx];
    //uu[PRS][k][j][i] =   uu[PRS][k][j][i0+nptsx] -fabs (uu[PRS][k][j][i0+nptsx]-uu[PRS][k][j][i0+nptsx+1]) ;
	 if (j > j0+nptsy-3 && i < i0+nptsx-2)
	 {
    uu[RHO][k][j][i] =  uu[RHO][k][j0+nptsy][i];
    uu[VX1][k][j][i] =  dmin(0.0,uu[VX1][k][j0+nptsy][i]);
    uu[VX2][k][j][i] =  dmin(0.0,uu[VX2][k][j0+nptsy][i]);
    uu[VX3][k][j][i] =  dmax(0.0,uu[VX3][k][j0+nptsy][i]);
    uu[PRS][k][j][i] =  uu[PRS][k][j0+nptsy][i];
	 }

  }}}

/* Assign internal boundary conditions for staggered field */
#ifdef STAGGERED_MHD

  for (k = k0; k <= k1; k++)           { x3 = grid[KDIR].x[k];
  for (j = j0-1; j <  j0+nptsy-1; j++) { x2 = grid[JDIR].xr[j];
  for (i =  i0 ; i <= i0+nptsx-1; i++) { x1 = grid[IDIR].x[i];
   scrh = dmin(grid[IDIR].x[i0+nptsx]/x1,1.2);
   uus[BX2s][k][j][i] = uus[BX2s][k][j][i0+nptsx]*pow(scrh,1.25); 
  }}}
 
  for (k = k0; k <= k1; k++)           { x3 = grid[KDIR].x[k];
  for (j =  j0 ; j <= j0+nptsy-1; j++) { x2 = grid[JDIR].x[j];
  for (i = i0-1; i <  i0+nptsx-1; i++) { x1 = grid[IDIR].xr[i];
   uus[BX1s][k][j][i] = uus[BX1s][k][j0+nptsy][i];
  }}}
     
  for (k = k0; k <= k1; k++)      { x3 = grid[KDIR].x[k];
  for (j = j0; j <= j0+nptsy-1; j++) { x2 = grid[JDIR].x[j];
  for (i = i0; i <= i0+nptsx-1; i++) { x1 = grid[IDIR].x[i];
   uu[BX1][k][j][i] = 0.5*(uus[BX1s][k][j][i]+uus[BX1s][k][j][i-1]);
   uu[BX2][k][j][i] = 0.5*(uus[BX2s][k][j][i]+uus[BX1s][k][j-1][i]);
  }}}
#endif    


  }
}

#endif

real mmod(real dbl,real dbr)
{
 real  db;

  if (dbl*dbr > 0.0){
      db = fabs(dbl) < fabs(dbr) ? dbl:dbr;
    }else{
      db = 0.0;
    }

  return db;
}


real vler(real dbl,real dbr)
{
 real s, db;

  s           = dbl*dbr;
  db = s > 0.0 ? 2.0*s/(dbl + dbr):0.0;

  return db;
}

void compute_pdisk(real p0d, real b0d, real mb, real hitst, int nptst, real *ztst, real *ptst, real *vtst)
{
 int i;
 real dz, pnow, paft, v2;
 real k1, k2, k3, k4;

 dz = hitst/(real)nptst;

 ptst[0] = 1.;
 ztst[0] = 0.;

 pnow = 1.0; 

 for (i=0 ; i < nptst-1 ; i++) {

  k1 = 0.0;
  k2 = 0.0;
  k3 = 0.0;
  k4 = 0.0;

  if (pnow > 0.0) {
   k1 = zforce(     ztst[i]   ,   pnow    ,p0d,b0d,mb);
   k1 = k1*dz;      
   k2 = zforce(ztst[i]+0.5*dz,pnow+0.5*k1,p0d,b0d,mb);
   k2 = k2*dz;
   k3 = zforce(ztst[i]+0.5*dz,pnow+0.5*k2,p0d,b0d,mb); 
   k3 = k3*dz;
   k4 = zforce(  ztst[i]+dz  ,  pnow+k3  ,p0d,b0d,mb); 
   k4 = k4*dz;
  }

  if (i == 0) {
   v2 = rforce(ztst[i],pnow,p0d,b0d,mb);
   vtst[i] = sqrt(dmax(v2,0.0));
  }

  paft = pnow + k1/6. + k2/3. + k3/3. + k4/6.;
  pnow = paft;

  ztst[i+1] = ztst[i]+dz;
  if ((paft > 0.0) && (paft == paft)) {
   ptst[i+1] = paft;
   v2 = rforce(ztst[i],pnow,p0d,b0d,mb);
   vtst[i+1] = sqrt(dmax(v2,0.0));
  } else {
   ptst[i+1] = 0.0;
   vtst[i+1] = 0.0;
  }
   
 }
 
}

real zforce(real z, real pnow, real p0d, real b0d, real mb)
{

 real  br, dbzr, dbrz, jp, grav;
 real  zforce;

 grav = -z/pow(1.+z*z,1.5);

 br = 5./4.*z/pow(z*z+mb*mb,13./8.);
 br = br*4./3.*b0d*pow(mb,5./4.);

 dbzr = (-5.*mb*mb*z*z-15./16.*pow(mb,4))/pow(z*z+mb*mb,21./8.);
 dbzr = dbzr*4./3.*b0d*pow(mb,5./4.);
 dbrz = (-45./16.*z*z+5./4.*mb*mb)/pow(z*z+mb*mb,21./8.); 
 dbrz = dbrz*4./3.*b0d*pow(mb,5./4.);

 jp = dbrz-dbzr;

 zforce = pow(pnow,3./5.)/p0d*grav-jp*br/p0d;

 return zforce;

}

real rforce(real z, real pnow, real p0d, real b0d, real mb)
{

 real bz, dbzr, dbrz, jp, grav;
 real rforce;

 grav = -1./pow(1.+z*z,1.5);

 bz = (2.*z*z+3./4.*mb*mb)/pow(z*z+mb*mb,13./8.);
 bz = bz*4./3.*b0d*pow(mb,5./4.);

 dbzr = (-5.*mb*mb*z*z-15./16.*pow(mb,4))/pow(z*z+mb*mb,21./8.);
 dbzr = dbzr*4./3.*b0d*pow(mb,5./4.);
 dbrz = (-45./16.*z*z+5./4.*mb*mb)/pow(z*z+mb*mb,21./8.);
 dbrz = dbrz*4./3.*b0d*pow(mb,5./4.);

 jp = dbrz-dbzr;

 rforce = -grav-jp*bz/pow(pnow,3./5.);
 rforce = rforce-5./2.*p0d*pow(pnow,2./5.)-p0d*zforce(z,pnow,p0d,b0d,mb)*z/pow(pnow,3./5.);

 return rforce;

}
/*
{
}
*/

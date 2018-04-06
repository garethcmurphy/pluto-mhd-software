#include "pluto.h"

#if INTERPOLATION != WENO5

/* *********************************************************** */
void GET_RHS (const State_1D *state, real dt, Grid *grid)
/* 
 *
 * PURPOSE:
 *
 *   compute right hand side of the PDE in different 
 *   geometries, taking contribution from one direction 
 *   at a time.
 *   Include geometrical source terms and forces.
 * 
 *
 *
 * LAST MODIFIED
 *
 *   Oct 1st, 2006 by A. Mignone.  e-mail: mignone@to.astro.it
 *
 ************************************************************* */
{
  int  i, nv, beg, end;
  int  nptsx=0, nptsy=0;
  real dtdx, dtdV, scrh;
  real **flux, **rhs, *p;

  if ( (  grid[IDIR].is_gbeg == 1 ) && (grid[JDIR].is_gbeg == 1  ) )  {

  nptsx = NPTSX;
  nptsy = NPTSY;
  }

  #if GEOMETRY != CARTESIAN
   static real *MTsrc_22, *MTsrc_33, *MTsrc_21;
   static real *ITsrc_31;
   static real **fA;

   real r1, r2, r3, r_1;
   real s1, s2, s_1, ct;
   real *r, *th;

   if (fA == NULL) {
     fA   = array_2D(NMAX_POINT, NVAR);
     MTsrc_22 = array_1D(NMAX_POINT);
     MTsrc_33 = array_1D(NMAX_POINT);
     MTsrc_21 = array_1D(NMAX_POINT);
     ITsrc_31 = array_1D(NMAX_POINT); /* -- used in cylindrical only -- */
   }
  #endif

  beg = grid[DIR].lbeg;
  end = grid[DIR].lend;

  rhs  = state->rhs;
  flux = state->flux;
  p    = state->press;

/* -------------------------------------------------
             Build Right Hand Side
   ------------------------------------------------- */

  #if GEOMETRY == CARTESIAN

   for (i = beg; i <= end; i++) {
     dtdx = dt/grid[DIR].dx[i];
     for (nv = NVAR; nv--;  ) {
       rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i - 1][nv]);
     }
     rhs[i][M1] -= dtdx*(p[i] - p[i - 1]);
   }
  
  #elif GEOMETRY == CYLINDRICAL

   r = grid[IDIR].x;

   if (DIR == IDIR) {         /* *********************************************
                                 *          Cylindrical, r                   *
                                 ********************************************* */

   /* ------------------------------------------------
              Compute FLUX X AREA
      ------------------------------------------------ */

     for (i = beg - 1; i <= end; i++) {

       r1 = grid[IDIR].A[i];
       for (nv = NVAR; nv--; ) fA[i][nv] = flux[i][nv]*r1;
       
       #if COMPONENTS == 3
        fA[i][MZ] *= r1;
       #endif
     }

   /* ------------------------------------------------
         Compute momentum and induction 
         tensor source terms
      ------------------------------------------------ */

     #if COMPONENTS == 3
      MOM_SOURCE (state, beg, end, KDIR, KDIR, MTsrc_33, grid);
      #if PHYSICS == MHD || PHYSICS == RMHD
       OMEGA_SOURCE (state, beg, end, KDIR, IDIR, ITsrc_31, grid); 
      #endif
     #endif

   /* ------------------------------------------------
            Construct rhs operators, r - dir
      ------------------------------------------------ */

     for (i = beg; i <= end; i++) {
       dtdV = dt/grid[IDIR].dV[i];
       dtdx = dt/grid[IDIR].dx[i];
       r_1  = 1.0/r[i];

       for (nv = NVAR; nv--;   ){
         rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i - 1][nv]);
       }
       EXPAND(rhs[i][iMR] -= dtdx*(p[i] - p[i - 1]);    ,
                                                        ,
              rhs[i][iMPHI] *= r_1;)

       #if COMPONENTS == 3 
        rhs[i][iMR] += dt*MTsrc_33[i]*r_1; 
        #if  PHYSICS == MHD || PHYSICS == RMHD
         rhs[i][iBPHI] -= dt*ITsrc_31[i]*r_1;
       /* -------------------------------------------------------
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           in resistive MHD, we use the formulation without the 
           source terms, just because we cannot compute them!
           The other formulation (more stable for toroidal jet 
           configurations) will be used for ideal MHD
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ------------------------------------------------------- */
         #if RESISTIVE_MHD == YES
          rhs[i][iBPHI] = - dtdx*(flux[i][iBPHI] - flux[i - 1][iBPHI]);
         #endif
        #endif
       #endif

     }
     
   } else if (DIR == JDIR) {    /* *********************************************
                                   *         Cylindrical, z                    *
                                   ********************************************* */

      for (i = beg; i <= end; i++) {
       dtdx = dt/grid[DIR].dx[i];
       for (nv = NVAR; nv--;  ) {
         rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i - 1][nv]);
       }
       rhs[i][iMZ] -= dtdx*(p[i] - p[i - 1]);
     }
   }

  #elif GEOMETRY == POLAR

   r   = grid[IDIR].x;

   if (DIR == IDIR) {      /* *********************************************
                              *          Polar, r                         *
                              ********************************************* */

   /* ------------------------------------------------
              Compute FLUX X AREA
      ------------------------------------------------ */

     for (i = beg - 1; i <= end; i++) {

       r1 = grid[IDIR].A[i];
       r2 = r1*r1;

       fA[i][DN] = flux[i][DN]*r1;
       EXPAND(fA[i][MX] = flux[i][MX]*r1;  ,
              fA[i][MY] = flux[i][MY]*r2;  ,  /* -- ang mom conservation -- */
              fA[i][MZ] = flux[i][MZ]*r1;)

       #if PHYSICS == MHD || PHYSICS == RMHD
        EXPAND(fA[i][BX] = flux[i][BX]*r1;  ,
               fA[i][BY] = flux[i][BY]   ;  ,
               fA[i][BZ] = flux[i][BZ]*r1;)
       #endif
       #if EOS != ISOTHERMAL
        fA[i][EN] = flux[i][EN]*r1;
       #endif
       for (nv = NFLX; nv < NVAR; nv++) fA[i][nv] = flux[i][nv]*r1;
     }

   /* ------------------------------------------------
         Compute momentum tensor source terms
      ------------------------------------------------ */

     #if COMPONENTS >= 2
      MOM_SOURCE (state, beg, end, JDIR, JDIR, MTsrc_22, grid);
     #endif
           
   /* ------------------------------------------------
            Construct rhs operators, r - dir
      ------------------------------------------------ */

     for (i = beg; i <= end; i++) {
       dtdV = dt/grid[IDIR].dV[i];
       dtdx = dt/grid[IDIR].dx[i];
       r_1  = grid[IDIR].r_1[i];
    
       scrh = EXPAND(0.0, + MTsrc_22[i],   );

       rhs[i][DN] = - dtdV*(fA[i][DN] - fA[i - 1][DN]);
       #if EOS != ISOTHERMAL
        rhs[i][EN] = - dtdV*(fA[i][EN] - fA[i - 1][EN]);
       #endif 
       EXPAND(
         rhs[i][MX] = - dtdV*(fA[i][MX] - fA[i - 1][MX]) - dtdx*(p[i] - p[i - 1])
                      + dt*scrh*r_1;                                                ,
         rhs[i][MY] = - dtdV*(fA[i][MY] - fA[i - 1][MY])*r_1;                       ,
         rhs[i][MZ] = - dtdV*(fA[i][MZ] - fA[i - 1][MZ]);
       )
       #if PHYSICS == MHD || PHYSICS == RMHD
        EXPAND(                                                     
          rhs[i][BX] = -dtdV*(fA[i][BX] - fA[i - 1][BX]);  ,
          rhs[i][BY] = -dtdx*(fA[i][BY] - fA[i - 1][BY]);  ,
          rhs[i][BZ] = -dtdV*(fA[i][BZ] - fA[i - 1][BZ]);
        )
       #endif

       /* ---- scalars ---- */

       for (nv = NFLX; nv < NVAR; nv++) rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i - 1][nv]);
     }
     
   } else if (DIR == JDIR) {   /* *********************************************
                                  *          Polar, phi                       *
                                  ********************************************* */

      for (i = beg; i <= end; i++) {
        dtdx = dt/(r[*NX_PT]*grid[JDIR].dx[i]);
        for (nv = NVAR; nv--;  ) {
          rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i - 1][nv]);
        }
        rhs[i][iMPHI] -= dtdx*(p[i] - p[i - 1]);
      }

   } else if (DIR == KDIR) {    /* *********************************************
                                   *          Polar, z                         *
                                   ********************************************* */

      for (i = beg; i <= end; i++) {
       dtdx = dt/grid[KDIR].dx[i];
       for (nv = NVAR; nv--;  ) {
         rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i - 1][nv]);
       }
       rhs[i][iMZ] -= dtdx*(p[i] - p[i - 1]);
     }
 
   }

  #elif GEOMETRY == SPHERICAL

   r  = grid[IDIR].x;
   th = grid[JDIR].x;

   if (DIR == IDIR) {    /* *********************************************
                            *          Spherical, r                     *
                            ********************************************* */

   /* ------------------------------------------------
              Compute FLUX X AREA
      ------------------------------------------------ */

     for (i = beg - 1; i <= end; i++) {

       r1 = grid[IDIR].xr[i];
       r2 = r1*r1;
       r3 = r2*r1;

    /* ------------------------------------------------
         fphi -> fphi*r*sin(theta) , 
         but the sin(theta) factor cancels out at the 
         end (so we do not include it here)  
      -------------------------------------------------- */

       fA[i][DN] = flux[i][DN]*r2;
       EXPAND(fA[i][MX] = flux[i][MX]*r2;  ,
              fA[i][MY] = flux[i][MY]*r2;  ,
              fA[i][MZ] = flux[i][MZ]*r3;)

       #if PHYSICS == MHD || PHYSICS == RMHD
        EXPAND(fA[i][BX] = flux[i][BX]*r2;  ,
               fA[i][BY] = flux[i][BY]*r1;  ,
               fA[i][BZ] = flux[i][BZ]*r1;)
       #endif

       #if EOS != ISOTHERMAL
        fA[i][EN] = flux[i][EN]*r2;
       #endif
       for (nv = NFLX; nv < NVAR; nv++) fA[i][nv] = flux[i][nv]*r2;
     } 

   /* ------------------------------------------------
         Compute momentum tensor source terms
      ------------------------------------------------ */

     MOM_SOURCE (state, beg, end, JDIR, JDIR, MTsrc_22, grid);
     #if COMPONENTS == 3
      MOM_SOURCE (state, beg, end, KDIR, KDIR, MTsrc_33, grid);
     #endif
           
   /* ------------------------------------------------
            Construct rhs operators, r - dir
      ------------------------------------------------ */

     for (i = beg; i <= end; i++) {
       dtdV = dt/grid[IDIR].dV[i];
       dtdx = dt/grid[IDIR].dx[i];
       r_1  = grid[IDIR].r_1[i];

     /* ---------------------------------
          add geometrical source terms 
        --------------------------------- */

       scrh = EXPAND(0.0, + MTsrc_22[i], + MTsrc_33[i]);

       rhs[i][DN] = -dtdV*(fA[i][DN] - fA[i - 1][DN]);
       #if EOS != ISOTHERMAL
        rhs[i][EN] = -dtdV*(fA[i][EN] - fA[i - 1][EN]);
       #endif
       EXPAND(
         rhs[i][MX] = - dtdV*(fA[i][MX] - fA[i - 1][MX]) - dtdx*(p[i] - p[i - 1])
                      + dt*scrh*r_1;                                               ,
         rhs[i][MY] = - dtdV*(fA[i][MY] - fA[i - 1][MY]);                          ,
         rhs[i][MZ] = - dtdV*(fA[i][MZ] - fA[i - 1][MZ])*r_1;
       )
       
       #if PHYSICS == MHD || PHYSICS == RMHD
        EXPAND(                                                     
          rhs[i][BX] = -dtdV*(fA[i][BX] - fA[i - 1][BX]);      ,
          rhs[i][BY] = -dtdx*(fA[i][BY] - fA[i - 1][BY])*r_1;  ,
          rhs[i][BZ] = -dtdx*(fA[i][BZ] - fA[i - 1][BZ])*r_1;
        )
       #endif

       /* ---- scalars ---- */

       for (nv = NFLX; nv < NVAR; nv++) rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i - 1][nv]);
     }

   } else if (DIR == JDIR) {   /* *********************************************
                                  *          Spherical, Theta                 *
                                  ********************************************* */

   /* ------------------------------------------------
              Compute FLUX X AREA
      ------------------------------------------------ */

     for (i = beg - 1; i <= end; i++) {

       s1 = grid[JDIR].A[i];
       s2 = s1*s1;

       fA[i][DN] = flux[i][DN]*s1;
       EXPAND(fA[i][MX] = flux[i][MX]*s1;  ,
              fA[i][MY] = flux[i][MY]*s1;  ,
              fA[i][MZ] = flux[i][MZ]*s2;)  /* --  fphi -> fphi*sin(theta+) -- */

       #if PHYSICS == MHD || PHYSICS == RMHD
        EXPAND(fA[i][BX] = flux[i][BX]*s1;  ,
               fA[i][BY] = flux[i][BY]*s1;  ,
               fA[i][BZ] = flux[i][BZ];)
       #endif

       #if EOS != ISOTHERMAL
        fA[i][EN] = flux[i][EN]*s1;
       #endif
       for (nv = NFLX; nv < NVAR; nv++) fA[i][nv] = flux[i][nv]*s1;
     }
    
   /* ------------------------------------------------
         Compute momentum tensor source terms
      ------------------------------------------------ */

     MOM_SOURCE (state, beg, end, JDIR, IDIR, MTsrc_21, grid);
     #if COMPONENTS == 3
      MOM_SOURCE (state, beg, end, KDIR, KDIR, MTsrc_33, grid);
     #endif

   /* ------------------------------------------------
          Construct rhs operators, theta - dir
      ------------------------------------------------ */

     r_1 = grid[IDIR].r_1[*NX_PT];
     for (i = beg; i <= end; i++) {

       dtdV = dt/grid[JDIR].dV[i]*r_1;
       dtdx = dt/grid[JDIR].dx[i]*r_1;      

       s_1 = 1.0/sin(th[i]);   
       ct  = grid[JDIR].ct[i];         /* = cot(theta)  */

       scrh = EXPAND(0.0, - MTsrc_21[i], + ct*MTsrc_33[i]);  /* -- m[th] source -- */

       rhs[i][DN] = - dtdV*(fA[i][DN] - fA[i - 1][DN]);
       #if EOS != ISOTHERMAL
        rhs[i][EN] = - dtdV*(fA[i][EN] - fA[i - 1][EN]);
       #endif
       EXPAND(
         rhs[i][MX] = - dtdV*(fA[i][MX] - fA[i - 1][MX]);                          , 
         rhs[i][MY] = - dtdV*(fA[i][MY] - fA[i - 1][MY]) - dtdx*(p[i] - p[i - 1])
                      + dt*scrh*r_1;                                                , 
         rhs[i][MZ] = - dtdV*(fA[i][MZ] - fA[i - 1][MZ])*s_1;
       )

       #if PHYSICS == MHD || PHYSICS == RMHD
        EXPAND(                                                     
          rhs[i][BX] = - dtdV*(fA[i][BX] - fA[i - 1][BX]);  ,
          rhs[i][BY] = - dtdV*(fA[i][BY] - fA[i - 1][BY]);  ,
          rhs[i][BZ] = - dtdx*(fA[i][BZ] - fA[i - 1][BZ]);
        )
       #endif

       /* ---- scalars ---- */

       for (nv = NFLX; nv < NVAR; nv++) rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i - 1][nv]);
     }

   } else if (DIR == KDIR) {    /* *********************************************
                                   *          Spherical, Phi                   *
                                   ********************************************* */

     for (i = beg; i <= end; i++) {
       dtdx = dt/(r[*NX_PT]*sin(th[*NY_PT])*grid[DIR].dx[i]);       
       for (nv = NVAR; nv--;  ) {
         rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i - 1][nv]);
       }       
       rhs[i][iMPHI] -= dtdx*(p[i] - p[i - 1]);
     }
   }

  #endif

/* ----------------------------------------------------------------
                Resistive and Viscous source terms...
   ---------------------------------------------------------------- */

    for (i = beg; i <= end; i++){
     rhs[i][MX] += dt*state->src[i][MX];
     rhs[i][EN] += dt*state->src[i][EN];
    }

/* ----------------------------------------------------------------
                Powell's source terms
   ---------------------------------------------------------------- */

  #if PHYSICS == MHD || PHYSICS == RMHD 
   #if MHD_FORMULATION == EIGHT_WAVES
    for (i = beg; i <= end; i++) {
      EXPAND(rhs[i][MX] += dt*state->src[i][MX];  ,
             rhs[i][MY] += dt*state->src[i][MY];  ,
             rhs[i][MZ] += dt*state->src[i][MZ];) 

      EXPAND(rhs[i][BX] += dt*state->src[i][BX];  ,
             rhs[i][BY] += dt*state->src[i][BY];  ,
             rhs[i][BZ] += dt*state->src[i][BZ];) 
      #if EOS != ISOTHERMAL
       rhs[i][EN] += dt*state->src[i][EN];
      #endif
    }
   #endif
  #endif
    
/* ------------------------------------------------------------------
                        Include gravity
   ------------------------------------------------------------------ */

  #if INCLUDE_GRAVITY == YES
  
   ADD_GRAVITY (state, dt, grid);
   
  #endif

  if (DIR == IDIR && *NY_PT < grid[JDIR].lbeg+nptsy) { 
   for (i = beg; i <  beg+nptsx ; i++) {
    for (nv = NVAR; nv--;  ) state->rhs[i][nv] = 0.0;
  }}

  if (DIR == JDIR && *NX_PT < grid[IDIR].lbeg+nptsx) {
   for (i = beg; i <  beg+nptsy ; i++) {
    for (nv = NVAR; nv--;  ) state->rhs[i][nv] = 0.0;
  }}
}
#endif

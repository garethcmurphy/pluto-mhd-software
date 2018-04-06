#include "pluto.h"

#define ALL_SIDES -1
void  TempCalc ( real ***V[], Grid *grid, real *SoundSpeed , real *AlfvenSpeed, int first);


/* ************************************************************ */
void UNSPLIT (const Data *d, Riemann_Solver *RIEMANN, 
              real *cmax, real *eta_max, Grid *grid)

/*
 *
 *
 * PURPOSE: 
 *     
 *   Handle x1, x2 and x3 sweep integrators
 *
 *
 * ARGUMENTS:
 *
 *
 *
 *    
 *************************************************************** */
{
  int  nx, ny, nz, ii, jj, kk;
  int  *i, *j, *k;
  int  in,  in_beg, in_end;
  int  it1, it1_beg, it1_end;
  int  it2, it2_beg, it2_end;
  int  np_tot, nv;
  
    int first=0;

  static real  one_third = 1.0/3.0;
  static State_1D state;
  static Data_Arr UU, UU_1;

  static real *SoundSpeed, *AlfvenSpeed;

  real dt;

/* ----------------------------------------------------
                   Allocate memory 
   ---------------------------------------------------- */

  if (state.rhs == NULL){

  first = 1;

    nx = grid[IDIR].np_tot;
    ny = grid[JDIR].np_tot;
    nz = grid[KDIR].np_tot;
 
    Make_State (&state);
 
    UU   = array_4D(nz, ny, nx, NVAR);
    UU_1 = array_4D(nz, ny, nx, NVAR);





    SoundSpeed = array_1D(nx);
    AlfvenSpeed = array_1D(nx);

    #ifdef STAGGERED_MHD
     AVERAGE_MAGNETIC_FIELD (d->Vs, UU, grid);
     for (nv = BX; nv < BX + DIMENSIONS; nv++){
     for (kk = KBEG; kk <= KEND; kk++) {
     for (jj = JBEG; jj <= JEND; jj++) {
     for (ii = IBEG; ii <= IEND; ii++) {
       d->Vc[nv][kk][jj][ii] = UU[kk][jj][ii][nv];
     }}}}
    #endif
  }

/* ----------------------------------------------------
                   Step 1: predictor
   ---------------------------------------------------- */

  ISTEP = 1;  

  BOUNDARY (d, ALL_SIDES, grid);


/**
 *
 * Calculate temperature (+ Alfven speed ) at midplane and pass to all
 * processors above in the cartesian grid
 *
 * Call this function again in the corrector step
 *
 *    TempCalc ( real ***V[], Grid *grid, real *SoundSpeed , real *AlfvenSpeed)
 */
   TempCalc ( d->Vc, grid, SoundSpeed , AlfvenSpeed, first);


  for (DIR = 0; DIR < DIMENSIONS; DIR++){
  
  /*  --------------------------------
       Set normal and tangent indexes  
      -------------------------------- */

    SET_INDEXES(&np_tot, &i, &j, &k,
                &in,  &in_beg,  &in_end,
                &it1, &it1_beg, &it1_end,
                &it2, &it2_beg, &it2_end, grid);

    for (it2 = it2_beg; it2 <= it2_end; it2++) {
    for (it1 = it1_beg; it1 <= it1_end; it1++) {

      for (in = 0; in < np_tot; in++) {
      for (nv = NVAR; nv--;  ) {
        state.v[in][nv] = d->Vc[nv][*k][*j][*i];
      }}

      CHECK_NAN (state.v,0,np_tot-1,0);
      PRIMTOCON (state.v, state.u, 0, np_tot - 1);
      RECONSTRUCT (&state, grid);

      #ifdef STAGGERED_MHD
       for (in = 0; in < np_tot - 1; in++) {
         state.vL[in][B1] = state.vR[in][B1] = d->Vs[DIR][*k][*j][*i];
       }
      #endif

      RIEMANN (&state, cmax + DIR, grid);

      #ifdef RESISTIVE_MHD
       RESISTIVE_FLUX1 (d->Vc, state.flux, state.src, grid, eta_max, SoundSpeed, AlfvenSpeed);
      #endif

      #ifdef STAGGERED_MHD
       EMF_SAVE (state.flux, grid, -1);
      #endif

      #if NVAR != NFLX
       ADVECT_FLUX (&state, grid);
      #endif

      #if ARTIFICIAL_VISCOSITY == YES
       VISC_FLUX   (&state, grid);
      #endif

      GET_RHS (&state, delta_t, grid);

      if (DIR == IDIR){
        for (in = in_beg; in <= in_end; in++) {
        for (nv = NVAR; nv--;  ) {
          UU_1[*k][*j][*i][nv] = state.u[in][nv] + state.rhs[in][nv];
          UU  [*k][*j][*i][nv] = state.u[in][nv]; /* -- save initial value -- */
        }}
      }else{
        for (in = in_beg; in <= in_end; in++) {
        for (nv = NVAR; nv--;  ) {
          UU_1[*k][*j][*i][nv] += state.rhs[in][nv];
        }}
      }
    }}
  }

  #ifdef STAGGERED_MHD
   FCT_UPDATE(1, d, grid);
   AVERAGE_MAGNETIC_FIELD (d->Vs, UU_1, grid);
  #endif

  /* -------------------------
       convert to primitive
     ------------------------- */
     
  DIR = IDIR;
  SET_INDEXES(&np_tot, &i, &j, &k,
              &in,  &in_beg,  &in_end,
              &it1, &it1_beg, &it1_end,
              &it2, &it2_beg, &it2_end, grid);

  for (it2 = KBEG; it2 <= KEND; it2++){
  for (it1 = JBEG; it1 <= JEND; it1++){
    CONTOPRIM (UU_1[*k][*j], state.v,  d->flag_x[it2][it1], IBEG, IEND);
    for (in = IBEG; in <= IEND; in++){
    for (nv = 0; nv < NVAR; nv++) {
      d->Vc[nv][*k][*j][*i] = state.v[in][nv];
    }} 
  }}
  
/* ----------------------------------------------------
                   STEP II  (or CORRECTOR)
   ---------------------------------------------------- */

  #if TIME_EVOLUTION_STEPS >= 2

   ISTEP = 2;
   BOUNDARY (d, ALL_SIDES, grid);

   /*
    * Second call of tempcalc
    */

   TempCalc ( d->Vc, grid, SoundSpeed , AlfvenSpeed, first);


   for (kk = KBEG; kk <= KEND; kk++){
   for (jj = JBEG; jj <= JEND; jj++){
   for (ii = IBEG; ii <= IEND; ii++){
   for (nv = NVAR; nv--;  ) {
     #if TIME_EVOLUTION == RK2 
      UU_1[kk][jj][ii][nv] = 0.5*(UU[kk][jj][ii][nv] + UU_1[kk][jj][ii][nv]);
     #elif TIME_EVOLUTION == RK3 
      UU_1[kk][jj][ii][nv] = 0.75*UU[kk][jj][ii][nv] + 0.25*UU_1[kk][jj][ii][nv];
     #endif
   }}}}

   #if TIME_EVOLUTION == RK2
    dt = 0.5*delta_t;
   #elif TIME_EVOLUTION == RK3
    dt = 0.25*delta_t;
   #endif

   for (DIR = 0; DIR < DIMENSIONS; DIR++){

  /* --------------------------------
      Set normal and tangent indexes  
     -------------------------------- */

     SET_INDEXES(&np_tot, &i, &j, &k,
                 &in,  &in_beg,  &in_end,
                 &it1, &it1_beg, &it1_end,
                 &it2, &it2_beg, &it2_end, grid);

     for (it2 = it2_beg; it2 <= it2_end; it2++) {
     for (it1 = it1_beg; it1 <= it1_end; it1++) {

       for (in = 0; in < np_tot; in++) {
       for (nv = NVAR; nv--;  ) {
         state.v[in][nv] = d->Vc[nv][*k][*j][*i];
       }}

       PRIMTOCON   (state.v, state.u, 0, np_tot - 1);
       RECONSTRUCT (&state, grid);

       #ifdef STAGGERED_MHD
        for (in = 0; in < np_tot - 1; in++) {
          state.vL[in][B1] = state.vR[in][B1] = d->Vs[DIR][*k][*j][*i];
        }
       #endif

       RIEMANN (&state, cmax + DIR, grid);
      
       #ifdef RESISTIVE_MHD
       RESISTIVE_FLUX1 (d->Vc, state.flux, state.src, grid, eta_max, SoundSpeed, AlfvenSpeed);
       #endif

       #ifdef STAGGERED_MHD
        EMF_SAVE (state.flux, grid, -1);
       #endif

       #if NVAR != NFLX
        ADVECT_FLUX (&state, grid);
       #endif

       #if ARTIFICIAL_VISCOSITY == YES
        VISC_FLUX   (&state, grid);
       #endif

       GET_RHS (&state, dt, grid);
      
       for (in = in_beg; in <= in_end; in++) {
       for (nv = NVAR; nv--;  ) {
          UU_1[*k][*j][*i][nv] += state.rhs[in][nv];
       }}
     }}
   }

   #ifdef STAGGERED_MHD
    FCT_UPDATE(2, d, grid);
    AVERAGE_MAGNETIC_FIELD (d->Vs, UU_1, grid);
   #endif

  /* -------------------------
       convert to primitive
     ------------------------- */
     
   DIR = IDIR;
   SET_INDEXES(&np_tot, &i, &j, &k,
               &in,  &in_beg,  &in_end,
               &it1, &it1_beg, &it1_end,
               &it2, &it2_beg, &it2_end, grid);

   for (it2 = KBEG; it2 <= KEND; it2++){
   for (it1 = JBEG; it1 <= JEND; it1++){
     CONTOPRIM (UU_1[*k][*j], state.v,  d->flag_x[it2][it1], IBEG, IEND);
     for (in = IBEG; in <= IEND; in++){
     for (nv = 0; nv < NVAR; nv++) {
       d->Vc[nv][*k][*j][*i] = state.v[in][nv];
     }} 
   }}

  #endif

/* ----------------------------------------------------
                  STEP   III   (or CORRECTOR)
   ---------------------------------------------------- */

  #if TIME_EVOLUTION_STEPS == 3

   ISTEP = 3;
   BOUNDARY (d, ALL_SIDES, grid);

   for (kk = KBEG; kk <= KEND; kk++){
   for (jj = JBEG; jj <= JEND; jj++){
   for (ii = IBEG; ii <= IEND; ii++){
   for (nv = NVAR; nv--;  ) {
     UU_1[kk][jj][ii][nv] = one_third*(    UU  [kk][jj][ii][nv] +
                                       2.0*UU_1[kk][jj][ii][nv]);
   }}}}

   dt = 2.0/3.0*delta_t;

   for (DIR = 0; DIR < DIMENSIONS; DIR++){

   /*  --------------------------------
        Set normal and tangent indexes  
       -------------------------------- */

     SET_INDEXES(&np_tot, &i, &j, &k,
                 &in,  &in_beg,  &in_end,
                 &it1, &it1_beg, &it1_end,
                 &it2, &it2_beg, &it2_end, grid);

     for (it2 = it2_beg; it2 <= it2_end; it2++) {
     for (it1 = it1_beg; it1 <= it1_end; it1++) {

       for (in = 0; in < np_tot; in++) {
       for (nv = NVAR; nv--;  ) {
         state.v[in][nv] = d->Vc[nv][*k][*j][*i];
       }}

       PRIMTOCON   (state.v, state.u, 0, np_tot - 1);
       RECONSTRUCT (&state, grid);
      
       #ifdef STAGGERED_MHD
        for (in = 0; in < np_tot - 1; in++) {
          state.vL[in][B1] = state.vR[in][B1] = d->Vs[DIR][*k][*j][*i];
        }
       #endif

       RIEMANN (&state, cmax + DIR, grid);
       
       #if NVAR != NFLX
        ADVECT_FLUX (&state, grid);
       #endif
       
       #if ARTIFICIAL_VISCOSITY == YES
        VISC_FLUX   (&state, grid);
       #endif
       
       #ifdef RESISTIVE_MHD
        RESISTIVE_FLUX1 (d->Vc, state.flux, grid, eta_max);
       #endif

       #ifdef STAGGERED_MHD
        EMF_SAVE (state.flux, grid, -1);
       #endif

       GET_RHS (&state, dt, grid);

       for (in = in_beg; in <= in_end; in++) {
       for (nv = NVAR; nv--;  ) {
         UU_1[*k][*j][*i][nv] += state.rhs[in][nv];
       }}
     }}
   }

   #ifdef STAGGERED_MHD
    FCT_UPDATE(3, d, grid);
    AVERAGE_MAGNETIC_FIELD (d->Vs, UU_1, grid);
   #endif

   /* -------------------------
        convert to primitive
      ------------------------- */
     
   DIR = IDIR;
   SET_INDEXES(&np_tot, &i, &j, &k,
               &in,  &in_beg,  &in_end,
               &it1, &it1_beg, &it1_end,
               &it2, &it2_beg, &it2_end, grid);

   for (it2 = KBEG; it2 <= KEND; it2++){
   for (it1 = JBEG; it1 <= JEND; it1++){
     CONTOPRIM (UU_1[*k][*j], state.v, d->flag_x[it2][it1], IBEG, IEND);
     for (in = IBEG; in <= IEND; in++){
     for (nv = 0; nv < NVAR; nv++) {
       d->Vc[nv][*k][*j][*i] = state.v[in][nv];
     }} 
   }}

  #endif
  
}
#undef ALL_SIDES





 void    TempCalc ( real ***V[], Grid *grid, real *SoundSpeed , real *AlfvenSpeed, int first)
 {

   int i, nx;
   static int proc_bottom;
   static MPI_Comm Sub_comm;
   MPI_Comm Cart_comm;

   int remain_dims[] = { 0, 1 };
   int myrank=99;
   int ngh = grid[JDIR].nghost;
   nx = grid[IDIR].np_tot;

   if ( first == 1) {
#ifdef PARALLEL
    AL_Get_cart_comm( SZ, &Cart_comm );
    MPI_Cart_sub (Cart_comm, remain_dims, &Sub_comm);
#endif
   }

   if  (grid[JDIR].is_gbeg == 1) {
   for (i=0 ; i<nx ; i++)  {
    //SoundSpeed[i] = V[PR][0][ngh][i]/V[DN][0][ngh][i];
	  SoundSpeed[i] =  0.01/grid[IDIR].x[i];
    AlfvenSpeed[i] = V[BY][0][ngh][i]/sqrt(V[DN][0][ngh][i]);
   }}



#ifdef PARALLEL
  int rc = MPI_Bcast (&SoundSpeed[0], nx, MPI_DOUBLE, 0, Sub_comm);
  rc = MPI_Bcast (&AlfvenSpeed[0], nx, MPI_DOUBLE, 0, Sub_comm);

#ifdef CHECK
    MPI_Comm_rank (Sub_comm, &myrank);
   printf("myrank , %d , %lf, %lf \n", myrank, SoundSpeed[30], 0.01/grid[IDIR].x[30]);
#endif
#endif
}


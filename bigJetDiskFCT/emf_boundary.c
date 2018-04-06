#include "pluto.h"

/* ------------------------------------------------------------------
           set boundary conditions on electric field 
   ------------------------------------------------------------------ */

#if DIMENSIONS == 3
 #define k0 1
#else
 #define k0 0
#endif
 
#define LOOP_X1_BEG(q)  i = IBEG - 1;                                 \
                        for (k = KBEG - k0; k <= KEND; k++) {         \
                        for (j = JBEG - 1 ; j <= JEND; j++) { q }}

#define LOOP_X1_END(q)  i = IEND + 1;                           \
                        for (k = KBEG - k0; k <= KEND; k++) {   \
                        for (j = JBEG - 1 ; j <= JEND; j++) {  q }}

#define LOOP_X2_BEG(q)  j = JBEG - 1;                           \
                        for (k = KBEG - k0; k <= KEND; k++) {   \
                        for (i = IBEG - 1 ; i <= IEND; i++) {  q }}

#define LOOP_X2_END(q)  j = JEND + 1;                           \
                        for (k = KBEG - k0; k <= KEND; k++) {   \
                        for (i = IBEG - 1 ; i <= IEND; i++) {  q }}

#define LOOP_X3_BEG(q)  k = KBEG - 1;                          \
                        for (j = JBEG - 1; j <= JEND; j++) {   \
                        for (i = IBEG - 1; i <= IEND; i++) {  q }}

#define LOOP_X3_END(q)  k = KEND + 1;                          \
                        for (j = JBEG - 1; j <= JEND; j++) {   \
                        for (i = IBEG - 1; i <= IEND; i++) {  q }}

void EMF_REFLECTIVE_BOUND (EMF *emf, int side, int type);
void EMF_OUTFLOW_BOUND    (EMF *emf, int side);
void EMF_PERIODIC_BOUND   (EMF *emf, int side);

/* ******************************************************************* */
void EMF_BOUNDARY (EMF *emf, Grid *grid)
/*
 *
 * PURPOSE:
 *
 *  Set lower and upper boundary on uu[n], where 
 *  n = 0, ... , NVAR-1. 
 *
 *
 *
 *
 * HISTORY:
 *
 *      Last modified by Andrea Mignone, 25 Oct 2003
 *
 ******************************************************************* */
{
  int   idim;
  int   this_side_type;   /* tells what kind of boundary are applied on this side */
  int   side, iside;      /* lower or upper */
  int   proc_has_side;    /* tells whether the current processor is on the boundary */
  int   i, j, k, nptsx, nptsy;
  static real **u1, **v1;
  int dimx[3] = {1, 0, 0};
  int dimy[3] = {0, 1, 0};
  int dimz[3] = {0, 0, 1};


  if ( (  grid[IDIR].is_gbeg == 1 ) && (grid[JDIR].is_gbeg == 1  ) )  {
     nptsx = NPTSX;
     nptsy = NPTSY;
  }
  else {
     nptsx = 0;
     nptsy = 0;
  }



  if (u1 == NULL){
    u1 = array_2D (NMAX_POINT,NVAR);
    v1 = array_2D (NMAX_POINT,NVAR);
  }

  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);

   AL_Exchange_dim (emf->ezj[0][0], dimx, SZ);
   AL_Exchange_dim (emf->ezi[0][0], dimy, SZ);
   #if DIMENSIONS == 3
    AL_Exchange_dim (emf->exj[0][0], dimz, SZ);
    AL_Exchange_dim (emf->exk[0][0], dimy, SZ);
    AL_Exchange_dim (emf->eyi[0][0], dimz, SZ);
    AL_Exchange_dim (emf->eyk[0][0], dimx, SZ);
   #endif

/*
   AL_Exchange (emf->ezj[0][0], SZ);
   AL_Exchange (emf->ezi[0][0], SZ);
   #if DIMENSIONS == 3
    AL_Exchange (emf->exj[0][0], SZ);
    AL_Exchange (emf->exk[0][0], SZ);
    AL_Exchange (emf->eyi[0][0], SZ);
    AL_Exchange (emf->eyk[0][0], SZ);
   #endif
*/
   MPI_Barrier (MPI_COMM_WORLD);

  #endif


/* INTERNAL BOUNDARY */

   for (k =   KBEG   ; k <= KEND; k++) {         
   for (j = JBEG - 1 ; j < grid[JDIR].lbeg+nptsy ; j++) { 
   for (i =   IBEG   ; i < grid[IDIR].lbeg+nptsx ; i++) { 
  
    if (i == grid[IDIR].lbeg+nptsx-1) {
     emf->ezj[k][j][i] = emf->ezj[k][j][grid[IDIR].lbeg+nptsx];
     emf->ezj[k][j][i] = 0.0;
    } else {
     emf->ezj[k][j][i] = 0.0;
    }

   }}}

   for (k =   KBEG   ; k <= KEND; k++) {
   for (j =   JBEG   ; j < grid[JDIR].lbeg+nptsy ; j++) {
   for (i =  IBEG-1  ; i < grid[IDIR].lbeg+nptsx ; i++) {

    if (j == grid[JDIR].lbeg+nptsy-1) {
     emf->ezi[k][j][i] = emf->ezi[k][grid[JDIR].lbeg+nptsy][i];
     emf->ezi[k][j][i] = 0.0;
    } else {
     emf->ezi[k][j][i] = 0.0;
    }

   }}}
    
/* --------------------------------------------------------
                  Loop on directions       
   -------------------------------------------------------- */

  for (idim = 0; idim < DIMENSIONS; idim++){

    if      (idim == IDIR) side = X1_BEG;
    else if (idim == JDIR) side = X2_BEG;
    else if (idim == KDIR) side = X3_BEG;

  /*  -----------------------------------------
        check whether a given processor 
        abuts the physical boundary 
      ----------------------------------------- */
 
    proc_has_side  = grid[idim].is_gbeg == 1;
    this_side_type = grid[idim].lbound;  

  /* ----------------------------------
       For a given dimension (idim), 
       loop on lower and upper sides    
     ---------------------------------- */

    for (iside = 0; iside < 2; iside++) {

      if (proc_has_side) {

        if (this_side_type == OUTFLOW){

          EMF_OUTFLOW_BOUND (emf, side);            

        }else if (   (this_side_type == REFLECTIVE) 
                  || (this_side_type == AXISYMMETRIC)
                  || (this_side_type == EQTSYMMETRIC) ){

          EMF_REFLECTIVE_BOUND (emf, side, this_side_type);

        }else if (this_side_type == USERDEF) {

          /* -------------------------------------------------
              for userdef boundaries, ghost EMF must be computed
              from the time marching algorithm by extending 
              integration in the corresponding (transverse) 
              boundary zones.
             ------------------------------------------------- */
/*
            EMF_USERDEF_BOUNDARY (emf, side, FACE_CENTER, 
                                  IBEG - 1 , IEND, 
                                  JBEG - 1 , JEND, 
                                  KBEG - k0, KEND, grid);
*/
        }else if (this_side_type == PERIODIC){

          #ifndef PARALLEL   
           EMF_PERIODIC_BOUND (emf, side);
          #endif

        }else{

          print ("Error: Lower boundary %d is not defined. \n", DIR);
          QUIT_PLUTO(1);
        }
      }

  /*  ----  set indexes for upper side  ----  */

      if      (idim == IDIR) side = X1_END;
      else if (idim == JDIR) side = X2_END;
      else if (idim == KDIR) side = X3_END;

      proc_has_side  = grid[idim].is_gend == 1;
      this_side_type = grid[idim].rbound;  
    }
  }
}

/* ************************************************************** */
void EMF_OUTFLOW_BOUND (EMF *emf, int side)
/*
 *
 *
 *
 *
 *
 *
 **************************************************************** */
{
  int nv, i, j, k;

  if (side == X1_BEG) { 

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X1_BEG(D_EXPAND(                      ,
       emf->yl[k][j][i] = emf->yl[k][j][IBEG];  
       emf->yr[k][j][i] = emf->yr[k][j][IBEG];  ,
       emf->zl[k][j][i] = emf->zl[k][j][IBEG];  
       emf->zr[k][j][i] = emf->zr[k][j][IBEG]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X1_BEG(D_EXPAND(                        ,
       emf->ezj[k][j][i] = emf->ezj[k][j][IBEG];  , 
       emf->eyk[k][j][i] = emf->eyk[k][j][IBEG];
     ))

    #endif    

  }else if (side == X1_END){

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X1_END(D_EXPAND(                      ,
       emf->yl[k][j][i] = emf->yl[k][j][IEND];  
       emf->yr[k][j][i] = emf->yr[k][j][IEND];  ,
       emf->zl[k][j][i] = emf->zl[k][j][IEND];  
       emf->zr[k][j][i] = emf->zr[k][j][IEND]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X1_END(D_EXPAND(                        ,                           
       emf->ezj[k][j][i] = emf->ezj[k][j][IEND];  ,                                             
       emf->eyk[k][j][i] = emf->eyk[k][j][IEND];
     ))

    #endif
    
  }else if (side == X2_BEG){

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X2_BEG(D_EXPAND(                      ,
       emf->xl[k][j][i] = emf->xl[k][JBEG][i];  
       emf->xr[k][j][i] = emf->xr[k][JBEG][i];  ,
       emf->zl[k][j][i] = emf->zl[k][JBEG][i];  
       emf->zr[k][j][i] = emf->zr[k][JBEG][i]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X2_BEG(D_EXPAND(                        ,
       emf->ezi[k][j][i] = emf->ezi[k][JBEG][i];  ,                        
       emf->exk[k][j][i] = emf->exk[k][JBEG][i];
     ))

    #endif

  }else if (side == X2_END){

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X2_END(D_EXPAND(                      ,
       emf->xl[k][j][i] = emf->xl[k][JEND][i];  
       emf->xr[k][j][i] = emf->xr[k][JEND][i];  ,
       emf->zl[k][j][i] = emf->zl[k][JEND][i];  
       emf->zr[k][j][i] = emf->zr[k][JEND][i]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X2_END(D_EXPAND(                        ,
       emf->ezi[k][j][i] = emf->ezi[k][JEND][i];  ,
       emf->exk[k][j][i] = emf->exk[k][JEND][i];
     ))
 
    #endif
 
  }else if (side == X3_BEG){
  
    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X3_BEG(                               
       emf->xl[k][j][i] = emf->xl[KBEG][j][i];  
       emf->xr[k][j][i] = emf->xr[KBEG][j][i];  
       emf->yl[k][j][i] = emf->yl[KBEG][j][i];  
       emf->yr[k][j][i] = emf->yr[KBEG][j][i]; 
     )

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X3_BEG( emf->exj[k][j][i] = emf->exj[KBEG][j][i]; 
                  emf->eyi[k][j][i] = emf->eyi[KBEG][j][i]; )
    #endif

  }else if (side == X3_END){

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X3_END(                     
       emf->xl[k][j][i] = emf->xl[KEND][j][i];  
       emf->xr[k][j][i] = emf->xr[KEND][j][i];  
       emf->yl[k][j][i] = emf->yl[KEND][j][i];  
       emf->yr[k][j][i] = emf->yr[KEND][j][i]; 
     )

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X3_END( emf->exj[k][j][i] = emf->exj[KEND][j][i];
                  emf->eyi[k][j][i] = emf->eyi[KEND][j][i]; )
 
    #endif
  }
}

/* *********************************************************** */
void EMF_REFLECTIVE_BOUND (EMF *emf, int side, int type)
/*
 *
 *
 *
 *
 *
 *
 ************************************************************* */
{
  int nv, i, j, k;
  real s = -1.0;

  if (type == EQTSYMMETRIC) s = 1.0;

  if (side == X1_BEG) {   

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X1_BEG(D_EXPAND(                      ,
       emf->yl[k][j][i] = emf->yl[k][j][IBEG];  
       emf->yr[k][j][i] = emf->yr[k][j][IBEG];  ,
       emf->zl[k][j][i] = emf->zl[k][j][IBEG];  
       emf->zr[k][j][i] = emf->zr[k][j][IBEG]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X1_BEG(D_EXPAND(                         ,
       emf->ezj[k][j][i] = s*emf->ezj[k][j][IBEG];  ,
       emf->eyk[k][j][i] = s*emf->eyk[k][j][IBEG]; 
     ))

    #endif   

  }else if (side == X1_END){  

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X1_END(D_EXPAND(                      ,
       emf->yl[k][j][i] = emf->yl[k][j][IEND];  
       emf->yr[k][j][i] = emf->yr[k][j][IEND];  ,
       emf->zl[k][j][i] = emf->zl[k][j][IEND];  
       emf->zr[k][j][i] = emf->zr[k][j][IEND]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X1_END(D_EXPAND(                         ,
       emf->ezj[k][j][i] = s*emf->ezj[k][j][IEND];  ,
       emf->eyk[k][j][i] = s*emf->eyk[k][j][IEND]; 
     ))
  
    #endif
 
  }else if (side == X2_BEG){  

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X2_BEG(D_EXPAND(                      ,
       emf->xl[k][j][i] = emf->xl[k][JBEG][i];  
       emf->xr[k][j][i] = emf->xr[k][JBEG][i];  ,
       emf->zl[k][j][i] = emf->zl[k][JBEG][i];  
       emf->zr[k][j][i] = emf->zr[k][JBEG][i]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X2_BEG(D_EXPAND(                          ,
       emf->ezi[k][j][i] = s*emf->ezi[k][JBEG][i];  ,
       emf->exk[k][j][i] = s*emf->exk[k][JBEG][i];
     ))

    #endif


  }else if (side == X2_END){  

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X2_END(D_EXPAND(                      ,
       emf->xl[k][j][i] = emf->xl[k][JEND][i];  
       emf->xr[k][j][i] = emf->xr[k][JEND][i];  ,
       emf->zl[k][j][i] = emf->zl[k][JEND][i];  
       emf->zr[k][j][i] = emf->zr[k][JEND][i]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X2_END(D_EXPAND(                         ,
       emf->ezi[k][j][i] = s*emf->ezi[k][JEND][i];  ,
       emf->exk[k][j][i] = s*emf->exk[k][JEND][i]; 
     ))

    #endif

  }else if (side == X3_BEG){  

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X3_BEG(                               
       emf->xl[k][j][i] = emf->xl[KBEG][j][i];  
       emf->xr[k][j][i] = emf->xr[KBEG][j][i];  
       emf->yl[k][j][i] = emf->yl[KBEG][j][i];  
       emf->yr[k][j][i] = emf->yr[KBEG][j][i]; 
     )

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X3_BEG( emf->exj[k][j][i] = s*emf->exj[KBEG][j][i]; 
                  emf->eyi[k][j][i] = s*emf->eyi[KBEG][j][i];) 

    #endif

  }else if (side == X3_END){  

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X3_END(                      
       emf->xl[k][j][i] = emf->xl[KEND][j][i];  
       emf->xr[k][j][i] = emf->xr[KEND][j][i];  
       emf->yl[k][j][i] = emf->yl[KEND][j][i];  
       emf->yr[k][j][i] = emf->yr[KEND][j][i]; 
     )

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X3_END( emf->exj[k][j][i] = s*emf->exj[KEND][j][i]; 
                  emf->eyi[k][j][i] = s*emf->eyi[KEND][j][i];) 
    #endif

  }

}

/* *********************************************************** */
void EMF_PERIODIC_BOUND (EMF *emf, int side)
/*
 *
 *
 *
 *
 *
 *
 ************************************************************** */
{
  int nv, i, j, k;

  if (side == X1_BEG) {

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X1_BEG(D_EXPAND(                        ,
       emf->yl[k][j][i] = emf->yl[k][j][i + NX];  
       emf->yr[k][j][i] = emf->yr[k][j][i + NX];  ,
       emf->zl[k][j][i] = emf->zl[k][j][i + NX];  
       emf->zr[k][j][i] = emf->zr[k][j][i + NX]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X1_BEG(D_EXPAND(                          ,
       emf->ezj[k][j][i] = emf->ezj[k][j][i + NX];  ,
       emf->eyk[k][j][i] = emf->eyk[k][j][i + NX];
     ))

    #endif

  }else if (side == X1_END){

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X1_END(D_EXPAND(                        ,
       emf->yl[k][j][i] = emf->yl[k][j][i - NX];  
       emf->yr[k][j][i] = emf->yr[k][j][i - NX];  ,
       emf->zl[k][j][i] = emf->zl[k][j][i - NX];  
       emf->zr[k][j][i] = emf->zr[k][j][i - NX]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X1_END(D_EXPAND(                          ,
       emf->ezj[k][j][i] = emf->ezj[k][j][i - NX];  ,
       emf->eyk[k][j][i] = emf->eyk[k][j][i - NX];  
     ))

    #endif
 
  }else if (side == X2_BEG){

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X2_BEG(D_EXPAND(                        ,
       emf->xl[k][j][i] = emf->xl[k][j + NY][i];  
       emf->xr[k][j][i] = emf->xr[k][j + NY][i];  ,
       emf->zl[k][j][i] = emf->zl[k][j + NY][i];  
       emf->zr[k][j][i] = emf->zr[k][j + NY][i]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X2_BEG(D_EXPAND(                          ,
       emf->ezi[k][j][i] = emf->ezi[k][j + NY][i];  ,
       emf->exk[k][j][i] = emf->exk[k][j + NY][i];
     ))

    #endif
    
  }else if (side == X2_END){

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X2_END(D_EXPAND(                        ,
       emf->xl[k][j][i] = emf->xl[k][j - NY][i];  
       emf->xr[k][j][i] = emf->xr[k][j - NY][i];  ,
       emf->zl[k][j][i] = emf->zl[k][j - NY][i];  
       emf->zr[k][j][i] = emf->zr[k][j - NY][i]; 
     ))

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X2_END(D_EXPAND(                          ,
       emf->ezi[k][j][i] = emf->ezi[k][j - NY][i];  ,
       emf->exk[k][j][i] = emf->exk[k][j - NY][i];
     ))

    #endif

  }else if (side == X3_BEG){
 
    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X3_BEG(                               
       emf->xl[k][j][i] = emf->xl[k + NZ][j][i];  
       emf->xr[k][j][i] = emf->xr[k + NZ][j][i];  
       emf->yl[k][j][i] = emf->yl[k + NZ][j][i];  
       emf->yr[k][j][i] = emf->yr[k + NZ][j][i]; 
     )

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X3_BEG(emf->exj[k][j][i] = emf->exj[k + NZ][j][i]; 
                 emf->eyi[k][j][i] = emf->eyi[k + NZ][j][i];) 

    #endif

  }else if (side == X3_END){

    #if CT_EMF_AVERAGE == HLL_UPWIND

     LOOP_X3_END(                      
       emf->xl[k][j][i] = emf->xl[k - NZ][j][i];  
       emf->xr[k][j][i] = emf->xr[k - NZ][j][i];  
       emf->yl[k][j][i] = emf->yl[k - NZ][j][i];  
       emf->yr[k][j][i] = emf->yr[k - NZ][j][i]; 
     )

    #elif CT_EMF_AVERAGE == ARITHMETIC

     LOOP_X3_END(emf->exj[k][j][i] = emf->exj[k - NZ][j][i];
                 emf->eyi[k][j][i] = emf->eyi[k - NZ][j][i];) 

    #endif
  }
}

#undef LOOP_X1_BEG
#undef LOOP_X1_END
#undef LOOP_X2_BEG
#undef LOOP_X2_END
#undef LOOP_X3_BEG
#undef LOOP_X3_END
#undef k0

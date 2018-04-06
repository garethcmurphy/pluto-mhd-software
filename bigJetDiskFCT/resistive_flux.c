#include "pluto.h"
#include "eta.h"

#define MHD_COOL YES
#define VISC_COOL
//#undef VISC_COOL
#define VISCOSITY YES
//#undef VISCOSITY
//#define DEBUG_VISC
real kvisc_func  (real x1,real x2, real x3, real rho, real etal, real etar); 
/* *********************************************************************** */
void RESISTIVE_FLUX1 (Data_Arr V, real **flux, real **src, Grid *grid, real *eta_max,
							real * SoundSpeed,
							real * AlfvenSpeed)
/* 
 *
 * PURPOSE
 *
 *    Add resistive terms to the energy
 *    and induction equation. It is called in the
 *    during the sweep integrators.
 *
 *    Written by
 *    Matsakos Titos & Andrea Mignone
 *
 ************************************************************************* */
{
	int q=0;
  int i,j,k,nv;
  int loc_beg, nptsx, nptsy;
  real etal[3], etar[3], eta[3];
  real dx_1, dy_1, dz_1;
  real dx1, dx2, dx3;
  real x1, x2, x3;
  real dxBy, dxBz, dyBx, dyBz, dzBx, dzBy;
  real ***Bx, ***By, ***Bz;
  real ***Vr, ***Vz, ***Vphi, ***Rho, ***Tr;
  real vrf,vzf, vphif;
  real scrh;
  real r, dr_1, th;
  real rho;


	real drVr;
	real drVz;
	real drVphi;
	real dzVz;
	real dzVr;
	real dzVphi;

	real divV;

	real tau_rr;
	real tau_rz;
	real tau_rphi;

	real tau_zz;
	real tau_zphi;

	real tau_phiphi;

	real kvisc=1e-6;

  nptsx = NPTSX;
  nptsy = NPTSY;

  Vr = V[VX];
  Vz = V[VY];
  Vphi = V[VZ];
  Bx = V[BX];
  By = V[BY];
  Bz = V[BZ];
  Rho = V[DN];
  Tr = V[TR];

  dzBx = dzBy = dyBx = dyBz = 0.0;
  
  if (DIR == IDIR){     /* ------- X 1     S W E E P  ------------- */
   
    j = *NY_PT;
    k = *NZ_PT;
   
    dy_1 = 1.0/grid[JDIR].dx[j];
    dz_1 = 1.0/grid[KDIR].dx[k];

    for (i = IBEG - 1; i <= IEND; i++){
     for (nv = 0; nv < NVAR; nv++) src[i][nv] = 0.0;
    }

    loc_beg = IBEG-1;
    if   (  grid[IDIR].is_gbeg == 1 ) loc_beg = IBEG;


 if ( (  grid[IDIR].is_gbeg == 1 ) && (grid[JDIR].is_gbeg == 1  ) )  {
 
    if (j < grid[JDIR].lbeg+nptsy) loc_beg = grid[IDIR].lbeg+nptsx-1;
}

    for (i = loc_beg; i <= IEND; i++){
      
      #if GEOMETRY == CARTESIAN
      
       dx_1 = 1.0/grid[IDIR].dx[i];
      
       EXPAND(                                               ,
              dxBy = (By[k][j][i + 1] - By[k][j][i])*dx_1;   ,
              dxBz = (Bz[k][j][i + 1] - Bz[k][j][i])*dx_1; )

       #if DIMENSIONS >= 2
        dyBx = 0.25*(  Bx[k][j + 1][i]     - Bx[k][j - 1][i]
	             + Bx[k][j + 1][i + 1] - Bx[k][j - 1][i + 1])*dy_1;
        #if DIMENSIONS == 3
         dzBx = 0.25*(  Bx[k + 1][j][i]     - Bx[k - 1][j][i]
	              + Bx[k + 1][j][i + 1] - Bx[k - 1][j][i + 1])*dz_1;
        #endif
       #endif

      #elif GEOMETRY == POLAR 

       r  = grid[IDIR].x[i];
       dr_1 = 1.0/grid[IDIR].dx[i];

       EXPAND(                                                                                                ,
              dxBy = 0.5*(By[k][j][i + 1] + By[k][j][i])/(r+0.5/dr_1) + (By[k][j][i + 1] - By[k][j][i])*dr_1; ,
              dxBz = (Bz[k][j][i + 1] - Bz[k][j][i])*dr_1;                                                    )


       #if DIMENSIONS >= 2
        dyBx = 0.25*( Bx[k][j + 1][i] - Bx[k][j - 1][i] 
                    + Bx[k][j + 1][i + 1] - Bx[k][j - 1][i + 1])*dy_1/(r+0.5/dr_1);
        #if DIMENSIONS == 3
         dzBx = 0.25*(  Bx[k + 1][j][i]     - Bx[k - 1][j][i]
                      + Bx[k + 1][j][i + 1] - Bx[k - 1][j][i + 1])*dz_1;
        #endif
       #endif

      #elif GEOMETRY == CYLINDRICAL

       r  = grid[IDIR].x[i];
       dr_1 = 1.0/grid[IDIR].dx[i];

       EXPAND(                                                                                                ,
              dxBy = (By[k][j][i + 1] - By[k][j][i])*dr_1;                                                    ,
              dxBz = 0.5*(Bz[k][j][i + 1] + Bz[k][j][i])/(r+0.5/dr_1) + (Bz[k][j][i + 1] - Bz[k][j][i])*dr_1; )

       #if DIMENSIONS >= 2
        dyBx = 0.25*(  Bx[k][j + 1][i]     - Bx[k][j - 1][i]
 	            + Bx[k][j + 1][i + 1] - Bx[k][j - 1][i + 1])*dy_1;
       #endif

      #elif GEOMETRY == SPHERICAL

       r  = grid[IDIR].x[i];
       dr_1 = 1.0/grid[IDIR].dx[i];

       EXPAND(                                                                                                ,
              dxBy = 0.5*(By[k][j][i + 1] + By[k][j][i])/(r+0.5/dr_1) + (By[k][j][i + 1] - By[k][j][i])*dr_1; ,
              dxBz = 0.5*(Bz[k][j][i + 1] + Bz[k][j][i])/(r+0.5/dr_1) + (Bz[k][j][i + 1] - Bz[k][j][i])*dr_1; )


       #if DIMENSIONS >= 2
        dyBx = 0.25*( Bx[k][j + 1][i] - Bx[k][j - 1][i]
                    + Bx[k][j + 1][i + 1] - Bx[k][j - 1][i + 1])*dy_1/(r+0.5/dr_1);
        #if DIMENSIONS == 3
         th =  grid[JDIR].x[j];

         dzBx = 0.25*(  Bx[k + 1][j][i]     - Bx[k - 1][j][i]
                      + Bx[k + 1][j][i + 1] - Bx[k - 1][j][i + 1])*dz_1/(r+0.5/dr_1)/sin(th);
        #endif
       #endif
      #endif  /* -- end #if GEOMETRY -- */

      /* --------------------------------------------------------------------
                           compute resistivity 
         -------------------------------------------------------------------- */

      x1  = grid[IDIR].x[i];
      x2  = grid[JDIR].x[j];
      x3  = grid[KDIR].x[k];
      dx1 = grid[IDIR].dx[i];
      dx2 = grid[JDIR].dx[j];
      dx3 = grid[KDIR].dx[k];

/*
 *
 * */
      //eta_func(x1 + 0.5*dx1, x2, x3, dyBz - dzBy, dzBx - dxBz, dxBy - dyBx, eta);
      eta_func_timedep(  x1  ,x2,x3, i ,SoundSpeed,AlfvenSpeed,etal);
      eta_func_timedep(x1+dx1,x2,x3,i+1,SoundSpeed,AlfvenSpeed,etar);
		kvisc=2.0/3.0*0.5*(Rho[k][j][i]*Tr[k][j][i]* *etal + Rho[k][j][i+1]*Tr[k][j][i+1]* *etar);

		for (q=0; q<3; ++q) { eta[q] = 0.5*(etal[q] + etar[q]); }
      /* --------------------------------------------------------------------
                             compute fluxes
         -------------------------------------------------------------------- */

      EXPAND(                                       ,
             flux[i][BY] += -eta[2]*(dxBy - dyBx);  ,
             flux[i][BZ] +=  eta[1]*(dzBx - dxBz);)
      
      scrh = EXPAND(0.0                                                     ,
                    - (By[k][j][i] + By[k][j][i + 1])*(dxBy - dyBx)*eta[2]  ,
                    + (Bz[k][j][i] + Bz[k][j][i + 1])*(dzBx - dxBz)*eta[1]);
      
      flux[i][EN] += 0.5*scrh;

#ifdef MHD_COOL
		// Calc cell centred values 
		// Calc currents
		// Calc Eta
      dzBx = 0.0;
		dzBy = 0.0;
		dxBy = (By[k][j][i+1] - By[k][j][i-1])/(grid[IDIR].x[i+1]-grid[IDIR].x[i-1]);
		dyBx = (Bx[k][j+1][i] - Bx[k][j-1][i])/(grid[JDIR].x[j+1]-grid[JDIR].x[j-1]);
      dxBz = (grid[IDIR].x[i+1]*Bz[k][j][i+1] - grid[IDIR].x[i-1]*Bz[k][j][i-1])/(grid[IDIR].x[i+1]-grid[IDIR].x[i-1])/x1;
		dyBz = (Bz[k][j+1][i] - Bz[k][j-1][i])/(grid[JDIR].x[j+1]-grid[JDIR].x[j-1]);
      //eta_func(x1 , x2, x3, dyBz - dzBy, dzBx - dxBz, dxBy - dyBx, eta);
		src[i][EN] -= (etal[0]*dyBz*dyBz+etal[1]*dxBz*dxBz+etal[2]*(dxBy-dyBx)*(dxBy-dyBx));
#endif


#ifdef VISCOSITY



/* x sweep*/
       r  = grid[IDIR].x[i];
       dr_1 = 1.0/grid[IDIR].dx[i];
    dy_1 = 1.0/grid[JDIR].dx[j];

/*
 * r derivatives on r face
 * */
	drVr  =(   Vr[k][j][i+1] -    Vr[k][j][i] ) * dr_1;
	drVz  =(   Vz[k][j][i+1] -    Vz[k][j][i] ) * dr_1;
	drVphi=( Vphi[k][j][i+1] -  Vphi[k][j][i] ) * dr_1;


/*
 * z derivatives on r face
 * */
	dzVr  = 0.25*((   Vr[k][j+1][i+1] +    Vr[k][j+1][i] ) - (   Vr[k][j-1][i+1] +    Vr[k][j-1][i] ))*dy_1; 
	dzVz  = 0.25*((   Vz[k][j+1][i+1] +    Vz[k][j+1][i] ) - (   Vz[k][j-1][i+1] +    Vz[k][j-1][i] ))*dy_1; 
	dzVphi= 0.25*(( Vphi[k][j+1][i+1] +  Vphi[k][j+1][i] ) - ( Vphi[k][j-1][i+1] +  Vphi[k][j-1][i] ))*dy_1; 


/* Face values velocity*/
			vrf = 0.5*( Vr[k][j][i]  + Vr[k][j][i+1] );
			vzf = 0.5*( Vz[k][j][i]  + Vz[k][j][i+1] );
			vphif = 0.5*( Vphi[k][j][i]  + Vphi[k][j][i+1] );
/*
 * Divergence on r face
 * */

	divV= ( (r +1./dr_1)*  Vr[k][j][i+1] - (r  )*(Vr[k][j][i] )) *dr_1 /(r+0.5/dr_1)   +dzVz;
	
	tau_rr   = kvisc * (2. *( drVr) - 2./3. * divV); 
	tau_rz   = kvisc * ( drVz + dzVr); 
	tau_rphi = kvisc * ( drVphi  -   vphif /(r +0.5/dr_1) ); /* + SOURCE TERM ABSORBED INTO tau_rphi*/

/*
	if (j==200) {
	printf("\ni=%d j=%d, %lf %lf %lf ", i, j,  drVphi,   Vphi[k][j+1][i] ,  Vphi[k][j][i]);
	}
	*/
	//tau_zz   = kvisc * ( 2. *( dzVz) -2./3.* divV ); 
	//tau_zphi = kvisc * ( dzVphi ); 

	//tau_phiphi = kvisc * (2. * ( Vr /r ) -2./3.*divV ); 


#ifdef DEBUG_VISC
		if (isnan(tau_rr))
		{
			printf("\n i,j = %d,%d",  i,j);
			printf("\nr = %lf",  r);
			printf("\ndx_1 = %lf",  dx_1);
			printf("\ndr_1 = %lf",  dr_1);
			printf("\ndy_1 = %lf",  dy_1);
			printf("\nVr = %lf",  Vr[k][j][i+1]  );
			printf("\nVz [%d][%d][%d]= %lf",k,j+1,i+1,  Vz[k][j+1][i+1]  );
			printf("\nVz = %lf",  Vz[k][j+1][i  ]  );
			printf("\nVz = %lf",  Vz[k][j-1][i+1]  );
			printf("\nVz = %lf",  Vz[k][j-1][i  ]  );
			printf("\ndzVz = %lf", dzVz);
			printf("\ndivV = %lf", divV);
			printf("\nCaught NAN tau_rr at line %d\n", __LINE__);
			exit(0);
			}
#endif 
             flux[i][VX] -=  tau_rr ;
             flux[i][VY] -=  tau_rz ;
             flux[i][VZ] -=  tau_rphi ;

             flux[i][EN] -=  vrf*tau_rr + vzf*tau_rz + vphif *tau_rphi  ;
			/*
			 * Compute source term for energy
			 */ 
/*
 * div V on cell centre
 * */
	divV= 0.5*( (r +1./dr_1)*  Vr[k][j][i+1] - (r-1./dr_1  )*(Vr[k][j][i-1] )) *dr_1 /r   +dzVz;
      eta_func_timedep(  x1  ,x2,x3, i ,SoundSpeed,AlfvenSpeed,etal);

		kvisc=2.0/3.0*(Rho[k][j][i]*Tr[k][j][i]* *etal);
	tau_phiphi = kvisc * (2. * ( Vr[k][j][i] /r ) -2./3.*divV ); 
             src [i][VX] +=  - tau_phiphi/r   ;


#ifdef VISC_COOL
		// Calc cell centred values 

/* x sweep*/
       r  = grid[IDIR].x[i];
       dr_1 = 1.0/grid[IDIR].dx[i];
    dy_1 = 1.0/grid[JDIR].dx[j];

/*
 * r derivatives on cell centre
 * */
	drVr  =0.5*(   Vr[k][j][i+1] -    Vr[k][j][i-1] ) * dr_1;
	drVz  =0.5*(   Vz[k][j][i+1] -    Vz[k][j][i-1] ) * dr_1;
	drVphi=0.5*( Vphi[k][j][i+1] -  Vphi[k][j][i-1] ) * dr_1;


/*
 * z derivatives on cell centre
 * */
	dzVr  = 0.5*(  Vr[k][j+1][i]  -    Vr[k][j-1][i] )*dy_1; 
	dzVz  = 0.5*(  Vz[k][j+1][i]  -    Vz[k][j-1][i] )*dy_1; 
	dzVphi= 0.5*(Vphi[k][j+1][i]  -  Vphi[k][j-1][i] )*dy_1; 

/*
 * div V on cell centre
 * */
	divV= 0.5*( (r +1./dr_1)*  Vr[k][j][i+1] - (r -1./dr_1 )*(Vr[k][j][i-1] )) *dr_1 /r   +dzVz;

/* kvisc on cell centre */
      eta_func_timedep(  x1  ,x2,x3, i ,SoundSpeed,AlfvenSpeed,etal);
		kvisc=2.0/3.0*(Rho[k][j][i]*Tr[k][j][i]* *etal);
	tau_rr   =  (2. *( drVr) - 2./3. * divV); 
	tau_rz   =  ( drVz + dzVr); 
	tau_rphi =  ( drVphi -  Vphi[k][j][i]/r ); /* + SOURCE TERM ABSORBED INTO tau_rphi*/

	tau_zz   =  ( 2. *( dzVz) -2./3.* divV ); 
	tau_zphi =  ( dzVphi ); 

	tau_phiphi =  (2. * ( Vr[k][j][i] /r ) -2./3.*divV ); 

		src[i][EN] -= kvisc*(0.5*( tau_rr* tau_rr +   tau_zz * tau_zz +  tau_phiphi* tau_phiphi)+ 
		 tau_rz* tau_rz +  tau_rphi * tau_rphi + tau_zphi *tau_zphi) ;
		
		if (isnan( src[i][EN]))
		{
			printf("\nnan");
		}
#endif
		
#ifdef DEBUG_VISC
				if (i==200 && j==200)
				{
					printf("\n flux %e" ,  flux[i][VZ]);
					printf(" drVphi %e %e %e" ,  drVphi,  Vphi[k][j][i] , Vphi[k][j][i+1] );
					printf(" tau_rphi %e" ,   tau_rphi);
				}
#endif 



#endif

      *eta_max = dmax(*eta_max,eta[0]);
      *eta_max = dmax(*eta_max,eta[1]);
      *eta_max = dmax(*eta_max,eta[2]);
	   *eta_max = dmax(*eta_max,kvisc / V[DN][k][j][i] );
    }    
    
  }else if (DIR == JDIR){    /* ------- X 2     S W E E P  ------------- */

    i = *NX_PT;
    k = *NZ_PT;
   
    dx_1 = 1.0/grid[IDIR].dx[i];
    dz_1 = 1.0/grid[KDIR].dx[k];

    for (j = JBEG - 1; j <= JEND; j++){
     for (nv = 0; nv < NVAR; nv++) src[j][nv] = 0.0;
    }

    loc_beg = IBEG-1;
 if ( (  grid[IDIR].is_gbeg == 1 ) && (grid[JDIR].is_gbeg == 1  ) )  {
  
    if (i < grid[IDIR].lbeg+nptsx) loc_beg = grid[JDIR].lbeg+nptsy-1;
}

    for (j = loc_beg; j <= JEND; j++){
   
      #if GEOMETRY == CARTESIAN
  
       dy_1 = 1.0/grid[JDIR].dx[j];
      
       EXPAND(dyBx = (Bx[k][j + 1][i] - Bx[k][j][i])*dy_1;                   ,
              dxBy = 0.25*(  By[k][j][i + 1]     - By[k][j][i - 1]
                           + By[k][j + 1][i + 1] - By[k][j + 1][i - 1])*dx_1;  ,
              dyBz = (Bz[k][j + 1][i] - Bz[k][j][i])*dy_1;)
       #if DIMENSIONS == 3
        dzBy = 0.25*(  By[k + 1][j][i]     - By[k - 1][j][i]
                     + By[k + 1][j + 1][i] - By[k - 1][j + 1][i])*dz_1;
       #endif

      #elif GEOMETRY == CYLINDRICAL
   
       dy_1 = 1.0/grid[JDIR].dx[j];
      
       EXPAND(dyBx = (Bx[k][j + 1][i] - Bx[k][j][i])*dy_1;                   ,
              dxBy = 0.25*(  By[k][j][i + 1]     - By[k][j][i - 1]
                           + By[k][j + 1][i + 1] - By[k][j + 1][i - 1])*dx_1;  ,
              dyBz = (Bz[k][j + 1][i] - Bz[k][j][i])*dy_1;)

      #elif GEOMETRY == POLAR 

       r    = grid[IDIR].x[i];
       dr_1 = 1.0/grid[IDIR].dx[i];
       dy_1 = 1.0/grid[JDIR].dx[j];

       EXPAND( dyBx = (Bx[k][j + 1][i] - Bx[k][j][i])*dy_1/r;                                                    ,
               dxBy = 0.25*( (r + 1.0/dr_1)*By[k][j][i + 1]     - (r - 1.0/dr_1)*By[k][j][i - 1]
                           + (r + 1.0/dr_1)*By[k][j + 1][i + 1] - (r - 1.0/dr_1)*By[k][j + 1][i - 1] )*dx_1/r;   ,
               dyBz = (Bz[k][j + 1][i] - Bz[k][j][i])*dy_1/r;                                                      )

       #if DIMENSIONS == 3
        dzBy = 0.25*(  By[k + 1][j][i]     - By[k - 1][j][i]
                     + By[k + 1][j + 1][i] - By[k - 1][j + 1][i])*dz_1;
       #endif

      #elif GEOMETRY == SPHERICAL

       r  = grid[IDIR].x[i];
       dr_1 = 1.0/grid[IDIR].dx[i];
       th  = grid[JDIR].x[j];
       dy_1 = 1.0/grid[JDIR].dx[j];

       EXPAND( dyBx = (Bx[k][j + 1][i] - Bx[k][j][i])*dy_1/r;                                                          ,
               dxBy = 0.25*( (r + 1.0/dr_1)*By[k][j][i + 1]     - (r - 1.0/dr_1)*By[k][j][i - 1]
                           + (r + 1.0/dr_1)*By[k][j + 1][i + 1] - (r - 1.0/dr_1)*By[k][j + 1][i - 1] )*dx_1/r;         ,
               dyBz = (sin(th + 1.0/dy_1)*Bz[k][j + 1][i] - sin(th)*Bz[k][j][i])*dy_1/r/sin(th+0.5/dy_1);  )

       #if DIMENSIONS == 3
        dzBy = 0.25*(  By[k + 1][j][i]     - By[k - 1][j][i]
                     + By[k + 1][j + 1][i] - By[k - 1][j + 1][i])*dz_1/(r*sin(th+0.5/dy_1));
       #endif
  
      #endif

      /* --------------------------------------------------------------------
                       compute resistivity 
         -------------------------------------------------------------------- */
    
      x1 = grid[IDIR].x[i];
      x2 = grid[JDIR].x[j];
      x3 = grid[KDIR].x[k];
      dx1 = grid[IDIR].dx[i];
      dx2 = grid[JDIR].dx[j];
      dx3 = grid[KDIR].dx[k];
      
      eta_func_timedep(x1,  x2  ,x3,i,SoundSpeed,AlfvenSpeed,etal);
      eta_func_timedep(x1,x2+dx2,x3,i,SoundSpeed,AlfvenSpeed,etar);
		kvisc= 2./3.*0.5*(*etal*Rho[k][j][i]*Tr[k][j][i] + *etar*Rho[k][j+1][i]*Tr[k][j+1][i]);

		for (q=0; q<3; ++q) { eta[q] = 0.5*(etal[q] + etar[q]); }

      /* --------------------------------------------------------------------
                          compute fluxes
         -------------------------------------------------------------------- */

      EXPAND(flux[j][BX] +=  eta[2]*(dxBy - dyBx);   ,
                                                     ,
             flux[j][BZ] += -eta[0]*(dyBz - dzBy);)
      
      scrh = EXPAND( + (Bx[k][j + 1][i] + Bx[k][j][i])*(dxBy - dyBx)*eta[2]  ,
                                                                            ,
                     - (Bz[k][j + 1][i] + Bz[k][j][i])*(dyBz - dzBy)*eta[0]);
      
      flux[j][EN] += 0.5*scrh;
      
#ifdef VISCOSITY
       r  = grid[IDIR].x[i];
       dr_1 = 1.0/grid[IDIR].dx[i];
    dy_1 = 1.0/grid[JDIR].dx[j];

/* y sweep*/
	drVr  = 0.25*((   Vr[k][j+1][i+1] +    Vr[k][j  ][i+1] ) - (   Vr[k][j+1][i-1] +    Vr[k][j][i-1] ))*dr_1; 
	drVz  = 0.25*((   Vz[k][j+1][i+1] +    Vz[k][j  ][i+1] ) - (   Vz[k][j+1][i-1] +    Vz[k][j][i-1] ))*dr_1; 
	drVphi= 0.25*(( Vphi[k][j+1][i+1] +  Vphi[k][j  ][i+1] ) - ( Vphi[k][j+1][i-1] +  Vphi[k][j][i-1] ))*dr_1; 

	dzVr  =(   Vr[k][j+1][i] -    Vr[k][j][i] ) * dy_1;
	dzVz  =(   Vz[k][j+1][i] -    Vz[k][j][i] ) * dy_1;
	dzVphi=( Vphi[k][j+1][i] -  Vphi[k][j][i] ) * dy_1;

/*
	if (j==200) {
	printf("\ni=%d j=%d, %lf %lf %lf ", i, j,  dzVphi,   Vphi[k][j+1][i] ,  Vphi[k][j][i]);
	}
	*/



	divV= ( (r +1./dr_1)*  (Vr[k][j][i+1] + Vr[k][j+1][i+1] ) - (r -1./dr_1)*(Vr[k][j][i-1]+Vr[k][j+1][i-1] ))*0.25 *dr_1 /r  +dzVz;
	
//	tau_rr   = kvisc * (2. *( drVr) - 2./3. * divV); 
	tau_rz   = kvisc * ( drVz + dzVr); 
//	tau_rphi = kvisc * ( drVphi  - Vphi/r); 

	tau_zz   = kvisc * ( 2. *( dzVz) -2./3.* divV ); 
	tau_zphi = kvisc * ( dzVphi ); 

	//tau_phiphi = kvisc * (2. * ( Vr /r ) -2./3.*divV ); 
      *eta_max = dmax(*eta_max,eta[0]);
      *eta_max = dmax(*eta_max,eta[1]);
      *eta_max = dmax(*eta_max,eta[2]);
	   *eta_max = dmax(*eta_max,kvisc / V[DN][k][j][i] );

#ifdef DEBUG_VISC
		if ( isnan(tau_zz))
		{
			printf("\ni,j = %d,%d",  i,j);
			printf("\nkvisc = %lf",  kvisc);
			printf("\nr = %lf",  r);
			printf("\ndx_1 = %lf",  dx_1);
			printf("\ndr_1 = %lf",  dr_1);
			printf("\ndy_1 = %lf",  dy_1);
			printf("\nVr = %lf",  Vr[k][j][i+1]  );
			printf("\nVz [%d][%d][%d]= %lf",k,j+1,i+1,  Vz[k][j+1][i+1]  );
			printf("\nVz = %lf",  Vz[k][j+1][i  ]  );
			printf("\nVz = %lf",  Vz[k][j-1][i+1]  );
			printf("\nVz = %lf",  Vz[k][j-1][i  ]  );
			printf("\ndzVz = %lf", dzVz);
			printf("\ndivV = %lf", divV);
			printf("\nCaught NAN tau_zz at line %d\n", __LINE__);
			exit(0);
			}
#endif 
			vrf   = 0.5*(   Vr[k][j][i]  +   Vr[k][j+1][i] );
			vzf   = 0.5*(   Vz[k][j][i]  +   Vz[k][j+1][i] );
			vphif = 0.5*( Vphi[k][j][i]  + Vphi[k][j+1][i] );

#ifdef DEBUG_VISC
				if (i==200 && j==200)
				{
					printf("\n flux %e" ,  flux[i][VZ]);
					printf(" dzVphi %e" ,  dzVphi);
					printf(" tau_zphi %e" ,   tau_zphi);
				}
#endif 

             flux[j][VX] -=  tau_rz ;
             flux[j][VY] -=  tau_zz ;
             flux[j][VZ] -=  tau_zphi ;
             flux[j][EN] -=  vrf*tau_rz + vzf*tau_zz + vphif *tau_zphi  ;

#ifdef DEBUG_VISC
				if (i==200 && j==200)
				{
					printf("\n flux %e" ,  flux[i][VZ]);
					printf(" dzVphi %e" ,  dzVphi);
					printf(" etal %e" ,  etal[0]);
					printf(" tau_zphi %e" ,   tau_zphi);
				}
#endif 


#endif
    }      

  }else{                 /* ------- X 3     S W E E P  ------------- */
   
    i = *NX_PT;
    j = *NY_PT;
   
    dx_1 = 1.0/grid[IDIR].dx[i];
    dy_1 = 1.0/grid[JDIR].dx[j];

    for (k = KBEG - 1; k <= KEND; k++){
   
      dz_1 = 1.0/grid[KDIR].dx[k];
  
      #if GEOMETRY == CARTESIAN
      
       dxBz = 0.25*(  Bz[k][j][i + 1]     - Bz[k][j][i - 1]
                    + Bz[k + 1][j][i + 1] - Bz[k + 1][j][i - 1])*dx_1;
       dyBz = 0.25*(  Bz[k][j + 1][i]     - Bz[k][j - 1][i]
                    + Bz[k + 1][j + 1][i] - Bz[k + 1][j - 1][i])*dy_1;
       dzBx = (Bx[k + 1][j][i] - Bx[k][j][i])*dz_1;
       dzBy = (By[k + 1][j][i] - By[k][j][i])*dz_1;

      #elif GEOMETRY == POLAR
      
       r  = grid[IDIR].x[i];
       dr_1 = 1.0/grid[IDIR].dx[i];
      
       dxBz = 0.25*( Bz[k][j][i + 1]     - Bz[k][j][i - 1]
                   + Bz[k + 1][j][i + 1] - Bz[k + 1][j][i - 1])*dx_1;
       dyBz = 0.25*(  Bz[k][j + 1][i]     - Bz[k][j - 1][i]
                    + Bz[k + 1][j + 1][i] - Bz[k + 1][j - 1][i])*dy_1/r;
       dzBx = (Bx[k + 1][j][i] - Bx[k][j][i])*dz_1;
       dzBy = (By[k + 1][j][i] - By[k][j][i])*dz_1;
      
      #elif GEOMETRY == SPHERICAL

       r    = grid[IDIR].x[i];
       dr_1 = 1.0/grid[IDIR].dx[i];
       th   =  grid[JDIR].x[j];
       dy_1 = 1.0/grid[JDIR].dx[j];

       dxBz = 0.25*( (r + 1.0/dr_1)*Bz[k][j][i + 1]     - (r - 1.0/dr_1)*Bz[k][j][i - 1]
                   + (r + 1.0/dr_1)*Bz[k + 1][j][i + 1] - (r - 1.0/dr_1)*Bz[k + 1][j][i - 1])*dx_1/r;
       dyBz = 0.25*(  sin(th + 1.0/dy_1)*Bz[k][j + 1][i]     - sin(th - 1.0/dy_1)*Bz[k][j - 1][i]
                    + sin(th + 1.0/dy_1)*Bz[k + 1][j + 1][i] - sin(th - 1.0/dy_1)*Bz[k + 1][j - 1][i])*dy_1/r/sin(th);
       dzBx = (Bx[k + 1][j][i] - Bx[k][j][i])*dz_1/r/sin(th);
       dzBy = (By[k + 1][j][i] - By[k][j][i])*dz_1/r/sin(th);

      #endif









      /* --------------------------------------------------------------------
                       compute resistivity 
         -------------------------------------------------------------------- */
      
      x1 = grid[IDIR].x[i];
      x2 = grid[JDIR].x[j];
      x3 = grid[KDIR].x[k];
      dx1 = grid[IDIR].dx[i];
      dx2 = grid[JDIR].dx[j];
      dx3 = grid[KDIR].dx[k];

      //eta_func(x1, x2, x3 + 0.5*dx3, dyBz - dzBy, dzBx - dxBz, dxBy - dyBx, eta);
      eta_func_timedep(x1,x2,  x3  ,i,SoundSpeed,AlfvenSpeed,etal);
      eta_func_timedep(x1,x2,x3+dx3,i,SoundSpeed,AlfvenSpeed,etar);
		for (q=0; q<3; ++q) { eta[q] = 0.5*(etal[q] + etar[q]); }
      
      /* --------------------------------------------------------------------
                          compute fluxes
         -------------------------------------------------------------------- */

      flux[k][BX] += -eta[1]*(dzBx - dxBz);
      flux[k][BY] +=  eta[0]*(dyBz - dzBy);

      scrh = + (By[k + 1][j][i] + By[k][j][i])*(dyBz - dzBy)*eta[0]
             - (Bx[k + 1][j][i] + Bx[k][j][i])*(dzBx - dxBz)*eta[1];

      flux[k][EN] += 0.5*scrh;

      *eta_max = dmax(*eta_max,eta[0]);
      *eta_max = dmax(*eta_max,eta[1]);
      *eta_max = dmax(*eta_max,eta[2]);
    }      
  }
}

real kvisc_func  (real x1,real x2,real x3 ,real rho, real etal, real etar) 
{
		real kvisc=0.0;

		kvisc=2.0/3.0/(1.e-1)*rho*0.5* (etal+etar);
		/* Kulzniak Kita */

	return kvisc;
}

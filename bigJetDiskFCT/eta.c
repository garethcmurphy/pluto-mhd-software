#include "pluto.h"

/*-------------------------------------------------------------------------------------------*/
/*-----x,y,z do not refer to cartesian coordinates but to the 1st,2nd and 3rd coordinate-----*/
/*-------------------------------------------------------------------------------------------*/

void eta_func_timedep(double x, double y, double z, int i, real * SoundSpeed , real * AlfvenSpeed , real eta_ar[])
{


 real alphap = 0.9;
 real alphat = alphap;
 real eps = 0.1;
 real b0d = eps*sqrt(2.*0.3); 
 real scrh1, scrh2;
 real height;

 height = sqrt(SoundSpeed[i]*x*x*x);
 scrh1 = y/height;
 //scrh2 = b0d/sqrt(x)*height*exp(-1.*scrh1*scrh1);
 scrh2 = AlfvenSpeed[i]*height*exp(-1.*scrh1*scrh1);
 // Sound speed for Purely hydro
 scrh2 = SoundSpeed[i]*height*exp(-1.*scrh1*scrh1);
 // Kluzniak Kita
 real r=x;
 real R=sqrt(x*x+y*y);
 if (R > 0 && r>0 ) scrh2 = SoundSpeed[i] + (gmm-1)/gmm * (1/R -1/r) ;
 if (x>0) scrh2 = scrh2 * sqrt (x*x*x);
 scrh2 = dmax(scrh2,0);

// scrh2 = SoundSpeed[i]*height*exp(-1.*scrh1*scrh1);
/*
  if ( y > height){
	  scrh2=0;
  }
  */
 //if (x> 10) scrh2=0;
 eta_ar[0] = alphat*scrh2;       
 eta_ar[1] = alphat*scrh2;       
 eta_ar[2] = alphap*scrh2;   



}

void eta_func(double x, double y, double z, double Jx, double Jy, double Jz, double eta_ar[])
{

 real alphap = 1.0;
 real alphat = 1.0;
 real eps = 0.1;
 real b0d = eps*sqrt(2.*0.3); 
 real scrh1, scrh2;

 scrh1 = y/eps/x;
 scrh2 = b0d*eps*sqrt(x)*exp(-1.*scrh1*scrh1);

 eta_ar[0] = alphat*scrh2;       
 eta_ar[1] = alphat*scrh2;       
 eta_ar[2] = alphap*scrh2;   

/*-----EXAMPLES---------------------------------*/
/*-----1)-Current triggered (and dependent)-----*/

/*
eta[0] = 0.01 + 0.001*(max(2.0,Jz))
*/

/*-----2)-Space dependent-----------------------*/

/*
eta[0] = 0.01*exp(-(20.0*x)*(20.0*x)-(10.0*y)*(10.0*y));
*/

}



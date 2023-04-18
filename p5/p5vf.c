#include "p5.h"

////The Verify algorithm
int Verify(unsigned char *ga, fq_t *y, fq_t *z, pubpar *par,fq_mat_t L0)
{
   fq_poly_t phi;
   fq_poly_init(phi,par->Fp);   
   // IntPoly(phi, y, par);
      
   fq_poly_t psi;
   fq_poly_init(psi,par->Fp);   
   // IntPoly(psi,z, par);
  
   int b;
   
   fq_t zl, zr, zo;
   fq_init(zl,par->Fp);
   fq_init(zr,par->Fp);
   fq_init(zo,par->Fp);
   fq_zero(zo,par->Fp);
   
   // fq_poly_evaluate_fq(zl,phi,zo,par->Fp);
   // fq_poly_evaluate_fq(zr,psi,zo,par->Fp);

   // interpolate y=F(0)
   fq_t temp;
	fq_init(temp,par->Fp);

	fq_zero(zl,par->Fp);
   fq_zero(zr,par->Fp);
	for(int j = 0; j < par->k; j++){
		fq_mul(temp,y[j],fq_mat_entry(L0,j,0),par->Fp);
		fq_add(zl,zl,temp,par->Fp);
      fq_mul(temp,z[j],fq_mat_entry(L0,j,0),par->Fp);
      fq_add(zr,zr,temp,par->Fp);
	}
   
////compare between two results   
   fmpz_t Zf;
   fmpz_init(Zf);
   unsigned char L[32];
   unsigned char R[32];
   unsigned char Temp[32];
   
   fq2fmpz(Zf,zl,par->Fp);
   fmpz2chars(Temp,Zf);
   int ms=crypto_scalarmult_ristretto255(L,Temp,ga);
     
   fq2fmpz(Zf,zr,par->Fp);
   fmpz2chars(Temp,Zf);
   ms=crypto_scalarmult_ristretto255_base(R,Temp);
    
 
   b=sodium_memcmp(L, R, 32);
   return b;   
}




#include "p5.h"

////This is the problem generation algorithm
void ProbGen(unsigned char *ga, fq_mat_t *c, fq_t *b, fq_mat_t X, pubpar *par)
{
 flint_rand_t state;
 flint_randinit(state);
 
 //choose the alpha
 fq_t alpha;
 fq_init(alpha,par->Fp);
 fq_randtest(alpha,state,par->Fp);
 while (fq_is_zero(alpha,par->Fp))
       fq_randtest(alpha,state,par->Fp);

//calculate the first element
 fmpz_t zalpha;
 fmpz_init(zalpha);
 fq2fmpz(zalpha,alpha,par->Fp);
 unsigned char zch[32]; 
 fmpz2chars(zch,zalpha);
 int ms=crypto_scalarmult_ristretto255_base(ga,zch);

//calculate the queries
 
 fq_t one, varq;
 fq_init(one,par->Fp);
 fq_one(one,par->Fp);
 fq_init(varq,par->Fp);
 
 //choose t random vectors r1,r2,...,rt and t random field elements gamma_1,...,gamma_t
 fq_mat_t *r;
 r=malloc(sizeof(fq_mat_t)*par->t);
 
 fq_t *gamma;
 gamma=malloc(sizeof(fq_t)*par->t);
 
 for (int i=0;i<par->t;i++)
 {
    fq_mat_init(r[i],par->m,1,par->Fp);
    fq_mat_randtest(r[i],state,par->Fp);
    fq_init(gamma[i],par->Fp);
    fq_randtest(gamma[i],state,par->Fp);
 }

 //calculate the curve CU
 fq_poly_t x;
 fq_poly_init(x,par->Fp);
 fq_poly_gen(x,par->Fp);

 fq_poly_t varp1, varp2;
 fq_poly_init(varp1,par->Fp);
 fq_poly_init(varp2,par->Fp);

 fq_poly_t *CU;
 CU=malloc(sizeof(fq_poly_t)*par->m);
 
 for (int i=0;i<par->m;i++)
 {
   fq_poly_init(CU[i],par->Fp);
   fq_poly_set_fq(CU[i],fq_mat_entry(X,i,0),par->Fp);
   for (int j=0;j<par->t;j++)
   {
      fq_poly_pow(varp1,x,j+1,par->Fp);
      fq_poly_set_fq(varp2,fq_mat_entry(r[j],i,0),par->Fp);
      fq_poly_mul(varp1,varp1,varp2,par->Fp);
      fq_poly_add(CU[i],CU[i], varp1,par->Fp);
   }
 }
 
 //calculate the polynomial BU
 fq_poly_t BU;
 fq_poly_init(BU,par->Fp);
 fq_poly_set_fq(BU,alpha,par->Fp);
 
 for (int j=0;j<par->t;j++)
 {
   fq_poly_pow(varp1,x,j+1,par->Fp);
   fq_poly_set_fq(varp2,gamma[j],par->Fp);
   fq_poly_mul(varp1,varp1,varp2,par->Fp);
   fq_poly_add(BU,BU,varp1,par->Fp);
 }  
 
 //calculate the queries
 for (int i=0;i<par->k;i++)
 {
    fq_mat_init(c[i],par->m,1,par->Fp);
    fq_init(b[i],par->Fp);
    
    for (int j=0;j<par->m;j++)
    {
       fq_poly_evaluate_fq(varq,CU[j],par->SID[i],par->Fp);
       fq_mat_entry_set(c[i],j,0,varq,par->Fp);
    }
    fq_poly_evaluate_fq(varq,BU,par->SID[i],par->Fp);
    fq_set(b[i],varq,par->Fp);
 }
}




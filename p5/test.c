#include "p5.h"
 
int main(int argc, char **argv)
{
////The public parameters 
pubpar *par;
par=malloc(sizeof(pubpar));
//p is the number of elements in the coefficient ring of F
fmpz_set_str(par->p,"7237005577332262213973186563042994240857116359379907606001950938285454250989",10);
fq_ctx_init(par->Fp,par->p,1,"t");
//The number of monomials in F
fmpz_t N;
fmpz_init(N);
//The state for generating random numbers
flint_rand_t state;
flint_randinit(state);
//The field element 1 in F_q
fq_t one;
//The function F
fq_mat_t F;
//The input X
fq_mat_t X;
//The query curve c(u)
fq_mat_t *c;
//The ver curve b(u)
fq_t *b;
//The function value
fq_t* y;
//The ver value
fq_t* z;
//The ver result
int ver;
//ver key
unsigned char *ga;
//The time counters
clock_t s0,t0,s1,t1,s2,t2,s3,t3,s4,t4;
double T_kg, T_pg, T_ct, T_vf, T_us,T_im, T_mx, mx;
//The native result
fq_t zy;
fq_init(zy,par->Fp);
//The time counters
clock_t sn,tn;
double T_n;


//parameters for choosing F and X
fmpz_t zf, bound;
fmpz_init(zf);
fmpz_init(bound);
fmpz_set_str(bound,"7237005577332262213973186563042994240857116359379907606001950938285454250989",10); 
//F_q variables
fq_t zq;
fq_init(zq,par->Fp);

//set d and t
par->d=2;
int RP=10;

for (int zt=3;zt<4;zt++)
for (int zm=1;zm<2;zm++)
{
par->t=zt;
par->k=(par->d+1)*par->t+1;
par->m=2000;

//N=C(m+d,d) is the number of monomials in F
fmpz_bin_uiui(N,par->m+par->d,par->d);
if (fmpz_cmp_ui(N,50000000)>0)
   continue;
par->iN=fmpz_get_ui(N); 

par->iN=200000;


//SID is the set of IDs of servers
fq_init(one,par->Fp);
fq_one(one,par->Fp);
par->SID=malloc(sizeof(fq_t)*par->k);
for (int i=0;i<par->k;i++)
{   
   fq_init(par->SID[i],par->Fp);
   fq_mul_ui(par->SID[i],one,i+1,par->Fp);   
}

////choose F
fq_mat_init(F,par->iN,1,par->Fp);
fmpz_set_str(zf,"6939831890859379238386960392965905689913722297935665337785291190938693596201",10);
//can be simplified if the bound is p
for (int i=0;i<par->iN;i++)
{
    fmpz_add_ui(zf,zf,1);
    fq_set_fmpz(zq,zf,par->Fp);
    fq_mat_entry_set(F,i,0,zq,par->Fp);
}


//choose X
fq_mat_init(X,par->m,1,par->Fp);
fmpz_set_str(zf,"5630933909596291259927897519658151851995601860579757376333565665993969197191",10);
//can be simplified if the bound is p
for (int i=0;i<par->m;i++)
{
    fmpz_add_ui(zf,zf,1);
    fq_set_fmpz(zq,zf,par->Fp);
    fq_mat_entry_set(X,i,0,zq,par->Fp);
}

T_n=(double)0;
sn=clock();
for (int rp=0;rp<RP;rp++)
{
  //  Eval(zy,F,X,par);
}
tn=clock();
// output the running time of all algorithms    
T_n=T_n+(double) (tn-sn)/CLOCKS_PER_SEC;



T_kg=(double)0;
T_pg=(double)0;
T_ct=(double)0;
T_vf=(double)0;
T_us=(double)0;
T_mx=(double)0;

for (int rep=0;rep<RP;rep++)
{
////Test the KeyGen algorithm//////////////////////////////////////////
//l is the random line
s1=clock();
  KeyGen(F);
t1=clock();

////Test the ProbGen algorithm/////////////////////////////////////////
c=malloc(sizeof(fq_mat_t)*par->k);
for (int i=0;i<par->k;i++)
    fq_mat_init(c[i],par->m,1,par->Fp);

b=malloc(sizeof(fq_t)*par->k);
for (int i=0;i<par->k;i++)
    fq_init(b[i],par->Fp);
//malloc ga
ga=malloc(sizeof(unsigned char)*32);
//execute the ProbGen algorithm
s2=clock();
  ProbGen(ga, c, b, X, par);
t2=clock();

////Test the Compute algorithm
y=malloc(sizeof(fq_t)*par->k);
z=malloc(sizeof(fq_t)*par->k);
for (int i=0;i<par->k;i++)
{
  fq_init(y[i],par->Fp);
  fq_init(z[i],par->Fp);
  fq_rand_not_zero(y[i],state,par->Fp);
  fq_rand_not_zero(z[i],state,par->Fp);
}

mx=0;
for (int i=0;i<par->k;i++)
{
  s0=clock();
  // Compute(y[i], z[i], F, c[i], b[i], par);
  t0=clock();
  if ((t0-s0)> mx) 
     mx=t0-s0;  
  T_ct=T_ct+(double) (t0-s0)/CLOCKS_PER_SEC;
}

fq_mat_t L0;
			fq_mat_init(L0,par->k,1,par->Fp);
			// fq_mat_entry(L0,j,0): L_j(0)= \prod_{k=1,k!=j}^m k/(k-j)
			fq_t ksubj;
	fq_init(ksubj,par->Fp);

	for (int j = 0; j < par->k; j ++){
		fq_one(fq_mat_entry(L0,j,0),par->Fp);

		for (int i = 1; i <= par->k; i++){
			if(i != (j + 1)){
				fq_mul_ui(fq_mat_entry(L0,j,0),fq_mat_entry(L0,j,0),i,par->Fp);
				fq_set_si(ksubj,i-(j+1),par->Fp);
				fq_div(fq_mat_entry(L0,j,0),fq_mat_entry(L0,j,0),ksubj,par->Fp);
			}
		}
	}

////Test the Verify algorithm
s4=clock();
  ver=Verify(ga, y,z,par,L0);
t4=clock();

// output the running time of all algorithms       
T_kg=T_kg+(double) (t1-s1)/CLOCKS_PER_SEC;
T_pg=T_pg+(double) (t2-s2)/CLOCKS_PER_SEC;
T_vf=T_vf+(double) (t4-s4)/CLOCKS_PER_SEC;
T_mx=T_mx+(double) mx/CLOCKS_PER_SEC;
}
T_us=T_pg+T_vf;
printf("m=%4d, d=%2d, N=%10d, t=%d, k=%3d, T_kg=%12.6f, T_pg=%12.6f, T_vf=%12.6f, T_ct=%12.6f, T_n=%12.6f, T_us=%12.6f\n", par->m, 
par->d, par->iN, par->t, par->k, T_kg/RP, T_pg/RP, T_vf/RP, T_ct/RP, T_n/RP, T_us/RP);


fq_mat_clear(F,par->Fp);
fq_mat_clear(X,par->Fp);

}


flint_randclear(state); 
 
 
printf("\n");
return 0;
}

#include <time.h>
#include "fmpz.h"
#include "fq.h"
#include "fmpz_poly.h"
#include "fq_poly.h"
#include "stdio.h"
#include "gmp.h"
#include "stdlib.h"
#include "string.h"
#include "fq_vec.h"
#include "fq_mat.h"
#include <sodium.h>

////The public parameters used by the protocol
typedef struct
{
  int m;
  int d;
  int t;
  int k;
  int iN;
  int ***I;
  fmpz_t p;
  fq_ctx_t Fp;
  fq_t * SID;  
} pubpar;


////compute y=F(X)
void Eval(fq_t y, fq_mat_t F, fq_mat_t X, pubpar *par);
////Interpolate f such that f(SID)=Y
void IntPoly(fq_poly_t f, fq_t *Y, pubpar *par);
////compute op=c*op
void fq_mat_scal_mul(fq_mat_t op, fq_t c, fq_ctx_t Fp);
////convert from fq to fmpz
void fq2fmpz(fmpz_t out, fq_t in, fq_ctx_t ctx);
////convert unsigned char [] to fmpz
void chars2fmpz(fmpz_t zo, unsigned char r[]);
////convert fmpz to unsigned char []
void fmpz2chars(unsigned char r[], fmpz_t zo);



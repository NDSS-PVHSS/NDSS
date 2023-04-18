#include "p5.h"

////The Compute algorithm
void Compute(fq_t yi, fq_t zi, fq_mat_t F, fq_mat_t ci, fq_t bi, pubpar* par)
{
   Eval(yi,F,ci,par);
   fq_mul(zi,yi,bi,par->Fp);
}


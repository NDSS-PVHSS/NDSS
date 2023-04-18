#include "mpoly.h"
 

////convert F to pk=l=l(u), vk=f(u)=F(l(u))
void KeyGen(fq_mat_t F);
////use pk=l(u) to convert X to (vk,sigma)
void ProbGen(unsigned char *ga, fq_mat_t *c, fq_t *b, fq_mat_t X, pubpar *par);
////compute pii=F(sigmai)
void Compute(fq_t yi, fq_t zi, fq_mat_t F, fq_mat_t ci, fq_t bi, pubpar* par);
////use (f,vk) to verify the servers' responses pi
int Verify(unsigned char *ga, fq_t *y, fq_t *z, pubpar *par,fq_mat_t L0);
 

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <unistd.h>
#include <malloc.h>
#include <gmp.h>
#include <string.h>
#include <flint/fq.h>
#include <flint/fq_mat.h>
#include <flint/fq_vec.h>
#include <flint/fq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_mpoly.h>

int n,t,d,m;
ulong N;

fq_ctx_t Fp;

flint_rand_t state;

bool DEBUG;

typedef struct{
	fq_mat_t s;
	fq_mat_t * a;
} input_struct;

// typedef struct{

// } output_struct;

// ******* Initialize public parameters ********
void init(int nn, int dd){

	// DEBUG = true;
	DEBUG = false;
	n = nn;   // length of data
	t = 2;    // privacy threshold
	d = dd;   // degree of polynomial f
	m = ceil(((d+1)*t+1) / 2.0);    // number of servers
	if (m <= 2 * t + 1) {
		m = 2 * t + 1;
	}
}

void L_Share(fq_mat_t x, fq_mat_t * s);
void L_Eval(int j,fq_mat_t f,fq_mat_t s_j,fq_t out_j);
int L_Ver(fq_mat_t L0, fq_t * out, fq_t y);
void L_interpolation(fmpz_t * locs, fmpz_t * vals, fmpz_t res);

void next(int _length, int a[]);
void Eval(fq_mat_t f, fq_mat_t x, fq_t y);
void Eval_diff(fq_mat_t f, fq_mat_t x, int _index,fq_t y);
void La_intpoly_coeff(fq_mat_t L0);

// ************ Share algorithm ********************
void H_Share(fq_mat_t x, fq_mat_t * s, fq_mat_t ** a){
	// x[n];
	// s[m][n];
	// a[t][m][n]

	fq_t temp;
	fq_init(temp, Fp);
	
	/*
	  Generating phi(u)

	  phi(u) is a degree-t polynomial with phi(0) = x
	  Here, we split phi(u)= [phi_0(u),...,phi_{n-1}(u)]
	  phi_i(0) = x_i
	*/
	fq_poly_t phi[n];
	for (int i = 0; i < n; i ++){
		fq_poly_init(phi[i],Fp);

		fq_poly_set_coeff(phi[i],0,fq_mat_entry(x,i,0),Fp);
		for (int k = 1; k <= t; k ++){
			fq_set_ui(temp, 0, Fp);
			while (fq_is_zero(temp, Fp)){
				fq_randtest(temp, state, Fp);
			}
			fq_poly_set_coeff(phi[i], k, temp, Fp);
		}	
	}

	/*
	  Computing s = [s_0,...s_{m-1}]	 
	  s_j = phi(j) 
	*/
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < m; j ++){
            fq_set_ui(temp,j+1,Fp);
            fq_poly_evaluate_fq(fq_mat_entry(s[j],i,0),phi[i],temp,Fp);
        }
    }  

    // Set r_i and use Pi_L to share them (L_share)
    // r[t][n]

    fq_mat_t r[t];  
    // r_1,...,r_{t} ;  r[k] = r_{k+1}
    
    for (int k = 0; k < t; k ++){
    	fq_mat_init(r[k], n, 1, Fp);
    	for (int i = 0; i < n; i ++){
    		fq_poly_get_coeff(fq_mat_entry(r[k], i, 0), phi[i], k + 1, Fp);
    	}
    	L_Share(r[k], a[k]);
    }

	if (DEBUG) {
		printf("**Debug in H_share**\n");
		printf("x = ");
		fq_mat_print_pretty(x, Fp);
		printf("\n");
		for (int i = 0; i < n; i ++){
			printf("phi_%d=",i);
			fq_poly_print_pretty(phi[i],  "X", Fp);
			printf("\n");
		}
	}
}

// // ************ Eval algorithm *********************
void Evalpro(fq_mat_t f, fq_mat_t x, fq_t y, fq_mat_t sigma){

	int index=0; //index of f

    // constant term
	if(d>=0){
		fq_set(y,fq_mat_entry(f,0,0),Fp);
		fq_mat_zero(sigma,Fp);
	}

	fq_t y1;
   	fq_init(y1,Fp);

	fq_t y2;
   	fq_init(y2,Fp);

	fq_t y3;
   	fq_init(y3,Fp);

	// degree-1 terms
	if(d>=1){
		for(int i1=0;i1<n;i1++){
			index=index+1;
			fq_set(y1,fq_mat_entry(f,index,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i1,0),Fp);
			fq_add(y,y,y1,Fp);
			
			fq_add(fq_mat_entry(sigma,i1,0),fq_mat_entry(sigma,i1,0),fq_mat_entry(f,index,0),Fp);
		}
	}

	// degree-2 terms
	if(d>=2){
		for(int i1=0;i1<n;i1++){ 
			for(int i2=i1;i2<n;i2++){
				index=index+1;
				if(index<N){
					fq_set(y1,fq_mat_entry(f,index,0),Fp);
					fq_mul(y2,y1,fq_mat_entry(x,i1,0),Fp); 
					fq_mul(y3,y2,fq_mat_entry(x,i2,0),Fp);
					fq_add(y,y,y3,Fp);

					if(i1==i2){
						fq_mul_ui(y3,y2,2,Fp);
						fq_add(fq_mat_entry(sigma,i1,0),fq_mat_entry(sigma,i1,0),y3,Fp);
					} 
					else{
						fq_mul(y3,y1,fq_mat_entry(x,i2,0),Fp);
						fq_add(fq_mat_entry(sigma,i1,0),fq_mat_entry(sigma,i1,0),y3,Fp);
						fq_add(fq_mat_entry(sigma,i2,0),fq_mat_entry(sigma,i2,0),y2,Fp);
					}
				}

			}
		}
	}
}

void H_Eval(int j,fq_mat_t f,fq_mat_t s_j,fq_t y_j, fq_mat_t sigma_j){
	// f: input length-n
	//    output length-1
	//    degree-d
	// sigma_j[n]


    // Eval f(s_j)
	Eval(f,s_j,y_j);
	// Eval df/dX_i (s_j)
	// for (int index = 0; index < n; index ++){
	// 	Eval_diff(f, s_j, index ,fq_mat_entry(sigma_j, index, 0));
	// }

	Evalpro(f,s_j,y_j,sigma_j);


	if (DEBUG){
		printf("**Debuging in H_Eval** \n");
		printf("f = ");
		fq_mat_print_pretty(f, Fp);
		printf("\n");
		
		printf("s_j = ");
		fq_mat_print_pretty(s_j, Fp);
		printf("\n");

		printf("y_j = ");
		fq_print_pretty(y_j,Fp);
		printf("\n");

		printf("sigma_j = ");
		fq_mat_print_pretty(sigma_j, Fp);
		printf("\n");
	}
}




// // ************ Ver algorithm **********************


			// L_j(u) = \prod_{k=1,k!=j}^m (k-u)/(k-j)
			// fq_mat_entry(H0,j,0): H_j(0)= (1-2(u-j) L'_j(j))L_j^2(0)
			// fq_mat_entry(_H0,j,0): \bar H_j(0) = ( -j) L_j^2(0)



void He_intpoly_coeff(fq_mat_t L0, fq_mat_t H0, fq_mat_t _H0){
	// compute H_j(0) and \bar H_j(0)

	fq_t temp;
	fq_t jsubk;
	fq_t Lsquare;
	fq_t Lj;
	fq_t onefq;
	fq_init(onefq,Fp);
	fq_one(onefq,Fp);

	fq_init(temp, Fp);
	fq_init(jsubk, Fp);
	fq_init(Lsquare, Fp);
	fq_init(Lj, Fp);

	for (int j = 0; j < m; j ++){
		fq_mul(Lsquare, fq_mat_entry(L0, j, 0), fq_mat_entry(L0, j, 0), Fp); 
		
		fq_zero(Lj, Fp);  // fq_set_ui(Lj, 2 * (j + 1), Fp);
		for (int k = 0; k < m; k ++){
			if (k != j){
				fq_set_si(jsubk, j - k, Fp);
				fq_inv(temp, jsubk, Fp);
				fq_add(Lj, Lj, temp, Fp);// fq_mul(Lj, Lj, temp, Fp);
			}
		}

		// printf("L'%d(%d)=",j+1,j+1);
		// fq_print_pretty(Lj,Fp);
		// printf("\n");

		fq_mul_ui(Lj,Lj,2 * (j + 1),Fp);
		fq_add(Lj, Lj, onefq, Fp);

		fq_mul(fq_mat_entry(H0, j, 0), Lj, Lsquare, Fp); // fq_mul(fq_mat_entry(H0, j, 0), Lj, fq_mat_entry(L0, j, 0), Fp);
		// fq_add(fq_mat_entry(H0, j, 0), fq_mat_entry(H0, j, 0), Lsquare, Fp);
		// H_j(0)=(1+2j L'_j(j))(L_j(0))^2
		fq_set_si(temp, -(j + 1) , Fp);
		fq_mul(fq_mat_entry(_H0, j, 0), temp, Lsquare, Fp);
		// \bar H_j(0) = (-j)(L_j(0))^2
	}

	// printf("first H'j(0)=");
	// 	fq_mat_print_pretty(H0,Fp);
	// 	printf("\n");

	fq_clear(onefq,Fp);
	fq_clear(temp, Fp);
	fq_clear(jsubk, Fp);
	fq_clear(Lsquare, Fp);
	fq_clear(Lj, Fp);
}



int H_Ver(fq_mat_t L0, fq_mat_t H0, fq_mat_t _H0, fq_t * out, fq_mat_t * out2, fq_mat_t ** out3, fq_t y){
	// out[m]
	// out2[m][n]
	// out3[t][m][n]
	// return 1 if y is correct; otherwise return 0
	
	fq_t zero;
	fq_init(zero,Fp);
	fq_zero(zero,Fp);


	fq_t temp;
	fq_t temp2;
	fq_init(temp,Fp);
	fq_init(temp2,Fp);

	fq_t ** r = malloc(sizeof(fq_t *) * t);
	for (int k = 0; k < t; k ++){
		r[k] = malloc(sizeof(fq_t) * n);
		for (int i = 0; i < n; i ++){
			fq_init(r[k][i], Fp);
		}
	}
	fq_t *  temp_vec = malloc(sizeof(fq_t) * m);
	int IsCorrect;
	for (int k = 0; k < t; k ++){
		for (int i = 0; i < n; i ++){
			for (int j = 0; j < m; j ++){
				fq_init(temp_vec[j], Fp);
				fq_set(temp_vec[j], fq_mat_entry(out3[k][j], i, 0), Fp);
			}
			IsCorrect = L_Ver(L0, temp_vec, r[k][i]);
			// printf("R Is correct: %d\n",IsCorrect);
		}
	}
	


	fq_t ** dphi = malloc(sizeof(fq_t *) * n);
	for (int i = 0; i < n; i ++){
		dphi[i] = malloc(sizeof(fq_t) * m);
		for (int j = 0; j < m; j ++){
			fq_init(dphi[i][j], Fp);
			fq_zero(dphi[i][j], Fp);
			for (int k = 0; k < t; k ++){     // k = 1..t in the paper
				fq_set_ui(temp, j+1, Fp);       // j
				fq_pow_ui(temp, temp, k, Fp); // j^{k-1}

				fq_mul_ui(temp, temp, k + 1, Fp); // kj^{k-1}
				fq_mul(temp, temp, r[k][i], Fp);

				fq_add(dphi[i][j], dphi[i][j], temp, Fp);

			}
			// printf("dphi_%d(%d)=",i+1,j+1);
			// 	fq_print_pretty(dphi[i][j],Fp);
			// 	printf("\n");
		}
	}

    // Generating psi(u)
	// psi: input: length-1
	//     output: length-1
	//     degree-t
	//     psi(0) = 0
	// dpsi: derivation of psi
	fq_poly_t psi;
	fq_poly_t dpsi;
	fq_poly_init(psi, Fp);
	fq_poly_init(dpsi, Fp);

	fq_poly_randtest(psi, state, t + 1, Fp);
	fq_poly_set_coeff(psi, 0, zero, Fp);
	fq_poly_derivative(dpsi, psi, Fp);

	fq_t * dif_vals_y = malloc(sizeof(fq_t) * m);
	fq_t * dif_vals_z = malloc(sizeof(fq_t) * m);


	for (int j = 0; j < m; j ++){
		fq_init(dif_vals_y[j], Fp);
		fq_zero(dif_vals_y[j], Fp);
		for (int i = 0; i < n; i ++){
			fq_mul(temp, dphi[i][j], fq_mat_entry(out2[j], i, 0), Fp);
			fq_add(dif_vals_y[j], dif_vals_y[j], temp, Fp);
		}		

		// printf("F'(%d)=",j+1);
		// fq_print_pretty(dif_vals_y[j],Fp);
		// printf("\n");
	}

	for (int j = 0; j < m; j ++){
		fq_init(dif_vals_z[j], Fp);
		// fq_zero(dif_vals_z[j], Fp);
        fq_set_ui(temp, j + 1,Fp);							// j

		fq_poly_evaluate_fq(temp2,psi,temp,Fp);			// psi(j)
		fq_mul(dif_vals_z[j], temp2, dif_vals_y[j], Fp);			// psi(j)F'(j)
		// for (int i = 0; i < n; i ++){
        //     fq_poly_evaluate_fq(temp2,psi,temp,Fp);			// psi(j)
        //     fq_mul(temp2, temp2, dphi[i][j], Fp);			// phi'_i(j)psi(j)
        //     fq_mul(temp2, temp2, fq_mat_entry(out2[j], i, 0), Fp);// phi'_i(j)psi(j)f'_i(s_j)
        //     fq_add(dif_vals_z[j], dif_vals_z[j], temp2, Fp);
		// }		
	    fq_poly_evaluate_fq(temp2,dpsi,temp,Fp);			// psi'(j)
        fq_mul(temp2, temp2, out[j], Fp);					// psi'(j)f(s_j)
        fq_add(dif_vals_z[j], dif_vals_z[j], temp2, Fp);
	}

	//compute psiy: psi(j)y_j
	fq_t psiy[m];
	for(int j = 0; j < m; j ++){
		fq_init(psiy[j],Fp);
		fq_set_ui(temp,j + 1,Fp);
		fq_poly_evaluate_fq(psiy[j],psi,temp, Fp);
		fq_mul(psiy[j],psiy[j],out[j], Fp);
	}

	// printf("Lj(0):");
	// fq_mat_print_pretty(L0,Fp);
	// printf("\n");

	// printf("Hj(0):");
	// fq_mat_print_pretty(H0,Fp);
	// printf("\n");

	// printf("\bar Hj(0):");
	// fq_mat_print_pretty(_H0,Fp);
	// printf("\n");


	//interpolate y=F(0)
	fq_zero(y,Fp);
	for(int j = 0; j < m; j++){
		fq_mul(temp, out[j],fq_mat_entry(H0,j,0),Fp);
		fq_add(y,y,temp,Fp);

		fq_mul(temp, dif_vals_y[j] ,fq_mat_entry(_H0,j,0),Fp);
		fq_add(y,y,temp,Fp);
	}


	//interpolate z=G(0)
	fq_t z;
	fq_init(z,Fp);
	fq_zero(z,Fp);	
	for(int j = 0; j < m; j++){
		fq_mul(temp, psiy[j],fq_mat_entry(H0,j,0),Fp);
		fq_add(z,z,temp,Fp);

		fq_mul(temp, dif_vals_z[j] ,fq_mat_entry(_H0,j,0),Fp);
		fq_add(z,z,temp,Fp);
	}
		
	int res = fq_is_zero(z,Fp);



	fq_poly_clear(dpsi, Fp);
	fq_poly_clear(psi, Fp);
	fq_clear(temp, Fp);
	fq_clear(temp2, Fp);
	fq_clear(zero,Fp);
	fq_clear(z, Fp);
	for (int k = 0; k < t; k ++){
		for (int i = 0; i < n; i ++){
			fq_clear(r[k][i], Fp);
		}
	}
	for(int j = 0; j < m; j ++){
		fq_clear(psiy[j],Fp);
		fq_clear(dif_vals_z[j], Fp);
		fq_clear(dif_vals_y[j], Fp);
	}
	for (int i = 0; i < n; i ++){
		for (int j = 0; j < m; j ++){
			fq_clear(dphi[i][j], Fp);
		}
	}
	for (int j = 0; j < m; j ++){
		fq_clear(temp_vec[j], Fp);
	}			
	return res;
}


int main(){
    flint_randinit(state);

    fmpz_t p;
	fmpz_init(p);
    fmpz_set_str(p,"340282366920938463463374607431768211297",10);  // 128 bit prime
    // fmpz_set_str(p,"11",10);
	fq_ctx_init(Fp,p,1,"gen");
	fmpz_clear(p);

	int n_list[9];
	n_list[0] = 1413;
	n_list[1] = 180;
	n_list[2] = 68;
	n_list[3] = 39;
	n_list[4] = 27;
	n_list[5] = 21;
	n_list[6] = 17;
	n_list[7] = 15;
	n_list[8] = 13;

	int new_n_list[7];
	int new_d_list[7];
	new_n_list[0] = 11; new_d_list[0] = 12;
	new_n_list[1] = 10; new_d_list[1] = 13;
	new_n_list[2] = 9; new_d_list[2] = 15;
	new_n_list[3] = 8; new_d_list[3] = 17;
	new_n_list[4] = 7; new_d_list[4] = 21;
	new_n_list[5] = 6; new_d_list[5] = 27;
	new_n_list[6] = 5; new_d_list[6] = 39;

	ulong NN = 200000;
	FILE *infp = fopen("input_size.out","w");
	FILE *outfp = fopen("output_size.out","w");


	for (int index = 1; index <= 10; index ++){
		// int nn = n_list[index];
		// int dd = index + 2;
		int nn = 2000;
    	int dd = 2;
    	N = NN * index;
    	// N = nn * (nn + 1) / 2;


		printf("n = %d, d = %d, N = %ld\n", nn, dd, N);

	    for (int times = 1; times <= 5; times ++){
	    	printf("times = %d\n", times);

			init(nn,dd);  // Initialize public parameters

			fq_t temp;
			fq_init(temp, Fp);
			
			// Data x
		    fq_mat_t x;
			fq_mat_init(x,n,1,Fp);
		    // for (int i = 0; i < n ; i ++){
		    // 	fq_set_ui(temp, 0, Fp);
		    // 	while (fq_is_zero(temp, Fp)){
		    // 		fq_randtest(temp, state, Fp);
		    // 	}
		    // 	fq_set(fq_mat_entry(x, i, 0), temp, Fp);
		    // }
			
		    // fq_mat_randtest(x,state,Fp);

			fmpz_t mv;
			fmpz_init_set(mv,p);
    		
			for (int i=0;i<n;i++)
			{
    			fmpz_add_ui(mv,mv,1);
    			fq_set_fmpz(fq_mat_entry(x,i,0),mv,Fp);
			}
		     

		    if (DEBUG) {
			    printf("x: ");
			    fq_mat_print_pretty(x,Fp);
			    printf("\n");
		    }

		    // Shares s_1,...,s_j (in_j)
			fq_mat_t * in;
		    in=malloc(sizeof(fq_mat_t) * m);
		    for(int j = 0;j < m; j ++)
				fq_mat_init(in[j],n,1,Fp);
			fq_mat_t ** in2;
			in2 = malloc(sizeof(fq_mat_t *) * t);
			for (int k = 0; k < t; k ++){
				in2[k] = malloc(sizeof(fq_mat_t) * m);
				for (int j = 0; j < m; j ++){
					fq_mat_init(in2[k][j], n, 1, Fp);
				}
			}


		    // n-variate function f of degree d
			// N=C(n+d,d) is the number of monomials in f
			// fmpz_t Nm;
			// fmpz_init(Nm);
			// fmpz_bin_uiui(Nm,n+d,d);
			// N=fmpz_get_ui(Nm);

		    fq_mat_t f;
			fq_mat_init(f,N,1,Fp);
			// fq_mat_randtest(f,state,Fp); //random choose
			
			for (int i=0;i<N;i++)
			{
    			fmpz_sub_ui(mv,mv,1);
    			fq_set_fmpz(fq_mat_entry(f,i,0),mv,Fp);
			}
			fmpz_clear(mv);

			// fq_mat_print_pretty(f,Fp);




		    // Output of servers (out_j)
		    fq_t * out=malloc(sizeof(fq_t)*m);
		    for(int j=0;j<m;++j){
		        fq_init(out[j],Fp);
				fq_zero(out[j],Fp);
		    }

		    fq_mat_t * out2 = malloc(sizeof(fq_mat_t) * m);
		    for (int j = 0; j < m; j ++){
		    	fq_mat_init(out2[j], n, 1, Fp);
		    }

		//**********************************
			
		    // Share
			clock_t share_start,share_end;
			double share_time = 0;

			share_start = clock();
			H_Share(x,in, in2);
			share_end = clock();
			share_time = (double) (share_end- share_start)/CLOCKS_PER_SEC;
		    printf("Time of Share: %f ms\n",share_time*1000);


			// Eval
			// printf("***Evaling****\n");	

			clock_t eval_start,eval_end;
			double eval_time;
		    eval_time = 0;

		    double total_eval_time = 0;
			for (int j = 0;j < m; j ++){
			    eval_start=clock();
				H_Eval(j, f, in[j], out[j], out2[j]);
			    eval_end=clock();
			    total_eval_time += (double) (eval_end-eval_start)/CLOCKS_PER_SEC;
			    if ( (double) (eval_end-eval_start)/CLOCKS_PER_SEC > eval_time ) {
			    	eval_time = (double) (eval_end-eval_start)/CLOCKS_PER_SEC;
			    }
		    }
		    printf("Time of Max_Eval: %f ms\n",eval_time*1000);
		    printf("Time of total Eval: %f ms\n",total_eval_time*1000);


			
			// // Verify
			fq_mat_t L0, H0, _H0;
			fq_mat_init(L0,m,1,Fp);
			fq_mat_init(H0,m,1,Fp);
			fq_mat_init(_H0,m,1,Fp);
			// L_j(u) = \prod_{k=1,k!=j}^m (k-u)/(k-j)
			// fq_mat_entry(H0,j,0): H_j(0)= (1-2(u-j) L'_j(j))L_j^2(0)
			// fq_mat_entry(_H0,j,0): \bar H_j(0) = ( -j) L_j^2(0)
						
			La_intpoly_coeff(L0);
			He_intpoly_coeff(L0,H0,_H0); 

			// if (DEBUG) {
			// 	printf("L_1(0)..L_m(0):");
			// 	fq_mat_print_pretty(L0,Fp);
			// }  
			
			

		    int IsCorrect;
		    fq_t y;
		    fq_init(y,Fp);
		    clock_t Ver_start,Ver_end;
		    Ver_start=clock();
		    IsCorrect = H_Ver(L0, H0,_H0, out, out2, in2 ,y);
		    Ver_end=clock();
			
			if (DEBUG) {
				printf("Ac? %d, y=",IsCorrect);
				fq_print_pretty(y,Fp);
				printf("\n");
			}
			printf("Ac? %d, y=",IsCorrect);
				fq_print_pretty(y,Fp);
				printf("\n");

		    double Ver_time= (double) (Ver_end-Ver_start)/CLOCKS_PER_SEC;
		    printf("Time of Ver: %f ms\n",Ver_time*1000);


			fq_t fx;
		    fq_init(fx,Fp);
		    // Eval f(x)
		    clock_t fx_start,fx_end;
		    fx_start=clock();
		    Eval(f,x,fx);
		    fx_end=clock();

			if (DEBUG) {
				printf("f(x)=");
				fq_print_pretty(fx,Fp);
				printf("\n");
			}
			printf("f(x)=");
				fq_print_pretty(fx,Fp);
				printf("\n");

		    double fx_time= (double) (fx_end-fx_start)/CLOCKS_PER_SEC;
		    printf("Time of directly eval f(x): %f ms\n",fx_time*1000);

		    // communication cost
		    for (int j = 0; j < m; j ++){
		    	fq_mat_fprint(infp, in[j], Fp);
		    	for (int k = 0; k < t; k ++){
		    		fq_mat_fprint(infp, in2[k][j], Fp);
		    	}
		    	fq_fprint(outfp, out[j], Fp);
				fq_mat_fprint(outfp, out2[j], Fp);
				for (int k = 0; k < t; k ++){
					fq_mat_fprint(outfp, in2[k][j], Fp);
				}		    	
		    }


			// clear memory
		    fq_mat_clear(f,Fp);
			fq_mat_clear(x,Fp);
			for(int j=0;j<m;++j){
		        fq_clear(out[j],Fp);
		    }
			// fq_mat_clear(L0,Fp);
			// fq_clear(y,Fp);
		    fq_clear(fx,Fp);
		    
		    for(int j = 0;j < m; j ++){
				fq_mat_clear(in[j],Fp);
				for (int k = 0; k < t ; k ++){
					fq_mat_clear(in2[k][j], Fp);
				}
		    }


	    }
		fprintf(infp,"\n\n\n");
		fprintf(outfp,"\n\n\n");
	}

	fclose(infp);
	fclose(outfp);

    return 0;
}


































// ************ Share algorithm ********************
void L_Share(fq_mat_t x, fq_mat_t * s){
	// x[1..n];
	// s[1..m][1..n];

	fq_t temp;
	fq_init(temp,Fp);

	/*
	  Generating phi(u)

	  phi(u) is a degree-t polynomial with phi(0) = x
	  Here, we split phi(u)= [phi_0(u),...,phi_{n-1}(u)]
	  phi_i(0) = x_i
	*/
	fq_poly_t phi[n];
	for (int i = 0; i < n; i ++){
		fq_poly_init(phi[i],Fp);

		fq_poly_set_coeff(phi[i],0,fq_mat_entry(x,i,0),Fp);
		for (int k = 1; k <= t; k ++){
			fq_set_ui(temp, 0, Fp);
			while (fq_is_zero(temp, Fp)){
				fq_randtest(temp, state, Fp);
			}
			fq_poly_set_coeff(phi[i], k, temp, Fp);
		}	
	}

	/*
	  Computing s = [s_1,...s_m]	 
	  s_j = phi(j) 
	*/
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < m; j ++){
            fq_set_ui(temp,j+1,Fp);
            fq_poly_evaluate_fq(fq_mat_entry(s[j],i,0),phi[i],temp,Fp);
        }
    }  


	if (DEBUG) {
		printf("**Debuging in L_share**:\n");
		printf("x = ");
		fq_mat_print_pretty(x, Fp);
		printf("\n");
		for (int i = 0; i < n; i ++){
			printf("phi_%d=",i);
			fq_poly_print_pretty(phi[i],"X",Fp);
			printf("\n");
		}
	}


	for (int i = 0; i < n; i ++){
		fq_poly_clear(phi[i], Fp);
	}
	fq_clear(temp,Fp);
}

// ************ Eval algorithm *********************
void next(int _length, int a[]){
	a[_length - 1] += 1;
	for (int i = _length - 1; i >= 1; i --){
		if (a[i] >= n) {
			a[i - 1] += 1;
			for (int j = i; j <= _length - 1; j ++){
				a[j] = a[i - 1];
			}
		}
	}
}

void Eval(fq_mat_t f, fq_mat_t x, fq_t y){

	int index=0; //index of f

    // constant term
	if(d>=0){
		fq_set(y,fq_mat_entry(f,0,0),Fp);
	}

	fq_t y1;
   	fq_init(y1,Fp);

	// degree-1 terms
	if(d>=1){
		for(int i1=0;i1<n;i1++){
			index=index+1;
			fq_set(y1,fq_mat_entry(f,index,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i1,0),Fp);
			fq_add(y,y,y1,Fp);
		}
	}

	// degree-2 terms
	if(d>=2){
		for(int i1=0;i1<n;i1++)
		for(int i2=i1;i2<n;i2++)
		{
			index=index+1;
			if(index<N){
			fq_set(y1,fq_mat_entry(f,index,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i1,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i2,0),Fp);
			fq_add(y,y,y1,Fp);
			}
		}
	}

	// degree-3 terms
	if(d>=3){
		for(int i1=0;i1<n;i1++)
		for(int i2=i1;i2<n;i2++)
		for(int i3=i2;i3<n;i3++)
		{
			index=index+1;
			if(index<N){
			fq_set(y1,fq_mat_entry(f,index,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i1,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i2,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i3,0),Fp);
			fq_add(y,y,y1,Fp);
			}
		}
	}

	// degree-4 terms
	if(d>=4){
		for(int i1=0;i1<n;i1++)
		for(int i2=i1;i2<n;i2++)
		for(int i3=i2;i3<n;i3++)
		for(int i4=i3;i4<n;i4++)
		{
			index=index+1;
			if(index<N){
			fq_set(y1,fq_mat_entry(f,index,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i1,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i2,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i3,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i4,0),Fp);
			fq_add(y,y,y1,Fp);
			}
		}
	}

	// degree-5 terms
	if(d>=5){
		for(int i1=0;i1<n;i1++)
		for(int i2=i1;i2<n;i2++)
		for(int i3=i2;i3<n;i3++)
		for(int i4=i3;i4<n;i4++)
		for(int i5=i4;i5<n;i5++)
		{
			index=index+1;
			if(index<N){
			fq_set(y1,fq_mat_entry(f,index,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i1,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i2,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i3,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i4,0),Fp);
			fq_mul(y1,y1,fq_mat_entry(x,i5,0),Fp);
			fq_add(y,y,y1,Fp);
			}
		}
	}

}

// void Eval(fq_mat_t f, fq_mat_t x, fq_t y){

// 	int index=0; //index of f
// 	fq_t y1;
//    	fq_init(y1,Fp);

// 	int a[d]; 
// 	for (int i = 0; i < d; i ++){
// 		a[i] = 0;
// 	}

// 	while ( (a[0] < n) & (index < N) ){
// 		fq_set(y1,fq_mat_entry(f,index,0),Fp);
// 		for (int j = 0; j < d; j ++){
// 			fq_mul(y1,y1,fq_mat_entry(x,a[j],0),Fp);
// 		}
// 		fq_add(y,y,y1,Fp);
// 		next(d,a);
// 		index += 1;
// 	}

// }

void Eval_diff(fq_mat_t f, fq_mat_t x, int _index,fq_t y){
	// set y = df(X)/d(X[_index]) | X=x

	int index=0; //index of f
	fq_t y1;
   	fq_init(y1,Fp);

	int a[d];
	for (int i = 0; i < d; i ++){
		a[i] = 0;
	}
	while ( ( a[0] < n ) & ( index < N ) ){
		// if (index % 100000 == 0){printf("N = %d, index = %d\n",N, index);}
		fq_set(y1,fq_mat_entry(f,index,0),Fp);
		ulong count = 0;
		for (int j = 0; j < d; j ++){
			if (a[j] == _index){
				count ++;
				if (count != 1){
					fq_mul(y1,y1,fq_mat_entry(x,a[j],0),Fp);				
				}
			}
			else{
				fq_mul(y1,y1,fq_mat_entry(x,a[j],0),Fp);
			}
		}
		fq_mul_ui(y1, y1, count, Fp);
		fq_add(y,y,y1,Fp);
		next(d,a);
		index += 1;
	}
}
void L_Eval(int j,fq_mat_t f,fq_mat_t s_j,fq_t out_j){
	// f: input length-n
	//    output length-1
	//    degree-d

    // Eval f(s_j)
	Eval(f,s_j,out_j);
}

// ************ Ver algorithm **********************

void La_intpoly_coeff(fq_mat_t L0){
	//interpolate f such that f(j)=y(j), return res=f(0)
	// locs[1..m]
	// vals[1..m]

	// compute L_j(0) (:=fq_mat_entry(L0,j,0))
	fq_t ksubj;
	fq_init(ksubj,Fp);

	for (int j = 0; j < m; j ++){
		fq_one(fq_mat_entry(L0,j,0),Fp);

		for (int k = 1; k <= m; k++){
			if(k != (j + 1)){
				fq_mul_ui(fq_mat_entry(L0,j,0),fq_mat_entry(L0,j,0),k,Fp);
				fq_set_si(ksubj,k-(j+1),Fp);
				fq_div(fq_mat_entry(L0,j,0),fq_mat_entry(L0,j,0),ksubj,Fp);
			}
		}
	}
}

int L_Ver(fq_mat_t L0, fq_t * out, fq_t y){
	// out[1..m]
	// return 1 if y is correct; otherwise return 0
	
    // Generating psi(u)
	// psi: input: length-1
	//     output: length-1
	//     degree-t
	//     psi(0) = 0
	fq_poly_t psi;
	fq_poly_init(psi,Fp);
	fq_poly_randtest(psi, state, t + 1, Fp);

	fq_t zero;
	fq_init(zero,Fp);
	fq_zero(zero,Fp);
	fq_poly_set_coeff(psi,0,zero,Fp);
	fq_clear(zero,Fp);


	fq_t temp;
	fq_init(temp,Fp);

	//interpolate y=F(0)
	fq_zero(y,Fp);
	for(int j = 0; j < m; j++){
		fq_mul(temp,out[j],fq_mat_entry(L0,j,0),Fp);
		fq_add(y,y,temp,Fp);
	}

	//interpolate z=G(0)
	//compute psiout: psi(j)f(s_j)
	fq_t psiout[m];
	for(int j = 0; j < m; j ++){
		fq_init(psiout[j],Fp);
		fq_set_ui(temp,j+1,Fp);
		fq_poly_evaluate_fq(psiout[j],psi,temp,Fp);
		fq_mul(psiout[j],psiout[j],out[j],Fp);
	}

	fq_t z;
	fq_init(z,Fp);
	fq_zero(z,Fp);	
	for(int j = 0; j < m; j++){
		fq_mul(temp,psiout[j],fq_mat_entry(L0,j,0),Fp);
		fq_add(z,z,temp,Fp);
	}
	

	int res = fq_is_zero(z,Fp);


	fq_poly_clear(psi, Fp);
	fq_clear(temp, Fp);
	fq_clear(z, Fp);
	for (int j = 0; j < m; j ++){
		fq_clear(psiout[j], Fp);
	}

	return res;
}






#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <unistd.h>
#include <malloc.h>
#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_mpoly.h>

int n,t,d,m;

fmpz_mod_mpoly_t f;

fmpz_mod_ctx_t univarctx;
fmpz_mod_mpoly_ctx_t nvarsctx;

fmpz_t p;

flint_rand_t state;

bool DEBUG;

// ******* Initialize public parameters ********
void init(int nn, int dd){

	// DEBUG = true;
	DEBUG = false;
	n = nn;           // length of data
	t = 3;            // privacy threshold
	d = dd;            // degree of polynomial f
	m = (d+1)*t+1;    // number of servers
    flint_randinit(state);

	fmpz_init(p);
    fmpz_set_str(p,"340282366920938463463374607431768211297",10);  // 128 bit prime
    // fmpz_set_str(p,"17",10);
    fmpz_mod_ctx_init(univarctx,p);
    fmpz_mod_mpoly_ctx_init(nvarsctx,n,ORD_LEX,p);
}



// ************ Share algorithm ********************
void L_Share(fmpz_t * x, fmpz_t ** s){
	// x[1..n];
	// s[1..m][1..n];

	fmpz_t temp;
	fmpz_init(temp);

	
	/*
	  Generating phi(u)

	  phi(u) is a degree-t polynomial with phi(0) = x
	  Here, we split phi(u)= [phi_1(u),...,phi_n(u)]
	  phi_i(0) = x_i
	*/
	fmpz_mod_poly_t phi[n + 1];
	for (int i = 1; i <= n; i ++){
		fmpz_mod_poly_init(phi[i],univarctx);
		fmpz_mod_poly_randtest(phi[i], state, t + 1, univarctx);
		fmpz_mod_poly_set_coeff_fmpz(phi[i],0,x[i],univarctx);	
	}

	/*
	  Computing s = [s_1,...s_m]	 
	  s_j = phi(j) 
	*/
    for(int i = 1; i <= n; i ++){
        for(int j = 1; j <= m; j ++){
            fmpz_set_ui(temp,j);
            fmpz_mod_poly_evaluate_fmpz(s[j][i],phi[i],temp,univarctx);
        }
    }  


	if (DEBUG) {
		for (int i = 1; i <= n; i ++){
			printf("phi_%d=",i);
			fmpz_mod_poly_print_pretty(phi[i],"X",univarctx);
			printf("\n");
		}
	}


	for (int i = 1; i <= n; i ++){
		fmpz_mod_poly_clear(phi[i], univarctx);
	}
	fmpz_clear(temp);
}

// ************ Eval algorithm *********************
void L_Eval(int j,fmpz_mod_mpoly_t f,fmpz_t * s_j,fmpz_t out_j){
	// f: input length-n
	//    output length-1
	//    degree-d

	// transform the tpye of s_j into "fmpz *const *" to evaluate

    fmpz ** vals;

    vals=FLINT_ARRAY_ALLOC(n,fmpz *);
    for(int i = 0;i < n; i ++){  // Evaluation require index start at 0
        vals[i]=FLINT_ARRAY_ALLOC(1,fmpz);
        fmpz_init_set(vals[i],s_j[i + 1]);
    }

    // Eval f(s_j)

    fmpz_mod_mpoly_evaluate_all_fmpz(out_j,f,vals,nvarsctx);

    if (DEBUG) {
    	printf("out_%d = ",j );
    	fmpz_print(out_j);
    	printf("\n");
    }

    for(int i = 0;i < n; i ++){  
        fmpz_clear(vals[i]);
    }
}

// ************ Ver algorithm **********************

void L_interpolation(fmpz_t * locs, fmpz_t * vals, fmpz_t res){
	// locs[1..m]
	// vals[1..m]

	fmpz_mod_poly_t L[m + 1];
	for (int i = 1; i <= m; i ++){
		fmpz_mod_poly_init(L[i],univarctx);
	}
	fmpz_t temp, _sum, zero;
	fmpz_init(temp);
	fmpz_init_set_ui(_sum, 0);
	fmpz_init_set_ui(zero, 0);

	// Compute L[i] (L_i(0))
	for (int i = 1; i <= m; i ++){
		// L[i] = L_i(0) = \prod_{k \neq i} (-u_k)/(u_i-u_k)
		fmpz_t roots[m - 1];
		for (int k = 1; k <= m; k ++){
			if (k < i) {
				fmpz_init_set(roots[k - 1],locs[k]);  // vector roots should start at index 0
			}
			else if (k > i) {
				fmpz_init_set(roots[k - 2], locs[k]); // vector roots should start at index 0
			}
		}
		if (DEBUG){
			printf("i = %d\n", i);
			printf("roots = ");
			for (int k = 0; k < m - 2; k ++){
				fmpz_print(roots[k]);
				printf(" ");
			}
			printf("\n");
		}
	
		fmpz_mod_poly_product_roots_fmpz_vec(L[i], roots, m - 1, univarctx);     
		
		if (DEBUG){
			printf("L_%d(X) = ", i);
			fmpz_mod_poly_print_pretty(L[i],"X", univarctx);
			printf("\n");

			for (int j = 1; j <= m; j ++){
				fmpz_mod_poly_evaluate_fmpz(temp, L[i], locs[j], univarctx);

				printf("L_%d(",i);
				fmpz_print(locs[j]);
				printf(")=");
				fmpz_print(temp);
				printf("\n");
			}
		}

		fmpz_mod_poly_evaluate_fmpz(temp, L[i], locs[i], univarctx);
		fmpz_mod_poly_scalar_div_fmpz(L[i], L[i], temp, univarctx); // normalize
	}

	// res = F(0) = \sum_{i} vals[i]*L_i[0]
	for (int i = 1; i <= m; i ++){
		fmpz_mod_poly_evaluate_fmpz(temp, L[i], zero, univarctx);
		fmpz_mod_addmul(_sum, _sum, temp, vals[i], univarctx);
	}
	fmpz_set(res, _sum);


	for (int i = 1; i <= m; i ++){
		fmpz_mod_poly_clear(L[i], univarctx);
	}
	fmpz_clear(temp);
	fmpz_clear(_sum);
	fmpz_clear(zero);
}
int L_Ver(fmpz_t * out, fmpz_t y){
	// out[1..m]
	// return 1 if y is correct; otherwise return 0
	
    // Generating psi(u)
	// psi: input: length-1
	//     output: length-1
	//     degree-t
	//     psi(0) = 0
	fmpz_mod_poly_t psi;
	fmpz_mod_poly_init(psi,univarctx);
	fmpz_mod_poly_randtest(psi, state, t + 1, univarctx);
	fmpz_mod_poly_set_coeff_ui(psi,0,0,univarctx);

	fmpz_t temp;
	fmpz_init(temp);


	fmpz_t locs[m + 1];
	for (int i = 1; i <= m; i ++){
		fmpz_init_set_ui(locs[i],i);
	}
	
	fmpz_t psiout[m + 1];
	for (int i = 1; i <= m; i ++){
		fmpz_mod_poly_evaluate_fmpz(temp,psi,locs[i],univarctx); //temp = psi(loc[i])
		fmpz_mod_mul(psiout[i],temp,out[i],univarctx);          //psiout[i] = psi(loc[i])*out[i]
	}

	fmpz_t z;
	fmpz_init(z);

	if (DEBUG) {
		printf("*Interpolating*\n");
	}
	L_interpolation(locs, out, y);
	L_interpolation(locs, psiout, z);

	int res = fmpz_equal_ui(z, 0);


	fmpz_mod_poly_clear(psi, univarctx);
	fmpz_clear(temp);
	fmpz_clear(z);
	for (int i = 1; i <= m; i ++){
		fmpz_clear(locs[i]);
	}

	return res;
}


int main(){
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

    for (int index = 0; index <= 8; index ++){
		int nn = n_list[index];
		int dd = index + 2;
		printf("n = %d, d = %d\n", nn,dd);

	    for (int times = 1; times <= 5; times ++){
			init(nn,dd);  // Initialize public parameters

			// Data x
		    fmpz_t * x=malloc(sizeof(fmpz_t) * (n+1));
		    for(int i = 1;i <= n; i ++){
		        fmpz_init(x[i]);
		        fmpz_randm(x[i],state,p);
		    }    
		    if (DEBUG) {
			    printf("x: ");
			    for (int i = 1; i <= n; i ++){
			        fmpz_print(x[i]);
			        printf(" ");
			    }
			    printf("\n");
		    }
		    
		    // Shares s_1,...,s_j (in_j)
		    fmpz_t ** in=malloc(sizeof(fmpz_t *)*(m+1));
		    for(int j = 1;j <= m; j ++){
		        in[j]=malloc(sizeof(fmpz_t)*(n+1));
		        for(int i = 1;i <= n; i ++){
		            fmpz_init(in[j][i]);
		        }
		    }


		    // Function f
		    // Here we use f = (X_1 + ..+ X_n + 1)^d
		    fmpz_mod_mpoly_t f;
			fmpz_mod_mpoly_init(f,nvarsctx);

		    ulong * exp = malloc(sizeof(ulong *) * n);
		    for(int i = 0; i < n; i ++){
		        for (int j = 0; j < n; j ++){
		        	exp[j] = 0;
		        }
		        exp[i] = 1;    
		        fmpz_mod_mpoly_set_coeff_ui_ui(f,1,exp,nvarsctx);
		    }
		    for (int j = 0; j < n; j ++){
		    	exp[j] = 0;
		    }
		    fmpz_mod_mpoly_set_coeff_ui_ui(f,1,exp,nvarsctx);
		    fmpz_mod_mpoly_pow_ui(f,f,d,nvarsctx);


		    if (DEBUG) {
		    	char ** str = (char **)malloc(sizeof(char *) * n);

		    	for (int i = 0; i < n; i ++){
		    		str[i] = malloc(sizeof(char) * 10);
		    		sprintf(str[i], "X_%d", i);
		    	}
		    	printf("f=");
		    	fmpz_mod_mpoly_print_pretty(f, (const char **)str, nvarsctx);
		    	printf("\n");
		    }


		    // Output of servers (out_j)
		    fmpz_t * out=malloc(sizeof(fmpz_t)*(m+1));
		    for(int j=1;j<m+1;++j){
		        fmpz_init_set_ui(out[j],0);
		    }



		//**********************************
			
		    // Share

		    // printf("***Sharing***\n");
			clock_t share_start,share_end;
			double share_time = 0;

			share_start = clock();
			L_Share(x,in);
			share_end = clock();
			share_time = (double) (share_end- share_start)/CLOCKS_PER_SEC;
		    printf("Time of Share: %f ms\n",share_time*1000);


			// Eval
			// printf("***Evaling****\n");	

			clock_t eval_start,eval_end;
			double eval_time;
		    eval_time = 0;

		    double total_eval_time = 0;
			for (int j = 1;j <= m; j ++){
			    eval_start=clock();
		        L_Eval(j,f,in[j],out[j]);
			    eval_end=clock();
			    total_eval_time += (double) (eval_end-eval_start)/CLOCKS_PER_SEC;
			    if ( (double) (eval_end-eval_start)/CLOCKS_PER_SEC > eval_time ) {
			    	eval_time = (double) (eval_end-eval_start)/CLOCKS_PER_SEC;
			    }
		    }
		    printf("Time of Max_Eval: %f ms\n",eval_time*1000);
		    printf("Time of total Eval: %f ms\n",total_eval_time*1000);

		    // Ver
		    // printf("***Verifying***\n");
		    int IsCorrect;
		    fmpz_t y;
		    fmpz_init(y);
		    clock_t dec_start,dec_end;
		    dec_start=clock();
		    IsCorrect = L_Ver(out,y);
		    dec_end=clock();
		    fmpz_clear(y);
		    double dec_time= (double) (dec_end-dec_start)/CLOCKS_PER_SEC;
		    printf("Time of Dec: %f ms\n",dec_time*1000);


			for(int j=1;j<m+1;++j){
		        fmpz_clear(out[j]);
		    }



			fmpz_t fx;
		    fmpz_init(fx);
		    fmpz ** vals;
		    vals = FLINT_ARRAY_ALLOC(n,fmpz *);
		    for(int i = 0;i < n; i ++){  // Evaluation require index start at 0
		        vals[i]=FLINT_ARRAY_ALLOC(1,fmpz);
		        fmpz_init_set(vals[i],x[i + 1]);
		    }
		    // Eval f(x)
		    clock_t fx_start,fx_end;
		    fx_start=clock();
		    fmpz_mod_mpoly_evaluate_all_fmpz(fx,f,vals,nvarsctx);
		    fx_end=clock();


		    double fx_time= (double) (fx_end-fx_start)/CLOCKS_PER_SEC;
		    printf("Time of directly eval f(x): %f ms\n",fx_time*1000);

		    fmpz_mod_mpoly_clear(f,nvarsctx);
		    fmpz_clear(fx);
		    for (int i = 0; i < n; i ++){
		    	fmpz_clear(vals[i]);
		    }
		    for(int j = 1;j <= m; j ++){
		        for(int i = 1;i <= n; i ++){
		            fmpz_clear(in[j][i]);
		        }
		    }
		    for(int i = 1;i <= n; i ++){
		        fmpz_clear(x[i]);
		    }
		}
	}
    return 0;
}
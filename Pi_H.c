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
void init(int nn){

	DEBUG = true;
	// DEBUG = false;
	n = nn;    // length of data
	t = 3;    // privacy threshold
	d = 2;    // degree of polynomial f
	m = ((d+1)*t+2) / 2;    // number of servers
	if (m <= 2 * t + 1) {
		m = 2 * t + 1;
	}
    flint_randinit(state);

	fmpz_init(p);
    // fmpz_set_str(p,"340282366920938463463374607431768211297",10);  // 128 bit prime
    fmpz_set_str(p,"97",10);
    fmpz_mod_ctx_init(univarctx,p);
    fmpz_mod_mpoly_ctx_init(nvarsctx,n,ORD_LEX,p);
}

void L_Share(fmpz_t * x, fmpz_t ** s);
void L_Eval(int j,fmpz_mod_mpoly_t f,fmpz_t * s_j,fmpz_t out_j);
int L_Ver(fmpz_t * out, fmpz_t y);
void L_interpolation(fmpz_t * locs, fmpz_t * vals, fmpz_t res);

// ************ Share algorithm ********************
void H_Share(fmpz_t * x, fmpz_t ** s, fmpz_t *** ss){
	// x[1..n];
	// s[1..m][1..n];
	// ss[1..t][1..m][1..n];

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

    // Set r_i and use Pi_L to share them (L_share)
    // ss[k][j] = s_{k,j}

    fmpz_t r[t + 1][n + 1];
    for (int k = 1; k <= t ; k ++){
    	for (int i = 1; i <= n; i ++){
	    	fmpz_init(r[k][i]);
    		fmpz_mod_poly_get_coeff_fmpz(r[k][i], phi[i], k, univarctx);
    	}
	    L_Share(r[k], ss[k]);
    } 


	if (DEBUG) {
		for (int i = 1; i <= n; i ++){
			printf("phi_%d=",i);
			fmpz_mod_poly_print_pretty(phi[i],"X",univarctx);
			printf("\n");
		}
	}

}

// ************ Eval algorithm *********************
void H_Eval(int j,fmpz_mod_mpoly_t f,fmpz_t * s_j,fmpz_t ** ss_j, fmpz_t out1_j, fmpz_t * out2_j, fmpz_t ** out3_j){
	// f: input length-n
	//    output length-1
	//    degree-d

	// s_j[1..n]
	// ss_j[1..t][1..n]
	// out1_j
	// out2_j[1..n]
	// out3_j[1..t][1..n]


	// transform the tpye of s_j into "fmpz *const *" to evaluate

    fmpz ** vals;
    vals=FLINT_ARRAY_ALLOC(n,fmpz *);
    for(int i = 0;i < n; i ++){  // Evaluation require index start at 0
        vals[i]=FLINT_ARRAY_ALLOC(1,fmpz);
        fmpz_init_set(vals[i],s_j[i + 1]);
    }


    // Eval f(s_j)
    fmpz_mod_mpoly_evaluate_all_fmpz(out1_j,f,vals,nvarsctx);


	// Eval df/dX_i (s_j)
	fmpz_mod_mpoly_t df;
    fmpz_mod_mpoly_init(df,nvarsctx);	
	for (int i = 0; i < n; i ++){
	    fmpz_mod_mpoly_derivative(df,f,i,nvarsctx);
        fmpz_mod_mpoly_evaluate_all_fmpz(out2_j[i],df,vals,nvarsctx);	    
	}



    // out3_j = ss_j
    for (int k = 1; k <= t; k ++){
    	for (int i = 1; i <= n; i ++){
    		fmpz_set(out3_j[k][i], ss_j[k][i]);
    	}
    }


    // if (DEBUG) {
    // 	printf("out_%d = ",j );
    // 	fmpz_print(out_j);
    // 	printf("\n");
    // }
}

// ************ Ver algorithm **********************


void H_interpolation(fmpz_t * locs, fmpz_t * vals, fmpz_t * dif_vals, fmpz_t res){
	// locs[1..m]
	// vals[1..m]
	// dif_vals[1..m]

	// Here, we set L[j] = L_j(u) a polynomial
	//              H[j] = H_j(0) a value
	//        H_tilde[j] = \tilde{H}_j(0) a value
	fmpz_mod_poly_t * L = malloc(sizeof(fmpz_mod_poly_t) * (m + 1));
	fmpz_t * H = malloc(sizeof(fmpz_t) * (m + 1));
	fmpz_t * H_tilde = malloc(sizeof(fmpz_t) * (m + 1));
	
	for (int j = 1; j <= m; j ++){
		fmpz_mod_poly_init(L[j],univarctx);
		fmpz_init(H[j]);
		fmpz_init(H_tilde[j]);
	}
	fmpz_t temp, _sum, zero;
	fmpz_init(temp);
	fmpz_init_set_ui(_sum, 0);
	fmpz_init_set_ui(zero, 0);

	// Compute polynomials L[j]  ( L_j(u) ), 
	//                     H[j]  ( H_j(u) ),
	//               H_tilde[j]  ( \tilde{H}_j(u) )



	for (int j = 1; j <= m; j ++){
		// L[j] = L_j(0) = \prod_{k \neq j} (-u_k)/(u_j-u_k)
		fmpz_t roots[m - 1];
		for (int k = 1; k <= m; k ++){
			if (k < j) {
				fmpz_init_set(roots[k - 1],locs[k]);  // vector roots should start at index 0
			}
			else if (k > j) {
				fmpz_init_set(roots[k - 2], locs[k]); // vector roots should start at index 0
			}
		}
		if (DEBUG){
			printf("j = %d\n", j);
			printf("roots = ");
			for (int k = 0; k <= m - 2; k ++){
				fmpz_print(roots[k]);
				printf(" ");
			}
			printf("\n");
		}
	
		fmpz_mod_poly_product_roots_fmpz_vec(L[j], roots, m - 1, univarctx);     
		
		if (DEBUG){
			printf("L_%d(X) = ", j);
			fmpz_mod_poly_print_pretty(L[j],"X", univarctx);
			printf("\n");

			for (int k = 1; k <= m; k ++){
				fmpz_mod_poly_evaluate_fmpz(temp, L[j], locs[k], univarctx);

				printf("L_%d(",j);
				fmpz_print(locs[k]);
				printf(")=");
				fmpz_print(temp);
				printf("\n");
			}
		}
		fmpz_mod_poly_evaluate_fmpz(temp, L[j], locs[j], univarctx);
		fmpz_mod_poly_scalar_div_fmpz(L[j], L[j], temp, univarctx); // normalize

	}


	// H_j(u) = (1 - 2(u - u_j) L'_j(u_j) ) ( L_j(u) )^2
	// H_j(0) = (1 - 2(  - u_j) L'_j(u_j) ) ( L_j(0) )^2

	// \tilde{H}_j(u) = (u - u_j) (L_j(u))^2
	// \tilde{H}_j(0) = (  - u_j) (L_j(0))^2
	for (int j = 1; j <= m; j ++){
		// H_j(u) = (1 - 2(u - u_j) L'_j(u_j) ) ( L_j(u) )^2
		// H_j(0) = (1 - 2(  - u_j) L'_j(u_j) ) ( L_j(0) )^2

		// \tilde{H}_j(u) = (u - u_j) (L_j(u))^2
		// \tilde{H}_j(0) = (  - u_j) (L_j(0))^2

		fmpz_t temp_L0;                // L_j(0)
		fmpz_mod_poly_t temp_Ldif;     // L'_j(u)
		fmpz_t temp_Ldif_uj;           // L'_j(u_j)

		fmpz_init(temp_L0);
		fmpz_mod_poly_init(temp_Ldif, univarctx);
		fmpz_init(temp_Ldif_uj);
		
		fmpz_mod_poly_evaluate_fmpz(temp_L0, L[j], zero, univarctx);
		fmpz_mod_poly_derivative(temp_Ldif, L[j], univarctx);
		fmpz_mod_poly_evaluate_fmpz(temp_Ldif_uj, temp_Ldif, locs[j], univarctx);


		// H_j(0) = (1 - 2(  - u_j) L'_j(u_j) ) ( L_j(0) )^2

		fmpz_set_ui(H[j], 1);
		fmpz_mod_mul(temp, locs[j], temp_Ldif_uj, univarctx); // u_j * L'_j(u_j)
		fmpz_mod_add_fmpz(H[j], H[j], temp, univarctx);
		fmpz_mod_add_fmpz(H[j], H[j], temp, univarctx);		
		fmpz_mod_mul(H[j], H[j], temp_L0, univarctx);
		fmpz_mod_mul(H[j], H[j], temp_L0, univarctx);

		// \tilde{H}_j(0) = (  - u_j) (L_j(0))^2		
		fmpz_mod_sub(H_tilde[j], zero, locs[j], univarctx);
		fmpz_mod_mul(H_tilde[j], H_tilde[j], temp_L0, univarctx);
		fmpz_mod_mul(H_tilde[j], H_tilde[j], temp_L0, univarctx);		
	}

	// res = F(0) = \sum_{j} vals[j]*H_j[0] + dif_vals[j]* \tilde{H}_j[0]
	for (int j = 1; j <= m; j ++){
		fmpz_mod_addmul(_sum, _sum, vals[j], H[j], univarctx);
		fmpz_mod_addmul(_sum, _sum, dif_vals[j], H_tilde[j], univarctx);
	}
	fmpz_set(res, _sum);

}
int H_Ver(fmpz_t * out1, fmpz_t ** out2, fmpz_t *** out3, fmpz_t y){
	// out1[1..m]
	// out2[1..m][1..n]
	// out3[1..m][1..t][1..n]

	// return 1 if y is correct; otherwise return 0


	fmpz_t temp;
	fmpz_init(temp);

	fmpz_t ** r = malloc(sizeof(fmpz_t *) * (t + 1));
	for (int k = 1; k <= t; k ++){
		r[k] = malloc(sizeof(fmpz_t) * (n + 1));
		for (int i = 1; i <= n; i ++){
			fmpz_init(r[k][i]);
		}
	}

	fmpz_t * temp_vec = malloc(sizeof(fmpz_t) * (m + 1));
	int IsCorrect;
	for (int k = 1; k <= t; k ++){
		for (int i = 1; i <= n; i ++){
			for (int j = 1; j <= m; j ++){
				fmpz_init_set(temp_vec[j], out3[j][k][i]);
			}
			IsCorrect = L_Ver(temp_vec, r[k][i]);
			if (IsCorrect == 0) {return 0;}
		}
	}


	// phi_ij[i][j] = phi'_i(j)
	fmpz_t ** phi_ij = malloc(sizeof(fmpz_t *) * (n + 1));

	for (int i = 1; i <= n; i ++){
		phi_ij[i] = malloc(sizeof(fmpz_t) * (m + 1));
		for (int j = 1; j <= m; j ++){
			fmpz_init_set_ui(phi_ij[i][j], 0);
			
			for (int k = 1; k <= t; k ++){
				fmpz_set_ui(temp,j);                            // j
				fmpz_mod_pow_ui(temp, temp, k - 1, univarctx);  // j^{k-1}
				fmpz_mul_ui(temp, temp, k);                     // kj^{k-1}
				fmpz_mod_mul(temp, temp, r[k][i],univarctx);    // r_{k,i}kj^{k-1}

				fmpz_mod_add(phi_ij[i][j], phi_ij[i][j], temp, univarctx);
			}
		}
	}


    // Generating psi(u)
	// psi: input: length-1
	//     output: length-1
	//     degree-t
	//     psi(0) = 0
	// dpsi: derivation of psi
	fmpz_mod_poly_t psi, dpsi;
	fmpz_mod_poly_init(psi,univarctx);
	fmpz_mod_poly_init(dpsi,univarctx);
	fmpz_mod_poly_randtest(psi, state, t + 1, univarctx);
	fmpz_mod_poly_set_coeff_ui(psi,0,0,univarctx);
	fmpz_mod_poly_derivative(dpsi, psi, univarctx);


	fmpz_t * locs = malloc(sizeof(fmpz_t) * (m + 1));
	for (int i = 1; i <= m; i ++){
		fmpz_init_set_ui(locs[i],i);
	}
	

	fmpz_t * psiout = malloc(sizeof(fmpz_t) * (m + 1));
	for (int i = 1; i <= m; i ++){
		fmpz_mod_poly_evaluate_fmpz(temp,psi,locs[i],univarctx); //temp = psi(loc[i])
		fmpz_mod_mul(psiout[i],temp,out1[i],univarctx);          //psiout[i] = psi(loc[i])*out[i]
	}

	fmpz_t * dif_vals_y = malloc(sizeof(fmpz_t) * (m + 1));
	fmpz_t * dif_vals_z = malloc(sizeof(fmpz_t) * (m + 1));


	for (int j = 1; j <= m; j ++){
		fmpz_init_set_ui(dif_vals_y[j], 0);
		for (int i = 1; i <= n; i ++){
			fmpz_mod_addmul(dif_vals_y[j], dif_vals_y[j], phi_ij[i][j], out2[j][i], univarctx);
		}		
	}

	for (int j = 1; j <= m; j ++){
		fmpz_init_set_ui(dif_vals_z[j], 0);
		for (int i = 1; i <= n; i ++){
			fmpz_set_ui(temp, j);                                    // j
			fmpz_mod_poly_evaluate_fmpz(temp, psi, temp, univarctx); // psi(j)
			fmpz_mod_mul(temp, temp, phi_ij[i][j], univarctx);       // phi'_i(j)psi(j)
			fmpz_mod_mul(temp, temp, out2[j][i], univarctx);         // phi'_i(j)psi(j)f'_i(s_j)

			fmpz_mod_add(dif_vals_z[j], dif_vals_z[j], temp, univarctx);
		}		
		fmpz_set_ui(temp, j);                                     // j
		fmpz_mod_poly_evaluate_fmpz(temp, dpsi, temp, univarctx); // psi'(j)
		fmpz_mod_mul(temp, temp, out1[j],univarctx);              // psi'(j)f(s_j)
		fmpz_mod_add(dif_vals_z[j], dif_vals_z[j], temp, univarctx);
	}



	if (DEBUG) {
		printf("*Interpolating*\n");
	}

	fmpz_t z;
	fmpz_init(z);
	H_interpolation(locs, out1, dif_vals_y,  y);
	H_interpolation(locs, psiout, dif_vals_z, z);


	return fmpz_equal_ui(z, 0);
}


int main(){
	for (int nn = 3; nn <= 200; nn = nn + 200){
	for (int nnn = 1; nnn <= 1; nnn ++){
    	printf("***************************************************************\n n = %d \n",nn);
    	printf("times = %d\n",nnn);

	init(nn);  // Initialize public parameters
	
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

	fmpz_t *** in2 = malloc(sizeof(fmpz_t **) * (t + 1));
    for(int k = 1; k <= t; k ++){
    	in2[k] = malloc(sizeof(fmpz_t *) * (m + 1));
	    for(int j = 1;j <= m; j ++){
	    	in2[k][j] = malloc(sizeof(fmpz_t) * (n + 1));
	        for(int i = 1;i <= n; i ++){
	            fmpz_init(in2[k][j][i]);
	        }
	    }
	}

	fmpz_t *** in3 = malloc(sizeof(fmpz_t **) * (m + 1));
	for (int j = 1; j <= m; j ++){
		in3[j] = malloc(sizeof(fmpz_t *) * (t + 1));
		for (int k = 1; k <= t; k ++){
			in3[j][k] = malloc(sizeof(fmpz_t) * (n + 1));
			for (int i = 1; i <= n; i ++){
				fmpz_init(in3[j][k][i]);
			}
		}
	}



    // Function f
    // Here we use f = (1 * X_1 + ..+ n * X_n + 1)^d
    fmpz_mod_mpoly_t f;
	fmpz_mod_mpoly_init(f,nvarsctx);

    ulong * exp = malloc(sizeof(ulong *) * n);
    for(int i = 0; i < n; i ++){
        for (int j = 0; j < n; j ++){
        	exp[j] = 0;
        }
        exp[i] = 1;    
        fmpz_mod_mpoly_set_coeff_ui_ui(f,i + 1,exp,nvarsctx);
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
    fmpz_t * out1 = malloc(sizeof(fmpz_t) * (m + 1));
    fmpz_t ** out2 = malloc(sizeof(fmpz_t *) * (m + 1));
    fmpz_t *** out3 = malloc(sizeof(fmpz_t **) * (m + 1));
    
    for (int j = 1; j <= m; j ++){
        fmpz_init_set_ui(out1[j],0);

        out2[j] = malloc(sizeof(fmpz_t) * (n + 1));
        for (int i = 1; i <= n; i ++ ){
        	fmpz_init_set_ui(out2[j][i], 0);
        }

        out3[j] = malloc(sizeof(fmpz_t *) * (t + 1));
        for (int k = 1; k <= t; k ++){
        	out3[j][k] = malloc(sizeof(fmpz_t) * (n + 1));
        	for (int i = 1; i <= n; i ++){
				fmpz_init_set_ui(out3[j][k][i],0);        		
        	}
        }
    }



//**********************************

    // Share
    printf("***Sharing***\n");
	clock_t share_start,share_end;
	double share_time = 0;

	share_start = clock();
	H_Share(x,in,in2);
	share_end = clock();
	share_time = (double) (share_end- share_start)/CLOCKS_PER_SEC;
    printf("Time of Share: %f ms\n",share_time*1000);



    for(int k = 1; k <= t; k ++){
	    for(int j = 1;j <= m; j ++){
	        for(int i = 1;i <= n; i ++){
	 			fmpz_set(in3[j][k][i], in2[k][j][i]);
	        }
	    }
	}

	// Eval
	printf("***Evaling***\n");	

	clock_t eval_start,eval_end;
	double eval_time;
    eval_time = 0;

    double total_eval_time = 0;
	for (int j = 1;j <= m; j ++){
	    eval_start=clock();
        H_Eval(j,f,in[j],in3[j], out1[j], out2[j], out3[j]);
	    eval_end=clock();
	    total_eval_time += (double) (eval_end-eval_start)/CLOCKS_PER_SEC;
	    if ( (double) (eval_end-eval_start)/CLOCKS_PER_SEC > eval_time ) {
	    	eval_time = (double) (eval_end-eval_start)/CLOCKS_PER_SEC;
	    }
    }
    printf("Time of Eval: %f ms\n",eval_time*1000);
    printf("Time of total Eval: %f ms\n",total_eval_time*1000);


    // Ver
    printf("***Verifying***\n");
    int IsCorrect;
    fmpz_t y;
    fmpz_init(y);
    clock_t dec_start,dec_end;
    dec_start=clock();
    IsCorrect = H_Ver(out1,out2,out3,y);
    dec_end=clock();
    double dec_time= (double) (dec_end-dec_start)/CLOCKS_PER_SEC;
    printf("Time of Dec: %f ms\n",dec_time*1000);


    if (DEBUG)
    {
    printf("****************\n");
    printf("x = ");
    for (int i = 1; i <= n; i ++){
        fmpz_print(x[i]);
        printf(" ");
    }
    printf("\n");


	char ** str = (char **)malloc(sizeof(char *) * n);
	for (int i = 0; i < n; i ++){
		str[i] = malloc(sizeof(char) * 10);
		sprintf(str[i], "X_%d", i);
	}
	printf("f = ");
	fmpz_mod_mpoly_print_pretty(f, (const char **)str, nvarsctx);
	printf("\n");

	fmpz_t res;
	fmpz_init(res);
    fmpz ** valx;
    valx=FLINT_ARRAY_ALLOC(n,fmpz *);
    for(int i=0;i<n;++i){
        valx[i]=FLINT_ARRAY_ALLOC(1,fmpz);
        fmpz_init_set(valx[i],x[i+1]);
    }	
	fmpz_mod_mpoly_evaluate_all_fmpz(res, f, valx, nvarsctx);
    printf("f(x) = ");
    fmpz_print(res);
    printf("\n");
    printf("y = ");
    fmpz_print(y);
    printf("\n");
    printf("IsCorrect = %d\n",IsCorrect);
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

	}
	}
    return 0;
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
	if (DEBUG) {
		printf("*Interpolating*\n");
	}
	L_interpolation(locs, psiout, z);


	return fmpz_equal_ui(z, 0);
}

























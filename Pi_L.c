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
void init(){

	// DEBUG = true;
	DEBUG = false;
	n = 2;    // length of data
	t = 4;    // privacy threshold
	d = 2;    // degree of polynomial f
	m = (d+1)*t+1;    // number of servers
    flint_randinit(state);

	fmpz_init(p);
    fmpz_set_str(p,"340282366920938463463374607431768211297",10);  // 128 bit prime
    // fmpz_set_str(p,"17",10);
    fmpz_mod_ctx_init(univarctx,p);
    fmpz_mod_mpoly_ctx_init(nvarsctx,n,ORD_LEX,p);
}



// ************ Share algorithm ********************
void Share(fmpz_t * x, fmpz_t ** s){
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
void Eval(int j,fmpz_mod_mpoly_t f,fmpz_t * s_j,fmpz_t out_j){
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
void H_interpolation(fmpz_t * locs, fmpz_t * vals, fmpz_t * dif_vals, fmpz_t res){

}
int Ver(fmpz_t * out, fmpz_t y){
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


int main(){
	init();  // Initialize public parameters
	
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
    printf("***Sharing***\n");
	Share(x,in);

	// Eval
	printf("***Evaling***\n");	
	for (int j = 1;j <= m; j ++){
        Eval(j,f,in[j],out[j]);
    }

    // Ver
    printf("***Verifying***\n");
    int IsCorrect;
    fmpz_t y;
    fmpz_init(y);
    IsCorrect = Ver(out,y);



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

    return 0;
}
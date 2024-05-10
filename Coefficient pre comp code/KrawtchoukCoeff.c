/*
    Author: R.A. Dwars
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3.0 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.
*/

#include "math.h"
#include "stdlib.h"
#include "flint/acb_poly.h"
#include "flint/acb_mat.h"
#include "flint/acb_calc.h"
#include "flint/acb_dirichlet.h"
#include "flint/arb_hypgeom.h"
#include "flint/acb_hypgeom.h"
#include "flint/profiler.h"

//Xi(1/2 + x*I)
void
Xi(arb_t result, arb_t x, slong prec)
{  
	arb_t ar;
	arb_init(ar);
	
	acb_t ac;
	acb_init(ac);
	
	arb_one(ar);
 	arb_mul_2exp_si(ar, ar, -1);
	acb_set_arb_arb(ac, ar, x);
	acb_dirichlet_xi(ac, ac, prec);
	acb_get_real(result, ac);
	
	arb_clear(ar);
	acb_clear(ac);
}

//sin(x/2)*2/x
void
Xin(arb_t result, arb_t x, slong prec)
{  
	arb_t ar, br;
	arb_init(ar);
	arb_init(br);
	
	//to avoid division by zero:
	arb_set_str(ar, "10e-400", prec);
	arb_add(ar, x, ar, prec);
	
	arb_mul_2exp_si(br, ar, -1);
	arb_sin(br, br, prec);
	arb_div(br, br, ar, prec);
	arb_mul_2exp_si(result, br, 1);
	
	arb_clear(ar);
	arb_clear(br);
}

// Compute the weight
void
wei(arb_t result, arb_t x, slong N, arb_t p, slong prec)
{  
	arb_t ar, br, cr, Nr;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	arb_init(Nr);
	
	arb_set_si(Nr, N);
	
	//GAMMA(N + 1)*p^x*(1 - p)^(N - x)/(GAMMA(x + 1)*GAMMA(N - x + 1))
	arb_add_si(ar, Nr, 1, prec);
	arb_gamma(ar, ar, prec);
	
	arb_pow(br, p, x, prec);
	arb_mul(ar, ar, br, prec);
	
	arb_sub_si(br, p, 1, prec);
	arb_neg(br, br);
	arb_sub(cr, Nr, x, prec);
	arb_pow(br, br, cr, prec);
	arb_mul(ar, ar, br, prec);

	arb_add_si(br, x, 1, prec);
	arb_gamma(br, br, prec);
	
	arb_add_si(cr, Nr, 1, prec);
	arb_sub(cr, cr, x, prec);
	arb_gamma(cr, cr, prec);
	arb_mul(br, br, cr, prec);
	
	arb_div(result, ar, br, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
	arb_clear(Nr);
}

// Compute the hypergeometric function
void
PolK(arb_t result, arb_t x, slong N, arb_t pr, arb_t n, slong prec)
{  
	arb_t ar, br, Nr, z;
	arb_init(ar);
	arb_init(br);
	arb_init(Nr);
	arb_init(z);
	
	arb_struct p[2];
	arb_init(p);
	arb_init(p + 1);
	
	arb_struct q[1];
	arb_init(q);
	
	arb_set_si(Nr, N);
	arb_set_str(ar, "10e-200", prec);
	arb_add(Nr, Nr, ar, prec);
	
	//fill the p values
	arb_neg(p, n);
	arb_neg(p + 1, x);
	
	//fill the q values
	arb_neg(q, Nr);
	
	//fill z value
	arb_inv(z, pr, prec);
	
	//hypergeom([-n, -x], [-N], 1/p)
	arb_hypgeom_pfq(result, p, 2, q, 1, z, 0, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(Nr);
	arb_clear(p);
	arb_clear(p + 1);
	arb_clear(q);
	arb_clear(z);
}

// Compute the nfactor
void
nfactor(arb_t result, slong N, arb_t p, arb_t n, slong prec)
{  
	arb_t ar, br, cr, Nr;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	arb_init(Nr);
	
	arb_set_si(Nr, N);
	
	//(-1)^(-2*n)*p^n*(-1 + p)^(-n)*pochhammer(-N, n)/n!
	arb_mul_2exp_si(ar, n, 1);
	arb_neg(ar, ar);
	arb_set_si(br, -1);
	arb_pow(ar, br, ar, prec);
	
	arb_pow(br, p, n, prec);
	arb_mul(ar, ar, br, prec);
	
	arb_sub_si(br, p, 1, prec);
	arb_neg(cr, n);
	arb_pow(br, br, cr, prec);
	arb_mul(ar, ar, br, prec);

	arb_neg(br, Nr);
	arb_hypgeom_rising(br, br, n, prec);
	arb_mul(ar, ar, br, prec);	
	
	arb_add_si(br, n, 1, prec);
	arb_gamma(br, br, prec);
	arb_div(result, ar, br, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
	arb_clear(Nr);
}

// Prepare the integrand Xi(x)*hyp(x,a,b,c,n)*wei(x,a,b,c)
void
summand(arb_t result, arb_t x, arb_t n, arb_t p, slong N, slong sw, slong prec)
{
	arb_t ar, br;
	arb_init(ar);
	arb_init(br);
	
	acb_t ac;
	acb_init(ac);

	//select the desired function to expand
    switch(sw)
    {
        case 0:
            Xi(ar, x, prec);
            break;

        case 1:
            Xin(ar, x, prec);
            break;

        default:
            break;
    }

	PolK(br, x, N, p, n, prec);
	arb_mul(ar, ar, br, prec);
	
	wei(br, x, N, p, prec);
	arb_mul(result, ar, br, prec);

	arb_clear(ar);	
	arb_clear(br);

	acb_clear(ac);
}

// Evaluate the integral and multiply output by the nfactor
void
sum_coeff(arb_t result, slong N, arb_t p, arb_t n, slong sw, slong prec)
{  
	arb_t ar, br, sumtot, x;
	arb_init(ar);
	arb_init(br);
	arb_init(sumtot);
	arb_init(x);
	
	slong i;
	
	arb_zero(sumtot);
	for (i = 0; i <= N; i++)
	{
		arb_set_si(x, i);
		summand(ar, x, n, p, N, sw, prec);
		arb_add(sumtot, sumtot, ar, prec);
	}
	
	//res * nfactor(N, a, n)
	nfactor(ar, N, p, n, prec);
	arb_mul(result, sumtot, ar, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(sumtot);
	arb_clear(x);
}

int main(int argc, char *argv[])
{
    arb_t p, ar, n, x, res;
    arb_init(p);
    arb_init(ar);
    arb_init(n);
    arb_init(x);
    arb_init(res);

	acb_t ac;
	acb_clear(ac);

	slong i, N, numthreads, prec, sw;

	int result = EXIT_SUCCESS;
	
	if (argc != 4)
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	//set precision and the integral upper bound
	prec = 30 * 3.32192809488736 + 1500;
	
	arb_set_str(p, argv[1], prec);
    numthreads = atol(argv[2]);
    sw = atol(argv[3]);
	
	arb_one(res);
	if ((sw < 0) || (sw > 1) || (numthreads < 1) || (arb_contains_nonpositive(p)) || (arb_ge(p, res)))
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	flint_set_num_threads(numthreads);
	
	N = 200;

	//loop through n = 0 to N and compute each coefficient
	for (i = 0; i <= N; i++)
    {
		arb_set_si(n, i);
		sum_coeff(res, N, p, n, sw, prec);
		arb_printn(res, 120, ARB_STR_NO_RADIUS);
		printf("\n");
    }

finish:

    if (result == EXIT_FAILURE)
    {
		flint_printf("\n");
		flint_printf("Required inputs:\n");
		flint_printf("%s p numthreads sw \n\n", argv[0]);
		flint_printf(
    "This script computes the first 250 coefficients for the expansion of the Riemann Xi-function\n"
    "into the orthogonal Dual Hahn polynomials as described in the paper:\n"
    ">> Visualising the flows of orthogonal polynomial expansions of the Riemann Xi-function <<\n"
    "-> 'numthreads' defines the number of desired threads/available cores to increase speed.\n"
    "-> 'p' is the required parameter that requires 0 < p < 1 \n"
    "-> 'sw' allows to select the function to expand: Xi(1/2+xi) (sw=0) or 2/x*sin(x/2) (sw=1)\n"
    "The user can push the output into a file by adding: >yyycoeffparms.txt \n");
    }

    arb_clear(p);
    arb_clear(ar);
    arb_clear(n);
    arb_clear(x);
    arb_clear(res);

	acb_clear(ac);
 
    flint_cleanup();

    return 0;  
} 


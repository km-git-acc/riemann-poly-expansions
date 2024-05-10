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
wei(arb_t result, arb_t x, arb_t b, arb_t c, slong prec)
{  
	arb_t ar, br;
	arb_init(ar);
	arb_init(br);
	
	//pochhammer(b, x)*c^x/x!
	arb_hypgeom_rising(ar, b, x, prec);
	
	arb_pow(br, c, x, prec);
	arb_mul(ar, ar, br, prec);
	
	arb_add_si(br, x, 1, prec);
	arb_gamma(br, br, prec);
	arb_div(result, ar, br, prec);

	arb_clear(ar);
	arb_clear(br);
}

// Compute the hypergeometric function
void
PolM(arb_t result, arb_t x, arb_t b, arb_t c, arb_t n, slong prec)
{  
	arb_t ar, br, z;
	arb_init(ar);
	arb_init(br);
	arb_init(z);
	
	arb_struct p[2];
	arb_init(p);
	arb_init(p + 1);
	
	arb_struct q[1];
	arb_init(q);
	
	//fill the p values
	arb_neg(p, n);
	arb_neg(p + 1, x);
	
	//fill the q values
	arb_set(q, b);
	
	//fill z value
	arb_one(ar);
	arb_inv(br, c, prec);
	arb_sub(z, ar, br, prec);
	
	//hypergeom([-n, -x], [b], 1 - 1/c)
	arb_hypgeom_pfq(result, p, 2, q, 1, z, 0, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(p);
	arb_clear(p + 1);
	arb_clear(q);
	arb_clear(z);
}

// Compute the nfactor
void
nfactor(arb_t result, arb_t b, arb_t c, arb_t n, slong prec)
{  
	arb_t ar, br;
	arb_init(ar);
	arb_init(br);
	
	//c^n*pochhammer(b, n)*(1 - c)^b/n!
	arb_pow(ar, c, n, prec);
	
	arb_hypgeom_rising(br, b, n, prec);
	arb_mul(ar, ar, br, prec);
	
	arb_sub_si(br, c, 1, prec);
	arb_neg(br, br);
	arb_pow(br, br, b, prec);
	arb_mul(ar, ar, br, prec);
	
	arb_add_si(br, n, 1, prec);
	arb_gamma(br, br, prec);
	arb_div(result, ar, br, prec);
	
	arb_clear(ar);
	arb_clear(br);
}

// Prepare the integrand Xi(x)*hyp(x,a,b,c,n)*wei(x,a,b,c)
void
summand(arb_t result, arb_t x, arb_t n, arb_t b, arb_t c, slong sw, slong prec)
{
	arb_t ar, br;
	arb_init(ar);
	arb_init(br);

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

	PolM(br, x, b, c, n, prec);
	arb_mul(ar, ar, br, prec);
	
	wei(br, x, b, c, prec);
	arb_mul(result, ar, br, prec);

	arb_clear(ar);	
	arb_clear(br);
}

// Evaluate the integral and multiply output by the nfactor
void
sum_coeff(arb_t result, arb_t b, arb_t c, arb_t n, slong sw, slong prec)
{  
	arb_t ar, br, sumtot, x;
	arb_init(ar);
	arb_init(br);
	arb_init(sumtot);
	arb_init(x);
	
	slong i;
	
	arb_zero(sumtot);
	for (i = 0; i <= 300; i++)
	{
		arb_set_si(x, i);
		summand(ar, x, n, b, c, sw, prec);
		arb_add(sumtot, sumtot, ar, prec);
	}
	
	//res * nfactor(N, a, n)
	nfactor(ar, b, c, n, prec);
	arb_mul(result, sumtot, ar, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(sumtot);
	arb_clear(x);
}

int main(int argc, char *argv[])
{
    arb_t b, c, ar, n, x, res;
    arb_init(b);
    arb_init(c);
    arb_init(ar);
    arb_init(n);
    arb_init(x);
    arb_init(res);

	acb_t ac;
	acb_clear(ac);

	slong i, N, numthreads, prec, sw;

	int result = EXIT_SUCCESS;
	
	if (argc != 5)
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	//set precision and the integral upper bound
	prec = 30 * 3.32192809488736 + 1500;
	
	arb_set_str(b, argv[1], prec);
	arb_set_str(c, argv[2], prec);
    numthreads = atol(argv[3]);
    sw = atol(argv[4]);

	if ((sw < 0) || (sw > 1) || (numthreads < 1))
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	flint_set_num_threads(numthreads);
	
	N = 400;

	//loop through n = 0 to N and compute each coefficient
	for (i = 0; i <= N; i++)
    {
		arb_set_si(n, i);
		sum_coeff(res, b, c, n, sw, prec);
		arb_printn(res, 120, ARB_STR_NO_RADIUS);
		printf("\n");
    }

finish:

    if (result == EXIT_FAILURE)
    {
		flint_printf("\n");
		flint_printf("Required inputs:\n");
		flint_printf("%s b c numthreads sw \n\n", argv[0]);
		flint_printf(
    "This script computes the first 250 coefficients for the expansion of the Riemann Xi-function\n"
    "into the orthogonal Dual Hahn polynomials as described in the paper:\n"
    ">> Visualising the flows of orthogonal polynomial expansions of the Riemann Xi-function <<\n"
    "-> 'numthreads' defines the number of desired threads/available cores to increase speed.\n"
    "-> 'b, c' are the required parameters.\n"
    "-> 'sw' allows to select the function to expand: Xi(1/2+xi) (sw=0) or 2/x*sin(x/2) (sw=1)\n"
    "The user can push the output into a file by adding: >yyycoeffparms.txt \n");
    }

    arb_clear(b);
    arb_clear(c);
    arb_clear(ar);
    arb_clear(n);
    arb_clear(x);
    arb_clear(res);

	acb_clear(ac);
 
    flint_cleanup();

    return 0;  
} 



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

// Compute tranformation factor
void
trans(arb_t result, arb_t x, arb_t c, arb_t d, slong prec)
{  
	//-c/2 - d/2 - 1/2 - sqrt(c^2 + 2*c*d + d^2 + 2*c + 2*d + 4*y + 1)/2
	arb_t ar, br;
	arb_init(ar);
	arb_init(br);

	arb_sqr(ar, c, prec);
	arb_mul(br, c, d, prec);
 	arb_mul_2exp_si(br, br, 1);
	arb_add(ar, ar, br, prec);
	
	arb_sqr(br, d, prec);
	arb_add(ar, ar, br, prec);
	
	arb_mul_2exp_si(br, c, 1);
	arb_add(ar, ar, br, prec);
	
	arb_mul_2exp_si(br, d, 1);
	arb_add(ar, ar, br, prec);
	
	arb_mul_2exp_si(br, x, 2);
	arb_add(ar, ar, br, prec);
	arb_add_si(ar, ar, 1, prec);
	
	arb_sqrt(ar, ar, prec);
	arb_mul_2exp_si(ar, ar, -1);
	
	arb_mul_2exp_si(br, c, -1);
	arb_add(ar, ar, br, prec);
	arb_mul_2exp_si(br, d, -1);
	arb_add(ar, ar, br, prec);
	arb_one(br);
	arb_mul_2exp_si(br, br, -1);
	arb_add(ar, ar, br, prec);
	
	arb_neg(ar, ar);

	//a slight augmentation of N is required when x = N
	arb_set_str(br, "10e-120", prec);
	arb_add(result, ar, br, prec);

	arb_clear(ar);
	arb_clear(br);
}

// Compute tranformation factor inverse
void
itrans(arb_t result, arb_t x, arb_t c, arb_t d, slong prec)
{  
	//x*(x + c + d + 1)
	arb_t ar;
	arb_init(ar);

	arb_add(ar, x, c, prec);
	arb_add(ar, ar, d, prec);
	arb_add_si(ar, ar, 1, prec);
	
	arb_mul(result, ar, x, prec);

	arb_clear(ar);
}
// Compute the weight
void
wei(arb_t result, arb_t x, slong N, arb_t c, arb_t d, slong prec)
{  
	arb_t ar, br, cr, dr, Nr;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	arb_init(dr);
	arb_init(Nr);
	
	arb_set_si(Nr, N);
	
	//(2*x + c + d + 1)*pochhammer(-N, x)*pochhammer(c + 1, x)*N!
	arb_add(ar, c, d, prec);
	arb_add_si(ar, ar, 1, prec);
	arb_mul_2exp_si(br, x, 1);
	arb_add(ar, ar, br, prec);
	
	arb_neg(br, Nr);
	arb_hypgeom_rising(br, br, x, prec);
	arb_mul(ar, ar, br, prec);
	
	arb_add_si(br, c, 1, prec);
	arb_hypgeom_rising(br, br, x, prec);
	arb_mul(ar, ar, br, prec);

	arb_add_si(br, Nr, 1, prec);
	arb_gamma(br, br, prec);
	arb_mul(ar, ar, br, prec);
	
	//(-1)^x*pochhammer(x + c + d + 1, N + 1)*pochhammer(d + 1, x)*x!
	arb_set_si(br, -1);
	arb_pow(br, br, x, prec);
 
	arb_add(cr, c, d, prec);
	arb_add(cr, cr, x, prec);
	arb_add_si(cr, cr, 1, prec);
	arb_add_si(dr, Nr, 1, prec);
	arb_hypgeom_rising(cr, cr, dr, prec);
	arb_mul(br, br, cr, prec);
	
	arb_add_si(cr, d, 1, prec);
	arb_hypgeom_rising(cr, cr, x, prec);
	arb_mul(br, br, cr, prec);
	
	arb_add_si(cr, x, 1, prec);
	arb_gamma(cr, cr, prec);
	arb_mul(br, br, cr, prec);
	
	//divide both of the above
	arb_div(result, ar, br, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
	arb_clear(dr);
	arb_clear(Nr);
}

// Compute the hypergeometric function
void
PolR(arb_t result, arb_t x, slong N, arb_t c, arb_t d, arb_t n, slong prec)
{  
	arb_t ar, br, Nr, z;
	arb_init(ar);
	arb_init(br);
	arb_init(Nr);
	arb_init(z);
	
	arb_struct p[3];
	arb_init(p);
	arb_init(p + 1);
	arb_init(p + 2);
	
	arb_struct q[2];
	arb_init(q);
	arb_init(q + 1);
	
	//a slight augmentation of N is required when n = N
	arb_set_si(Nr, N);
	arb_set_str(ar, "10e-120", prec);
	arb_add(Nr, Nr, ar, prec);
	
	//fill the p values
	arb_neg(p, n);
	
	trans(br, x, c, d, prec);
	arb_neg(p + 1, br);
	
	arb_add(ar, c, d, prec);
	arb_add(ar, ar, br, prec);
	arb_add_si(ar, ar, 1, prec);
	arb_set(p + 2, ar);
	
	//fill the q values
	arb_add_si(q, c, 1, prec);
	
	arb_neg(ar, Nr);
	arb_set(q + 1, ar);
	
	//fill z value
	arb_one(z);
	
	//hypergeom([-n, -trans(x, c, d), trans(x, c, d) + c + d + 1], [c + 1, -N], 1)
	arb_hypgeom_pfq(result, p, 3, q, 2, z, 0, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(Nr);
	arb_clear(p);
	arb_clear(p + 1);
	arb_clear(p + 2);
	arb_clear(q);
	arb_clear(q + 1);
	arb_clear(z);
}

// Compute the nfactor
void
nfactor(arb_t result, slong N, arb_t c, arb_t d, arb_t n, slong prec)
{  
	arb_t ar, br, cr, Nr;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	arb_init(Nr);
	
	arb_set_si(Nr, N);
	
	//GAMMA(n + c + 1)*GAMMA(d - n + N + 1)
	arb_add(ar, n, c, prec);
	arb_add_si(ar, ar, 1, prec);
	arb_gamma(ar, ar, prec);
	
	arb_sub(br, d, n, prec);
	arb_add(br, br, Nr, prec);
	arb_add_si(br, br, 1, prec);
	arb_gamma(br, br, prec);
	arb_mul(ar, ar, br, prec);

	//GAMMA(c + 1)*GAMMA(d + 1)*GAMMA(n + 1)*GAMMA(-n + N + 1)
	arb_add_si(br, c, 1, prec);
	arb_gamma(br, br, prec);
	
	arb_add_si(cr, d, 1, prec);
	arb_gamma(cr, cr, prec);
	arb_mul(br, br, cr, prec);

	arb_add_si(cr, n, 1, prec);
	arb_gamma(cr, cr, prec);
	arb_mul(br, br, cr, prec);
	
	arb_add_si(cr, Nr, 1, prec);
	arb_sub(cr, cr, n, prec);
	arb_gamma(cr, cr, prec);
	arb_mul(br, br, cr, prec);
	
	//divide both of the above
	arb_div(result, ar, br, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
	arb_clear(Nr);
}

// Prepare the integrand Xi(x)*hyp(x,a,b,c,n)*wei(x,a,b,c)
void
summand(arb_t result, arb_t x, arb_t n, arb_t c,arb_t d, slong N, slong sw, slong prec)
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

	itrans(br, x, c, d, prec);
	PolR(br, br, N, c, d, n, prec);
	arb_mul(ar, ar, br, prec);
	
	wei(br, x, N, c, d, prec);
	arb_mul(result, ar, br, prec);

	arb_clear(ar);	
	arb_clear(br);

	acb_clear(ac);
}

// Evaluate the integral and multiply output by the nfactor
void
sum_coeff(arb_t result, slong N, arb_t c, arb_t d, arb_t n, slong sw, slong prec)
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
		summand(ar, x, n, c, d, N, sw, prec);
		arb_add(sumtot, sumtot, ar, prec);
	}
	
	//res * nfactor(N, a, n)
	nfactor(ar, N, c, d, n, prec);
	arb_mul(result, sumtot, ar, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(sumtot);
	arb_clear(x);
}

int main(int argc, char *argv[])
{
    arb_t c, d, ar, n, x, res;
    arb_init(c);
    arb_init(d);
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
	
	arb_set_str(c, argv[1], prec);
	arb_set_str(d, argv[2], prec);
	numthreads = atol(argv[3]);
	sw = atol(argv[4]);
	
	N = 125;
	
	if ((sw < 0) || (sw > 1) || (numthreads < 1))
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	flint_set_num_threads(numthreads);

	//loop through n = 0 to N and compute each coefficient
	for (i = 0; i <= N; i++)
    {
		arb_set_si(n, i);
		sum_coeff(res, N, c, d, n, sw, prec);
		arb_printn(res, 120, ARB_STR_NO_RADIUS);
		printf("\n");
    }

finish:

    if (result == EXIT_FAILURE)
    {
		flint_printf("\n");
		flint_printf("Required inputs:\n");
		flint_printf("%s c d numthreads sw \n\n", argv[0]);
		flint_printf(
    "This script computes the first 250 coefficients for the expansion of the Riemann Xi-function\n"
    "into the orthogonal Dual Hahn polynomials as described in the paper:\n"
    ">> Visualising the flows of orthogonal polynomial expansions of the Riemann Xi-function <<\n"
    "-> 'numthreads' defines the number of desired threads/available cores to increase speed.\n"
    "-> 'c, d' are the required parameters. Both c and d must be > -1 OR c and d < -N -1. \n"
    "-> 'sw' allows to select the function to expand: Xi(1/2+xi) (sw=0) or 2/x*sin(x/2) (sw=1)\n"
    "The user can push the output into a file by adding: >yyycoeffparms.txt \n");
    }

    arb_clear(c);
    arb_clear(d);
    arb_clear(ar);
    arb_clear(n);
    arb_clear(x);
    arb_clear(res);

	acb_clear(ac);
 
    flint_cleanup();

    return 0;  
} 


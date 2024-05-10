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
	arb_set_str(ar, "10e-500", prec);
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
wei(arb_t result, arb_t x, slong prec)
{  
	arb_t ar;
	arb_init(ar);
	
	//exp(-x^2)
	arb_sqr(ar, x, prec);
	arb_neg(ar, ar);
	arb_exp(result, ar, prec);
	
	arb_clear(ar);
}

// Compute the hypergeometric function
void
hyp(arb_t result, arb_t x, arb_t n, slong prec)
{  
	arb_t ar, br, cr;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	
	arb_t z;
	arb_init(z);
	
	//fill the p values
 	arb_mul_2exp_si(ar, n, -1);
	arb_neg(ar, ar);
	
	//fill the q values
	arb_one(br);
 	arb_mul_2exp_si(br, br, -1);
	
	//fill z value
	arb_sqr(z, x, prec);
	
	//KummerU(-n/2, 1/2, x^2)
	arb_hypgeom_1f1(result, ar, br, z, 0, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
	arb_clear(z);
}

// Compute the nfactor
void
nfactor(arb_t result, arb_t n, slong prec)
{  
	arb_t ar, br, cr, pi;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	arb_init(pi);

	arb_const_pi(pi, prec);
	
	arb_mul_2exp_si(cr, n, -1);
	
	//gammma(n/2+1/2)/(pi*gamma(n/2+1))
	arb_one(ar);
	arb_mul_2exp_si(ar, ar, -1);
	arb_add(ar, ar, cr, prec);
	arb_gamma(ar, ar, prec);
	
	arb_add_si(br, cr, 1, prec);
	arb_gamma(br, br, prec);
	arb_mul(br, br, pi, prec);
	
	arb_div(result, ar, br, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
	arb_clear(pi);
}

// Prepare the integrand Xi*hyp*wei
int
f_integrand(acb_ptr result, const acb_t xc, void *param, slong order, slong prec)
{
	arb_t ar, br, x;
	arb_init(ar);
	arb_init(br);
	arb_init(x);
	
	acb_t ac;
	acb_init(ac);

	arb_ptr(n);
	
	slong sw;

	n = ((arb_ptr *)(param))[0];
	sw = ((slong *)(param))[1];
	
	acb_get_real(x, xc);

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

	//compute the hypergeom
	hyp(br, x, n, prec);
	arb_mul(ar, ar, br, prec);
	
	//compute the weight
	wei(br, x, prec);
	arb_mul(ar, ar, br, prec);

	acb_set_arb(result, ar);

	arb_clear(ar);
	arb_clear(br);	
	arb_clear(x);

	acb_clear(ac);

    return 0;
}

// Evaluate the integral and multiply output by the nfactor
void
integral_coeff(arb_t result, arb_t n, slong intlim, slong sw, slong prec)
{  
	arb_t ar, br, pi;
	arb_init(ar);
	arb_init(br);
	arb_init(pi);
	
	acb_t ac, ai, bi;
	acb_init(ac);
	acb_init(ai);
	acb_init(bi);

	arb_const_pi(pi, prec);

	slong goal;

    //evaluate the integrals
	acb_calc_integrate_opt_t options;
	acb_calc_integrate_opt_init(options);

	mag_t tol;
	mag_init(tol);

	goal = prec;
	mag_set_ui_2exp_si(tol, 1, -prec);

	//set limits of integration, note that the lower limit must be 0+
	arb_zero(ar);
	arb_set_si(br, intlim);
	acb_set_arb(ai, ar);
	acb_set_arb(bi, br);

	void *param[2];
	param[0] = (void *) n;
	param[1] = (void *) sw;

	//evaluate integral
	acb_calc_integrate(ac, f_integrand, param, ai, bi, goal, tol, NULL, prec);
		
	acb_get_real(result, ac);
	
	//res * nfactor(n) * 2
	nfactor(ar, n, prec);
 	arb_mul_2exp_si(ar, ar, 1);
	arb_mul(result, result, ar, prec);
	
	mag_clear(tol);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(pi);

	acb_clear(ac);
	acb_clear(ai);
	acb_clear(bi);
}

int main(int argc, char *argv[])
{
    arb_t a, n, res;
    arb_init(a);
    arb_init(n);
    arb_init(res);

	slong i, prec, numthreads, intlim, sw;

	int result = EXIT_SUCCESS;
	
	if (argc != 3)
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	//set precision and the integral upper bound
	prec = 30 * 3.32192809488736 + 2500;
	intlim = 80;

    numthreads = atol(argv[1]);
    sw = atol(argv[2]);
	
	if ((sw < 0) || (sw > 1) || (numthreads < 1))
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	flint_set_num_threads(numthreads);
	
	//loop through n = 0 to 200 and compute each coefficient
	for (i = 0; i <= 500; i++)
    {
		arb_set_si(n, i);
				//no need to compute the coefficient for n=odd
		if (i % 2)
			arb_zero(res);
		else 
			integral_coeff(res, n, intlim, sw, prec);
		arb_printn(res, 120, ARB_STR_NO_RADIUS);
		printf("\n");
    }

finish:

    if (result == EXIT_FAILURE)
    {
		flint_printf("\n");
		flint_printf("Required inputs:\n");
		flint_printf("%s numthreads sw \n\n", argv[0]);
		flint_printf(
    "This script computes the first 200 coefficients for the expansion of the Riemann Xi-function\n"
    "into the orthogonal Hermite polynomials as described in the paper:\n"
    ">> Visualising the flows of orthogonal polynomial expansions of the Riemann Xi-function <<\n"
    "-> 'numthreads' defines the number of desired threads/available cores to increase speed.\n"
    "-> 'sw' allows to select the function to expand: Xi(1/2+xi) (sw=0) or 2/x*sin(x/2) (sw=1)\n"
    "The user can push the output into a file by adding: >yyycoeffparms.txt \n");
    }

    arb_clear(a);
    arb_clear(n);
    arb_clear(res);
 
    flint_cleanup();

    return 0;  
} 

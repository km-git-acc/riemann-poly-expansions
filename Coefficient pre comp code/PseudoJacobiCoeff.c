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
	arb_set(ar, x);
	if (arb_contains_zero(ar))
	{
		arb_set_str(br, "10e-200", prec);
		arb_add(ar, ar, br, prec);
	}
	
	arb_mul_2exp_si(br, ar, -1);
	arb_sin(br, br, prec);
	arb_div(br, br, ar, prec);
	arb_mul_2exp_si(result, br, 1);
	
	arb_clear(ar);
	arb_clear(br);
}

// Compute the weight
void
wei(arb_t result, arb_t x, slong N, arb_t a, slong prec)
{  
	arb_t ar, br, cr;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	
	//(1+x^2)^(-N-1)*exp(2*a*arctan(x)
	arb_sqr(ar, x, prec);
	arb_add_si(ar, ar, 1, prec);
	arb_set_si(br, -N-1);
	arb_pow(ar, ar, br, prec);
	
 	arb_mul_2exp_si(br, a, 1);
	arb_atan(cr, x, prec);
	arb_mul(br, br, cr, prec);
	arb_exp(br, br, prec);
	
	arb_mul(result, ar, br, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
}

// Compute the hypergeometric function
void
hyp(acb_t result, arb_t x, slong N, arb_t a, arb_t n, slong prec)
{  
	arb_t ar, br;
	arb_init(ar);
	arb_init(br);
	
	acb_t q, z;
	acb_init(q);
	acb_init(z);
	
	acb_struct p[2];
	acb_init(p);
	acb_init(p + 1);
	
	//fill the p values
	arb_neg(ar, n);
	acb_set_arb(p, ar);
	arb_sub_si(ar, n, 2*N+1, prec);
	acb_set_arb(p + 1, ar);
	
	//fill the q value
	arb_set_si(ar, -N);
	acb_set_arb_arb(q, ar, a);
	
	//fill z value
	arb_one(ar);
	arb_mul_2exp_si(ar, ar, -1);
	arb_neg(br, x);
	arb_mul_2exp_si(br, br, -1);
	acb_set_arb_arb(z, ar, br);
	
	//Hypergeom2f1([-n,n-2*N-1],[-N+I*a],(1-I*x)/2)
	acb_hypgeom_pfq(result, p, 2, q, 1, z, 0, prec);
	
	arb_clear(ar);
	arb_clear(br);
	
	acb_clear(p);
	acb_clear(p + 1);
	acb_clear(q);
	acb_clear(z);
}

// Compute the nfactor
void
nfactor(acb_t result, slong N, arb_t a, arb_t n, slong prec)
{  
	arb_t ar, br, cr, pi;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	arb_init(pi);

	acb_t ac, bc;
	acb_init(ac);
	acb_init(bc);
	
	arb_const_pi(pi, prec);
	
	//4^N*(2*N+1-2*n)*gamma(n-N+a*I)^2*|gamma(N+1-n+a*I)|^2*(-1)^n
	arb_set_si(ar, 4);
	arb_set_si(br, N);
	arb_pow(ar, ar, br, prec);
	
	arb_set_si(br, 2*N+1);
	arb_mul_2exp_si(cr, n, 1);
	arb_sub(br, br, cr, prec);
	arb_mul(ar, ar, br, prec);	
	
	arb_set_si(br, N);
	arb_sub(br, n, br, prec);
	acb_set_arb_arb(ac, br, a);
	acb_gamma(ac, ac, prec);
	acb_sqr(ac, ac, prec);
	
	arb_set_si(br, N + 1);
	arb_sub(br, br, n, prec);
	acb_set_arb_arb(bc, br, a);
	acb_gamma(bc, bc, prec);
	acb_abs(br, bc, prec);
	arb_sqr(br, br, prec);
	arb_mul(ar, ar, br, prec);
	
	arb_set_si(br, -1);
	arb_pow(br, br, n, prec);
	arb_mul(ar, ar, br, prec);
	
	acb_mul_arb(ac, ac, ar, prec);
	
	//Pi*gamma(n+1)*gamma(2*N+2-n)*gamma(-N+a*I)^2
	arb_add_si(ar, n, 1, prec);
	arb_gamma(ar, ar, prec);
	arb_mul(ar, ar, pi, prec);
	
	arb_set_si(br, 2*N+2);
	arb_sub(br, br, n, prec);
	arb_gamma(br, br, prec);
	arb_mul(ar, ar, br, prec);
	
	arb_set_si(br, -N);
	acb_set_arb_arb(bc, br, a);
	acb_gamma(bc, bc, prec);
	acb_sqr(bc, bc, prec);
	
	acb_mul_arb(bc, bc, ar, prec);
	
	//divide both of the above
	acb_div(result, ac, bc, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
	arb_clear(pi);
	
	acb_clear(ac);
	acb_clear(bc);
}

// Prepare the integrand Xi(x)*hyp(x,a,b,c,n)*wei(x,a,b,c)
int
f_integrand(acb_ptr result, const acb_t xc, void *param, slong order, slong prec)
{
	arb_t ar, br, x;
	arb_init(ar);
	arb_init(br);
	arb_init(x);
	
	acb_t ac, bc;
	acb_init(ac);
	acb_init(bc);

	arb_ptr(a);
	arb_ptr(n);

	slong N, sw;

	N = ((slong *)(param))[0];
	a = ((arb_ptr *)(param))[1];
	n = ((arb_ptr *)(param))[2];
	sw = ((slong *)(param))[3];
	
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
	acb_set_arb(ac, ar);

	hyp(bc, x, N, a, n, prec);
	acb_mul(ac, ac, bc, prec);
	
	wei(ar, x, N, a, prec);
	acb_mul_arb(result, ac, ar, prec);

	arb_clear(ar);	
	arb_clear(br);	
	arb_clear(x);

	acb_clear(ac);
	acb_clear(bc);

    return 0;
}

// Evaluate the integral and multiply output by the nfactor
void
integral_coeff(acb_t result, slong N, arb_t a, arb_t n, slong intlim, slong sw, slong prec)
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

	//set limits of integration
	arb_set_si(ar, -intlim);
	arb_set_si(br, intlim);
	acb_set_arb(ai, ar);
	acb_set_arb(bi, br);

	void *param[4];
	param[0] = (void *) N;
	param[1] = (void *) a;
	param[2] = (void *) n;
	param[3] = (void *) sw;

	//evaluate integral
	acb_calc_integrate(result, f_integrand, param, ai, bi, goal, tol, NULL, prec);
	
	//res * nfactor(N, a, n)
	nfactor(ac, N, a, n, prec);
	acb_mul(result, result, ac, prec);
	
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
    arb_t a, ar, n, x;
    arb_init(a);
    arb_init(ar);
    arb_init(n);
    arb_init(x);
	
	acb_t res;
    acb_init(res);

	slong i, prec, N, numthreads, intlim, sw;

	int result = EXIT_SUCCESS;
	
	if (argc != 4)
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	//set precision and the integral upper bound
	prec = 30 * 3.32192809488736 + 2500;
	intlim = 200;
	
	arb_set_str(a, argv[1], prec);
    numthreads = atol(argv[2]);
    sw = atol(argv[3]);
	
	if ((sw < 0) || (sw > 1) || (numthreads < 1))
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	flint_set_num_threads(numthreads);
	
	N = 300;

	//loop through n = 0 to N and compute each coefficient
	for (i = 0; i <= N; i++)
    {
		arb_set_si(n, i);
		integral_coeff(res, N, a, n, intlim, sw, prec);
		acb_get_real(ar, res);
		arb_printn(ar, 60, ARB_STR_NO_RADIUS);
		printf(", ");
		acb_get_imag(ar, res);
		arb_printn(ar, 60, ARB_STR_NO_RADIUS);
		printf("\n");
    }

finish:

    if (result == EXIT_FAILURE)
    {
		flint_printf("\n");
		flint_printf("Required inputs:\n");
		flint_printf("%s a numthreads sw \n\n", argv[0]);
		flint_printf(
    "This script computes the first 250 coefficients for the expansion of the Riemann Xi-function\n"
    "into the orthogonal Pseudo Jacobi polynomials as described in the paper:\n"
    ">> Visualising the flows of orthogonal polynomial expansions of the Riemann Xi-function <<\n"
    "-> 'numthreads' defines the number of desired threads/available cores to increase speed.\n"
    "-> 'a' is the required parameter. Note that N is set in the code above.\n"
    "-> 'sw' allows to select the function to expand: Xi(1/2+xi) (sw=0) or 2/x*sin(x/2) (sw=1)\n"
    "The user can push the output into a file by adding: >yyycoeffparms.txt \n");
    }

    arb_clear(a);
    arb_clear(ar);
    arb_clear(n);
    arb_clear(x);
	
    acb_clear(res);
 
    flint_cleanup();

    return 0;  
} 


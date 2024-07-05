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
		arb_set_str(br, "10e-700", prec);
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
wei(arb_t result, acb_t xc, arb_t a, arb_t b, slong prec)
{  
	acb_t ac, bc, cc;
	acb_init(ac);
	acb_init(bc);
	acb_init(cc);
	
	//(1-x)^a
	acb_sub_si(ac, xc, 1, prec);
	acb_neg(ac, ac);
	acb_set_arb(cc, a);
	acb_pow_analytic(ac, ac, cc, 1, prec);

	//(1+x)^b
 	acb_add_si(bc, xc, 1, prec);
	acb_set_arb(cc, b);
	acb_pow_analytic(bc, bc, cc, 1, prec);
	
	acb_mul(ac, ac, bc, prec);
	acb_get_real(result, ac);

	acb_clear(ac);
	acb_clear(bc);	
	acb_clear(cc);	
}

// Compute the hypergeometric function
void
hyp(arb_t result, arb_t x, arb_t a, arb_t b, arb_t n, slong prec)
{  
	arb_t ar, z;
	arb_init(ar);
	arb_init(z);
	
	arb_struct p[2];
	arb_init(p);
	arb_init(p + 1);

	arb_struct q[1];
	arb_init(q);
	
	//fill the p values
	arb_neg(p + 0, n);
	
	arb_add_si(ar, n, 1, prec);
	arb_add(ar, ar, a, prec);
	arb_add(p + 1, ar, b, prec);
	
	//fill the q values
	arb_add_si(q, a, 1, prec);
	
	//fill z value
	arb_sub_si(ar, x, 1, prec);
	arb_neg(ar, ar);
	arb_mul_2exp_si(z, ar, -1);
	
	//hypergeom([-n, n + a + b + 1], [a + 1], (1 - x)/2)
	arb_hypgeom_pfq(result, p, 2, q, 1, z, 0, prec);
	
	arb_clear(ar);
	arb_clear(z);
	
	arb_clear(p);
	arb_clear(p + 1);

	arb_clear(q);
}

// Compute the nfactor
void
nfactor(arb_t result, arb_t a, arb_t b, arb_t n, slong prec)
{  
	arb_t ar, br, cr, a1, b1, ab1;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	arb_init(a1);
	arb_init(b1);
	arb_init(ab1);
	
	arb_add_si(a1, a, 1, prec);
	arb_add_si(b1, b, 1, prec);
	
	arb_add(ab1, a, b, prec);
	arb_add_si(ab1, ab1, 1, prec);

	//GAMMA(n + a + 1)*2^(-1 - a - b)*(2*n + a + b + 1)*GAMMA(n + a + b + 1)
	arb_add(ar, n, a1, prec);
	arb_gamma(ar, ar, prec);
	
	arb_add(br, n, ab1, prec);
	arb_gamma(br, br, prec);
	arb_mul(ar, ar, br, prec);
	
	arb_mul_2exp_si(br, n, 1);
	arb_add(br, br, ab1, prec);
	arb_mul(ar, ar, br, prec);
	
	arb_neg(br, ab1);
	arb_set_si(cr, 2);
	arb_pow(br, cr, br, prec);
	arb_mul(ar, ar, br, prec);
	
	//gamma(a+1)^2*gamma(n+1)*gamma(n+b+1)
	arb_gamma(br, a1, prec);
	arb_pow_ui(br, br, 2, prec);
	
	arb_add_si(cr, n, 1, prec);
	arb_gamma(cr, cr, prec);
	arb_mul(br, br, cr, prec);
	
	arb_add(cr, n, b1, prec);
	arb_gamma(cr, cr, prec);
	arb_mul(br, br, cr, prec);
	
	//divide both of the above
	arb_div(result, ar, br, prec);
	
	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
	arb_clear(a1);
	arb_clear(b1);
	arb_clear(ab1);
}

// Prepare the integrand Xi(x)*hyp(x,a,b,c,n)*wei(x,a,b,c)
int
f_integrand(acb_ptr result, const acb_t xc, void *param, slong order, slong prec)
{
	arb_t ar, br, cr, x;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	arb_init(x);
	
	acb_t ac;
	acb_init(ac);

	arb_ptr(a);
	arb_ptr(b);
	arb_ptr(n);
	
	slong sw;

	a = ((arb_ptr *)(param))[0];
	b = ((arb_ptr *)(param))[1];
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

	hyp(br, x, a, b, n, prec);
	arb_mul(ar, ar, br, prec);
	
//	//(1/(x+1))^(a+b+2)
//	arb_add_si(br, x, 1, prec);
//	arb_add(cr, a, b, prec);
//	arb_add_si(cr, cr, 2, prec);
//	arb_neg(cr, cr);
//	arb_pow(br, br, cr, prec);
//	arb_mul(ar, ar, br, prec);
//	
//	//x^b
//	arb_pow(br, x, b, prec);
//	arb_mul(ar, ar, br, prec);
	
	acb_set(ac, xc);
	wei(br, ac, a, b, prec);
	arb_mul(ar, ar, br, prec);
	
	acb_set_arb(result, ar);

	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);	
	arb_clear(x);

	acb_clear(ac);

    return 0;
}

// Evaluate the integral and multiply output by the nfactor
void
integral_coeff(arb_t result, arb_t a, arb_t b, arb_t n, slong sw, slong prec)
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
	//arb_set_str(ar, "10e-600", prec);
	arb_set_str(ar, "10e-400", prec);
	arb_one(br);
	arb_sub(ar, br, ar, prec);
	acb_set_arb(bi, ar);
	
	arb_neg(ar, ar);
	acb_set_arb(ai, ar);

	void *param[4];
	param[0] = (void *) a;
	param[1] = (void *) b;
	param[2] = (void *) n;
	param[3] = (void *) sw;

	//evaluate integral
	acb_calc_integrate(ac, f_integrand, param, ai, bi, goal, tol, NULL, prec);
		
	acb_get_real(result, ac);
	
	//res * nfactor(n)
	nfactor(ar, a, b, n, prec);
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
    arb_t a, b, n, res;
    arb_init(a);
    arb_init(b);
    arb_init(n);
    arb_init(res);

	slong i, prec, numthreads, sw;

	int result = EXIT_SUCCESS;
	
	if (argc != 5)
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	//set precision and the integral upper bound
	//prec = 30 * 3.32192809488736 + 3800;
	prec = 30 * 3.32192809488736 + 2800;
	
	arb_set_str(a, argv[1], prec);
	arb_set_str(b, argv[2], prec);
    numthreads = atol(argv[3]);
    sw = atol(argv[4]);
	
	if ((sw < 0) || (sw > 1) || (numthreads < 1))
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	flint_set_num_threads(numthreads);
	
	//loop through n = 0 to 300 and compute each coefficient
	for (i = 0; i <= 260; i++)
    {
		arb_set_si(n, i);
		integral_coeff(res, a, b, n, sw, prec);
		arb_printn(res, 120, ARB_STR_NO_RADIUS);
		printf("\n");
    }

finish:

    if (result == EXIT_FAILURE)
    {
		flint_printf("\n");
		flint_printf("Required inputs:\n");
		flint_printf("%s a b numthreads sw \n\n", argv[0]);
		flint_printf(
    "This script computes the first 300 coefficients for the expansion of the Riemann Xi-function\n"
    "into the orthogonal Jacobi polynomials as described in the paper:\n"
    ">> Visualising the flows of orthogonal polynomial expansions of the Riemann Xi-function <<\n"
    "-> 'numthreads' defines the number of desired threads/available cores to increase speed.\n"
    "-> 'a' and 'b' are the required parameters.\n"
    "-> 'sw' allows to select the function to expand: Xi(1/2+xi) (sw=0) or 2/x*sin(x/2) (sw=1)\n"
    "The user can push the output into a file by adding: >yyycoeffparms.txt \n");
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(n);
    arb_clear(res);
 
    flint_cleanup();

    return 0;  
} 


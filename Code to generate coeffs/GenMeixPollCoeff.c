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
	arb_set_str(ar, "10e-200", prec);
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
wei(acb_t result, arb_t x, arb_t a, slong prec)
{  
	acb_t ac, bc;
	acb_init(ac);
	acb_init(bc);
	
	//gammma(a+ix)*gammma(a-ix)
	acb_set_arb_arb(ac, a, x);
	acb_conj(bc, ac);
	acb_gamma(ac, ac, prec);
	acb_gamma(bc, bc, prec);

	acb_mul(result, ac, bc, prec);
	
	acb_clear(ac);
	acb_clear(bc);
}

// Compute the hypergeometric function
void
hyp(acb_t result, arb_t x, arb_t a, arb_t n, slong prec)
{  
	arb_t ar;
	arb_init(ar);
	
	acb_struct p[2];
	acb_init(p);
	acb_init(p + 1);

	acb_struct q[1];
	acb_init(q);
	
	acb_t z;
	acb_init(z);
	
	//fill the p values
	arb_neg(ar, n);
	acb_set_arb(p + 0, ar);
	acb_set_arb_arb(p + 1, a, x);
	
	//fill the q values
	arb_add(ar, a, a, prec);
	acb_set_arb(q + 0, ar);
	
	//fill z value
	acb_set_si(z, 2);
	
	//Kummer3F2([-n,a+i*x, a-I*x],[a+b, a+c],1)
	acb_hypgeom_pfq(result, p, 2, q, 1, z, 0, prec);
	
	arb_clear(ar);
	
	acb_clear(p);
	acb_clear(p + 1);
	acb_clear(q);
	acb_clear(z);
}

// Compute the nfactor
void
nfactor(arb_t result, arb_t a, arb_t n, slong prec)
{  
	arb_t aa, ar, br, cr, pi;
	arb_init(aa);
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	arb_init(pi);

	arb_const_pi(pi, prec);
	
	arb_add(aa, a, a, prec);
	
	//gammma(n+a+a)*4^a
	arb_set_si(ar, 4);
	arb_pow(ar, ar, a, prec);
	
	arb_add(br, n, aa, prec);
	arb_gamma(br, br, prec);
	arb_mul(ar, ar, br, prec);
	
	//2*Pi*gamma(a+a)^2*gamma(n+1)
	arb_gamma(cr, aa, prec);
	arb_pow_ui(cr, cr, 2, prec);
	arb_mul(br, pi, cr, prec);
	arb_mul_2exp_si(br, br, 1);
	
	arb_add_si(cr, n, 1, prec);
	arb_gamma(cr, cr, prec);
	arb_mul(br, br, cr, prec);
	
	//divide both of the above
	arb_div(result, ar, br, prec);
	
	arb_clear(aa);
	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
	arb_clear(pi);
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
	
	slong sw;

	a = ((arb_ptr *)(param))[0];
	n = ((arb_ptr *)(param))[1];
	sw = ((slong *)(param))[2];
	
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
	
	hyp(ac, x, a, n, prec);
	acb_mul_arb(ac, ac, ar, prec);
	
	wei(bc, x, a, prec);
	acb_mul(result, ac, bc, prec);

	arb_clear(ar);	
	arb_clear(br);	
	arb_clear(x);

	acb_clear(ac);
	acb_clear(bc);

    return 0;
}

// Evaluate the integral and multiply output by the nfactor
void
integral_coeff(arb_t result, arb_t a, arb_t n, slong intlim, slong sw, slong prec)
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
	arb_zero(ar);
	arb_set_si(br, intlim);
	acb_set_arb(ai, ar);
	acb_set_arb(bi, br);

	void *param[3];
	param[0] = (void *) a;
	param[1] = (void *) n;
	param[2] = (void *) sw;

	//evaluate integral
	acb_calc_integrate(ac, f_integrand, param, ai, bi, goal, tol, NULL, prec);
		
	acb_get_real(result, ac);
	
	//res * nfactor(n)
	nfactor(ar, a, n, prec);
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
	
	if (argc != 4)
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	//set precision and the integral upper bound
	prec = 30 * 3.32192809488736 + 2400;
	//intlim = 390;
	intlim = 300;
	
	arb_set_str(a, argv[1], prec);
    numthreads = atol(argv[2]);
    sw = atol(argv[3]);
	
	if ((sw < 0) || (sw > 1) || (numthreads < 1) || (arb_is_negative(a)))
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	flint_set_num_threads(numthreads);
	
	//loop through n = 0 to 200 and compute each coefficient
	for (i = 351; i <= 400; i++)
    {
		arb_set_si(n, i);
		//no need to compute the coefficient for n=odd
		if (i % 2)
			arb_zero(res);
		else 
			integral_coeff(res, a, n, intlim, sw, prec);
		arb_printn(res, 120, ARB_STR_NO_RADIUS);
		printf("\n");
    }

finish:

    if (result == EXIT_FAILURE)
    {
		flint_printf("\n");
		flint_printf("Required inputs:\n");
		flint_printf("%s a numthreads sw \n\n", argv[0]);
		flint_printf(
    "This script computes the first 200 coefficients for the expansion of the Riemann Xi-function\n"
    "into the orthogonal symmetric Meixner Pollaczek polynomials (phi = pi/2) as described in the paper:\n"
    ">> Visualising the flows of orthogonal polynomial expansions of the Riemann Xi-function <<\n"
    "-> 'numthreads' defines the number of desired threads/available cores to increase speed.\n"
    "-> 'a' is the required free parameter that needs be real and larger than 0\n"
    "-> 'sw' allows to select the function to expand: Xi(1/2+xi) (sw=0) or 2/x*sin(x/2) (sw=1)\n"
    "The user can push the output into a file by adding: >yyycoeffparms.txt \n");
    }

    arb_clear(a);
    arb_clear(n);
    arb_clear(res);
 
    flint_cleanup();

    return 0;  
} 

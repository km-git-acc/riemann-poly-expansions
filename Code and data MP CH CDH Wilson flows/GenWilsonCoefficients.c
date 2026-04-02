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
Poch(arb_t result, arb_t a1, arb_t b1, slong prec)
{  
	arb_t ar, br;
	arb_init(ar);
	arb_init(br);
	
	arb_add(ar, a1, b1, prec);
	arb_gamma(ar, ar, prec);
	
	arb_gamma(br, a1, prec);
	
	arb_div(result, ar, br, prec);

	arb_clear(ar);
	arb_clear(br);
}

void
OO(arb_t sum, const arb_t x, slong N, slong prec)

{   
    arb_t ar, br, exp52x, exp92x, exp2x, pi, pi2;
	arb_init(ar);
	arb_init(br);
	arb_init(exp52x);
	arb_init(exp92x);
	arb_init(exp2x);
	arb_init(pi);
	arb_init(pi2);
	arb_const_pi(pi, prec);
	
	slong n, np2, np4;
	
	//Pi^2
	arb_pow_ui(pi2, pi, 2, prec);
	
	//exp(5/2*x)
	arb_mul_si(ar, x, 5, prec);
	arb_mul_2exp_si(ar, ar, -1);
	arb_exp(exp52x, ar, prec);
	
	//exp(9/2*x)
	arb_mul_si(ar, x, 9, prec);
	arb_mul_2exp_si(ar, ar, -1);
	arb_exp(exp92x, ar, prec);
	
	//exp(2*x)
	arb_mul_2exp_si(ar, x, 1);
	arb_exp(exp2x, ar, prec);
	
	arb_zero(sum);
	for (n = 1; n <= N; n++)
	{
		np2 = pow(n, 2);
		np4 = pow(n, 4);
		
		//4*n^4*Pi^2*exp(9/2*x)
		arb_mul(ar, exp92x, pi2, prec);
		arb_mul_si(ar, ar, np4, prec);
		arb_mul_si(ar, ar, 4, prec);
		
		//6*n^2*Pi*exp(5/2*x)
		arb_mul(br, exp52x, pi, prec);
		arb_mul_si(br, br, np2, prec);
		arb_mul_si(br, br, 6, prec);
		
		arb_sub(ar, ar, br, prec);
		
		//exp(-n^2*Pi*exp(2*x))
		arb_mul(br, exp2x, pi, prec);
		arb_mul_si(br, br, np2, prec);
		arb_neg(br, br);
		arb_exp(br, br, prec);
		
		arb_mul(ar, ar, br, prec);
		
		arb_add(sum, sum, ar, prec);
	}

	arb_clear(ar);
	arb_clear(br);
	arb_clear(exp52x);
	arb_clear(exp92x);
	arb_clear(exp2x);
	arb_clear(pi);
	arb_clear(pi2);
}

// Compute the hypergeometric function
void
Hyp4f3(acb_t result, arb_t a, arb_t b, arb_t n, arb_t k, slong prec)
{  
	arb_t ar, half;
	arb_init(ar);
	arb_init(half);
	
	acb_struct p[4];
	acb_init(p);
	acb_init(p + 1);
	acb_init(p + 2);
	acb_init(p + 3);

	acb_struct q[3];
	acb_init(q);
	acb_init(q + 1);
	acb_init(q + 2);
	
	acb_t z;
	acb_init(z);
	
	arb_one(ar);
	arb_mul_2exp_si(half, ar, -1);
	
	//fill the p values
	arb_neg(ar, k);
	acb_set_arb(p + 0, ar);
	
	arb_add(ar, n, n, prec);
	arb_add(ar, ar, k, prec);	
	arb_add(ar, ar, a, prec);	
	arb_add(ar, ar, b, prec);	
	arb_sub(ar, ar, half, prec);	
	acb_set_arb(p + 1, ar); 
	
	arb_add(ar, n, a, prec);
	arb_add(ar, ar, b, prec);	
	acb_set_arb(p + 2, ar);
	
	arb_add(ar, n, a, prec);
	arb_add(ar, ar, a, prec);	
	acb_set_arb(p + 3, ar);
	
	
	//fill the q values
	arb_add(ar, n, a, prec);
	arb_add(ar, ar, b, prec);
	arb_mul_2exp_si(ar, ar, 1);
	acb_set_arb(q + 0, ar);
	
	arb_add(ar, n, a, prec);
	arb_add(ar, ar, half, prec);
	acb_set_arb(q + 1, ar);
	
	arb_add(ar, n, a, prec);
	acb_set_arb(q + 2, ar);
	
	//fill z value
	acb_set_si(z, 1);
	
	//Kummer3F2([-n,a+i*x, a-I*x],[a+b, a+c],1)
	acb_hypgeom_pfq(result, p, 4, q, 3, z, 0, prec);
	
	arb_clear(ar);
	arb_clear(half);
	
	acb_clear(p);
	acb_clear(p + 1);
	acb_clear(p + 2);
	acb_clear(p + 3);
	acb_clear(q);
	acb_clear(q + 1);
	acb_clear(q + 2);
	acb_clear(z);
}

// Compute the hypergeometric function
void
Hyp2f1(acb_t result, acb_t x, arb_t n, arb_t k, arb_t a, arb_t b, slong prec)
{  
	arb_t ar, half;
	arb_init(ar);
	arb_init(half);
	
	acb_struct p[2];
	acb_init(p);
	acb_init(p + 1);

	acb_struct q[1];
	acb_init(q);
	
	acb_t z;
	acb_init(z);
	
	arb_one(ar);
	arb_mul_2exp_si(half, ar, -1);
	
	//fill the p values
	arb_add(ar, a, n, prec);
	arb_add(ar, ar, k, prec);	
	acb_set_arb(p + 0, ar);
	
	arb_add(ar, ar, half, prec);
	acb_set_arb(p + 1, ar);
	
	
	//fill the q values
	arb_add(ar, n, n, prec);
	arb_add(ar, ar, k, prec);
	arb_add(ar, ar, k, prec);
	arb_add(ar, ar, a, prec);
	arb_add(ar, ar, b, prec);
	arb_add(ar, ar, half, prec);
	acb_set_arb(q + 0, ar);
	
	//fill z value
	acb_mul_2exp_si(z, x, -1);
	acb_tanh(z, z, prec);
	acb_pow_ui(z, z, 2, prec);
	
	//Kummer3F2([-n,a+i*x, a-I*x],[a+b, a+c],1)
	acb_hypgeom_pfq(result, p, 2, q, 1, z, 0, prec);
	
	arb_clear(ar);
	arb_clear(half);
	
	acb_clear(p);
	acb_clear(p + 1);
	acb_clear(q);
	acb_clear(z);
}

void
Pochfact(arb_t res, arb_t a, arb_t b, arb_t n, arb_t k, slong prec)
{   
    arb_t ar, br, cr, dr, half;
	arb_init(ar);
	arb_init(br);
	arb_init(cr);
	arb_init(dr);
	arb_init(half);
	
	arb_one(ar);
	arb_mul_2exp_si(half, ar, -1);
	
	arb_add(ar, a, n, prec);
	Poch(br, ar, k, prec);
	
	arb_add(ar, ar, half, prec);
	Poch(ar, ar, k, prec);
	
	arb_mul(ar, ar, br, prec);
	
	arb_add(br, a, b, prec);
	arb_add(br, br, n, prec);
	Poch(br, br, k, prec);
	
	arb_mul(ar, ar, br, prec);
	
	arb_add(br, a, b, prec);
	arb_add(br, br, n, prec);
	arb_add(br, br, n, prec);
	arb_sub(br, br, half, prec);
	Poch(br, br, k, prec);
	
	arb_mul(ar, ar, br, prec);
	
	arb_add(br, n, half, prec);
	Poch(br, br, k, prec);
	
	arb_add(cr, a, b, prec);
	arb_add(cr, cr, n, prec);
	arb_add(cr, cr, n, prec);
	arb_sub(cr, cr, half, prec);
	arb_add(dr, k, k, prec);
	Poch(cr, cr, dr, prec);
	
	arb_mul(br, br, cr, prec);
	
	arb_add_si(cr, k, 1, prec);
	arb_gamma(cr, cr, prec);
	
	arb_mul(br, br, cr, prec);
	
	arb_div(res, ar, br, prec);

	arb_clear(ar);
	arb_clear(br);
	arb_clear(cr);
	arb_clear(dr);
	arb_clear(half);
}

// Prepare the integrand OO(x) * hyp2F1(x) * tanh(x/2)^(2*k) * tanh(x/2)^(2*n)
int
f_integrand(acb_ptr result, const acb_t xc, void *param, slong order, slong prec)
{
	arb_t ar, x;
	arb_init(ar);
	arb_init(x);
	
	acb_t ac, bc;
	acb_init(ac);
	acb_init(bc);

	arb_ptr(a);
	arb_ptr(b);
	arb_ptr(n);
	arb_ptr(k);
	
	slong sw; 

	a = ((arb_ptr *)(param))[0];
	b = ((arb_ptr *)(param))[1];
	n = ((arb_ptr *)(param))[2];
	k = ((arb_ptr *)(param))[3];
	sw = ((slong *)(param))[4];
	
	acb_get_real(x, xc);
	
	if (sw == 0)
		OO(ar, x, 12, prec);
	else
		arb_one(ar);
	
	acb_set_arb(ac, x);
	Hyp2f1(ac, ac, n, k, a, b, prec);
	
	acb_mul_arb(ac, ac, ar, prec);
	
	acb_mul_2exp_si(bc, xc, -1);
	acb_tanh(bc, bc, prec);
	arb_add(ar, n, k, prec);
	arb_mul_2exp_si(ar, ar, 1);
	acb_pow_arb(bc, bc, ar, prec);

	acb_mul(ac, ac, bc, prec);
	
	acb_mul_2exp_si(bc, xc, -1);
	acb_sech(bc, bc, prec);
	arb_add(ar, a, a, prec);
	acb_pow_arb(bc, bc, ar, prec);
	
	acb_mul(result, ac, bc, prec);

	arb_clear(ar);	
	arb_clear(x);

	acb_clear(ac);
	acb_clear(bc);

    return 0;
}

// Evaluate the integral and multiply output by the nfactor
void
integral_coeff(arb_t result, arb_t a, arb_t b, arb_t n, arb_t k, slong intlim, slong sw, slong prec)
{  
	arb_t ar, br;
	arb_init(ar);
	arb_init(br);
	
	acb_t ac, ai, bi;
	acb_init(ac);
	acb_init(ai);
	acb_init(bi);

	slong goal;

    //evaluate the integrals
	acb_calc_integrate_opt_t options;
	acb_calc_integrate_opt_init(options);

	mag_t tol;
	mag_init(tol);

	goal = prec;
	mag_set_ui_2exp_si(tol, 1, -prec);

	//set limits of integration, note that the lower limit must be 0+
	arb_set_str(ar,"10e-80",prec);
	if (sw == 0)
		arb_set_si(br, intlim);
	else
		arb_set_str(br,"0.5",prec);
	acb_set_arb(ai, ar);
	acb_set_arb(bi, br);

	void *param[5];
	param[0] = (void *) a;
	param[1] = (void *) b;
	param[2] = (void *) n;
	param[3] = (void *) k;
	param[4] = (void *) sw;

	//evaluate integral
	acb_calc_integrate(ac, f_integrand, param, ai, bi, goal, tol, NULL, prec);
		
	acb_get_real(result, ac);
	
	//pochfactor
	Pochfact(ar, a, b, n, k, prec);
	
	arb_mul(result, result, ar, prec);
	
	Hyp4f3(ac, a, b, n, k, prec);
	acb_get_real(ar, ac);
	
	arb_mul(result, result, ar, prec);
	
	mag_clear(tol);
	
	arb_clear(ar);
	arb_clear(br);

	acb_clear(ac);
	acb_clear(ai);
	acb_clear(bi);
}

int main(int argc, char *argv[])
{
    arb_t a, b, ar, br, k, ksum, n, x, res;
    arb_init(a);
    arb_init(b);
    arb_init(ar);
    arb_init(br);
    arb_init(k);
    arb_init(ksum);
    arb_init(n);
    arb_init(x);
    arb_init(res);
	
	acb_t ac;
	acb_init(ac);

	slong i, ki, prec, numthreads, intlim, sw;

	int result = EXIT_SUCCESS;
	
	if (argc != 5)
    {
        result = EXIT_FAILURE;
        goto finish;
    }
	
	//set precision and the integral upper bound
	prec = 30 * 3.32192809488736 + 8402;
	intlim = 4;
	
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

//	arb_set_str(a,"1.1",prec);
//	arb_set_str(b,"2.1",prec);
//	arb_set_str(n,"3",prec);
//	arb_set_str(k,"1",prec);
//	arb_set_str(x,"3.1",prec);
//	Pochfact(ar, a, b, n, k, prec);
//	arb_printn(ar, 120, ARB_STR_NO_RADIUS);
//	printf("\n");
//	Hyp4f3(ac, a, b, n, k, prec);
//	acb_printn(ac, 120, ARB_STR_NO_RADIUS);
//	printf("\n");
	acb_clear(ac);

	//loop through n = 0 to 200 and compute each coefficient
	for (i = 1499; i <= 1501; i++)
    {
		arb_set_si(n, i);
		arb_zero(ksum);
		for (ki = 0; ki <= 75; ki++)
		{
			arb_set_si(k, ki);
			integral_coeff(res, a, b, n, k, intlim, sw, prec);
			arb_add(ksum, ksum, res, prec);
		}
	
		arb_one(ar);
		arb_mul_2exp_si(ar, ar, 2 * i + 1 );
		arb_mul(ksum, ksum, ar, prec);
		
		arb_add(ar, n, n, prec);
		arb_add_si(ar, ar, 1, prec);
		arb_gamma(ar, ar, prec);
		arb_div(ksum, ksum, ar, prec);
		
		arb_add(ar, a, a, prec);
		arb_add(ar, ar, b, prec);
		arb_add(ar, ar, b, prec);
		arb_add(ar, ar, n, prec);
		arb_sub_si(ar, ar, 1, prec);
		Poch(ar, ar, n, prec);
		
		arb_div(ksum, ksum, ar, prec);
		
		arb_printn(ksum, 120, ARB_STR_NO_RADIUS);
		printf("\n");
		
	}

finish:

    if (result == EXIT_FAILURE)
    {
		flint_printf("\n");
		flint_printf("Required inputs:\n");
		flint_printf("%s a b numthreads sw \n\n", argv[0]);
		flint_printf(
    "This script computes the first 150 coefficients for the expansion of the Riemann Xi-function\n"
    "into the orthogonal Jacobi polynomials as described in the paper:\n"
    ">> Visualising the flows of orthogonal polynomial expansions of the Riemann Xi-function <<\n"
    "-> 'numthreads' defines the number of desired threads/available cores to increase speed.\n"
    "-> 'a, b' are the required parameters. Note that a + b  > 1 and all real.\n"
    "-> 'sw' allows to select the function to expand: Xi(1/2+xi) (sw=0) or 2/x*sin(x/2) (sw=1)\n"
    "The user can push the output into a file by adding: >yyycoeffparms.txt \n");
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(ar);
    arb_clear(br);
    arb_clear(k);
    arb_clear(ksum);
    arb_clear(n);
    arb_clear(x);
    arb_clear(res);
 
    flint_cleanup();

    return 0;  
} 
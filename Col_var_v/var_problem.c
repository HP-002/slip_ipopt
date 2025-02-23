#include "ipInterface.h"

int main()
{
    ipindex n = -1;                      /* number of variables */
    ipindex m = -1;                      /* number of constraints */
    ipindex nele_jac;                    /* number of nonzeros in the Jacobian of the constraints */
    ipindex nele_hess;                   /* number of nonzeros in the Hessian of the Lagrangian (lower or upper triangular part only) */
    ipindex index_style;                 /* indexing style for matrices */
    ipnumber *x_L = NULL;                /* lower bounds on x */
    ipnumber *x_U = NULL;                /* upper bounds on x */
    ipnumber *g_L = NULL;                /* lower bounds on g */
    ipnumber *g_U = NULL;                /* upper bounds on g */
    IpoptProblem nlp = NULL;             /* IpoptProblem */
    enum ApplicationReturnStatus status; /* Solve return code */
    ipnumber *x = NULL;                  /* starting point and solution vector */
    ipnumber *mult_g = NULL;             /* constraint multipliers at the solution */
    ipnumber *mult_x_L = NULL;           /* lower bound multipliers at the solution */
    ipnumber *mult_x_U = NULL;           /* upper bound multipliers at the solution */
    ipnumber obj;                        /* objective value */
    struct MyUserData user_data;         /* our user data for the function evaluations */
    ipindex i;                           /* generic counter */
    /* set the number of variables and allocate space for the bounds */
    n = 4;
    x_L = (ipnumber *)malloc(sizeof(ipnumber) * n);
    x_U = (ipnumber *)malloc(sizeof(ipnumber) * n);
    /* set the values for the variable bounds */
    for (i = 0; i < n; i++)
    {
        x_L[i] = 1.0;
        x_U[i] = 5.0;
    }

    /* set the number of constraints and allocate space for the bounds */
    m = 2;
    g_L = (ipnumber *)malloc(sizeof(ipnumber) * m);
    g_U = (ipnumber *)malloc(sizeof(ipnumber) * m);
    /* set the values of the constraint bounds */
    g_L[0] = 25;
    g_U[0] = 2e19;
    g_L[1] = 40;
    g_U[1] = 40;

    /* set the number of nonzeros in the Jacobian and Hessian */
    nele_jac = 8;
    nele_hess = 10;

    /* set the indexing style to C-style (start counting of rows and column indices at 0) */
    index_style = 0;

    /* create the IpoptProblem */
    nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, index_style,
                             &eval_f, &eval_g, &eval_grad_f,
                             &eval_jac_g, &eval_h);

    /* We can free the memory now - the values for the bounds have been
     * copied internally in CreateIpoptProblem
     */
    free(x_L);
    free(x_U);
    free(g_L);
    free(g_U);

    /* Set some options.  Note the following ones are only examples,
     * they might not be suitable for your problem.
     */
    AddIpoptNumOption(nlp, "tol", 3.82e-6);
    AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
    AddIpoptStrOption(nlp, "output_file", "ipopt.out");

    /* allocate space for the initial point and set the values */
    x = (ipnumber *)malloc(sizeof(ipnumber) * n);
    x[0] = 1.0;
    x[1] = 5.0;
    x[2] = 5.0;
    x[3] = 1.0;

    /* allocate space to store the bound multipliers at the solution */
    mult_g = (ipnumber *)malloc(sizeof(ipnumber) * m);
    mult_x_L = (ipnumber *)malloc(sizeof(ipnumber) * n);
    mult_x_U = (ipnumber *)malloc(sizeof(ipnumber) * n);

    /* Initialize the user data */
    user_data.g_offset[0] = 0.;
    user_data.g_offset[1] = 0.;
    user_data.nlp = nlp;

    /* Set the callback method for intermediate user-control.
     * This is not required, just gives you some intermediate control in
     * case you need it.
     */
    /* SetIntermediateCallback(nlp, intermediate_cb); */

    /* solve the problem */
    status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, &user_data);

    if (status == Solve_Succeeded)
    {
        printf("\n\nSolution of the primal variables, x\n");
        for (i = 0; i < n; i++)
        {
            printf("x[%d] = %e\n", (int)i, x[i]);
        }

        printf("\n\nSolution of the constraint multipliers, lambda\n");
        for (i = 0; i < m; i++)
        {
            printf("lambda[%d] = %e\n", (int)i, mult_g[i]);
        }
        printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
        for (i = 0; i < n; i++)
        {
            printf("z_L[%d] = %e\n", (int)i, mult_x_L[i]);
        }
        for (i = 0; i < n; i++)
        {
            printf("z_U[%d] = %e\n", (int)i, mult_x_U[i]);
        }

        printf("\n\nObjective value\nf(x*) = %e\n", obj);
    }
    else
    {
        printf("\n\nERROR OCCURRED DURING IPOPT OPTIMIZATION.\n");
    }

    /* Now we are going to solve this problem again, but with slightly
     * modified constraints.  We change the constraint offset of the
     * first constraint a bit, and resolve the problem using the warm
     * start option.
     */
    user_data.g_offset[0] = 0.2;

    if (status == Solve_Succeeded)
    {
        /* Now resolve with a warmstart. */
        AddIpoptStrOption(nlp, "warm_start_init_point", "yes");
        /* The following option reduce the automatic modification of the
         * starting point done my Ipopt.
         */
        AddIpoptNumOption(nlp, "bound_push", 1e-5);
        AddIpoptNumOption(nlp, "bound_frac", 1e-5);
        status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, &user_data);

        if (status == Solve_Succeeded)
        {
            printf("\n\nSolution of the primal variables, x\n");
            for (i = 0; i < n; i++)
            {
                printf("x[%d] = %e\n", (int)i, x[i]);
            }

            printf("\n\nSolution of the constraint multipliers, lambda\n");
            for (i = 0; i < m; i++)
            {
                printf("lambda[%d] = %e\n", (int)i, mult_g[i]);
            }
            printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
            for (i = 0; i < n; i++)
            {
                printf("z_L[%d] = %e\n", (int)i, mult_x_L[i]);
            }
            for (i = 0; i < n; i++)
            {
                printf("z_U[%d] = %e\n", (int)i, mult_x_U[i]);
            }

            printf("\n\nObjective value\nf(x*) = %e\n", obj);
        }
        else
        {
            printf("\n\nERROR OCCURRED DURING IPOPT OPTIMIZATION WITH WARM START.\n");
        }
    }

    /* free allocated memory */
    FreeIpoptProblem(nlp);
    free(x);
    free(mult_g);
    free(mult_x_L);
    free(mult_x_U);

    return (status == Solve_Succeeded) ? EXIT_SUCCESS : EXIT_FAILURE;
}
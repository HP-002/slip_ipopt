/* Header files */
#include "IpStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

__uint64_t Ra_ = 1;
__uint64_t B_ = 2;
__uint64_t M_ = 3;
__uint64_t G_ = 9.8; // Gravity?
__uint64_t L0_ = 5;
__uint64_t K0_ = 6;
__uint64_t C_ = 7;
__uint8_t R_ = 8;
__uint64_t Kt_ = 9;
__uint64_t J_ = 10;
__uint64_t La_ = 11;
__uint64_t Kb_ = 12;

/* MyUserData contains our Problem */
struct MyUserData
{
    ipnumber g_offset[2]; /**< This is an offset for the constraints.  */
    IpoptProblem nlp;     /**< The problem to be solved. Required in intermediate_cb. */
};

/* Callback Functions */


/* ---------------bool eval_h---------------
 * Method to request either the sparsity structure or the values of the Hessian of the Lagrangian.
 *
 * The Jacobian is the matrix of derivatives where the derivative of constraint function g_i with respect
 * to variable x_j is placed in row i and column j.
 */
static bool eval_h(
    ipindex n,
    ipnumber *x,
    bool new_x,
    ipnumber obj_factor,
    ipindex m,
    ipnumber *lambda,
    bool new_lambda,
    ipindex nele_hess,
    ipindex *iRow,
    ipindex *jCol,
    ipnumber *values,
    UserDataPtr user_data)
{
    (void)n;
    (void)new_x;
    (void)m;
    (void)new_lambda;
    (void)user_data;

    if (values == NULL)
    {
        ipindex idx; /* nonzero element counter */
        ipindex row; /* row counter for loop */
        ipindex col; /* col counter for loop */

        /* return the structure. This is a symmetric matrix, fill the lower left
         * triangle only. */

        /* the hessian for this problem is actually dense */
        idx = 0;
        for (row = 0; row < 4; row++)
        {
            for (col = 0; col <= row; col++)
            {
                iRow[idx] = row;
                jCol[idx] = col;
                idx++;
            }
        }

        assert(idx == nele_hess);
        (void)nele_hess;
    }
    else
    {
        /* return the values. This is a symmetric matrix, fill the lower left
         * triangle only */

        /* fill the objective portion */
        values[0] = obj_factor * (2 * x[3]); /* 0,0 */

        values[1] = obj_factor * (x[3]); /* 1,0 */
        values[2] = 0;                   /* 1,1 */

        values[3] = obj_factor * (x[3]); /* 2,0 */
        values[4] = 0;                   /* 2,1 */
        values[5] = 0;                   /* 2,2 */

        values[6] = obj_factor * (2 * x[0] + x[1] + x[2]); /* 3,0 */
        values[7] = obj_factor * (x[0]);                   /* 3,1 */
        values[8] = obj_factor * (x[0]);                   /* 3,2 */
        values[9] = 0;                                     /* 3,3 */

        /* add the portion for the first constraint */
        values[1] += lambda[0] * (x[2] * x[3]); /* 1,0 */

        values[3] += lambda[0] * (x[1] * x[3]); /* 2,0 */
        values[4] += lambda[0] * (x[0] * x[3]); /* 2,1 */

        values[6] += lambda[0] * (x[1] * x[2]); /* 3,0 */
        values[7] += lambda[0] * (x[0] * x[2]); /* 3,1 */
        values[8] += lambda[0] * (x[0] * x[1]); /* 3,2 */

        /* add the portion for the second constraint */
        values[0] += lambda[1] * 2; /* 0,0 */

        values[2] += lambda[1] * 2; /* 1,1 */

        values[5] += lambda[1] * 2; /* 2,2 */

        values[9] += lambda[1] * 2; /* 3,3 */
    }

    return true;
} // End of eval_h
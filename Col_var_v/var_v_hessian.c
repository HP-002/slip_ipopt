#include "var_header.c";

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

        /* There are 55 variables, thus Hessian will be a 55x55 Matrix. 
         * Hence, there we need to fill in values for (55*56)/2 = 1540 values */

        /* the hessian for this problem is actually dense */
        idx = 0;
        for (row = 0; row < 55; row++)
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
        /* There are 55 variables, thus Hessian will be a 55x55 Matrix. 
         * Hence, there we need to fill in values for (55*56)/2 = 1540 values */

        /* Objective Hessian */
        for (int i = 0; i < 1540; i++)
        {
            values[i] = 0;
        }

        values[1] = obj_factor * x[1]/(8.*Ra_); // (1,0)
        values[2] = obj_factor * x[0]/(8.*Ra_); // (1,1)
        values[3] = obj_factor * x[2]/(4.*Ra_); // (2,0)
        values[5] = obj_factor * x[0]/(4.*Ra_); // (2,2)
        values[6] = obj_factor * x[3]/(4.*Ra_); // (3,0)
        values[9] = obj_factor * x[0]/(4.*Ra_); // (3,3)
        values[10] = obj_factor * x[4]/(4.*Ra_); // (4,0)
        values[14] = obj_factor * x[0]/(4.*Ra_); // (4,4)
        values[15] = obj_factor * x[5]/(4.*Ra_); // (5,0)
        values[20] = obj_factor * x[0]/(4.*Ra_); // (5,5)
        values[21] = obj_factor * x[6]/(4.*Ra_); // (6,0)
        values[27] = obj_factor * x[0]/(4.*Ra_); // (6,6)
        values[28] = obj_factor * x[7]/(4.*Ra_); // (7,0)
        values[35] = obj_factor * x[0]/(4.*Ra_); // (7,7)
        values[36] = obj_factor * x[8]/(4.*Ra_); // (8,0)
        values[44] = obj_factor * x[0]/(4.*Ra_); // (8,8)
        values[45] = obj_factor * x[9]/(8.*Ra_); // (9,0)
        values[54] = obj_factor * x[0]/(8.*Ra_); // (9,9)

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
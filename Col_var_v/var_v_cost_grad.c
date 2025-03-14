#include "var_header.c";

/* ---------------bool eval_grad_f---------------
 * Method to request the gradient of the objective function.
 */
static bool eval_grad_f(
    ipindex n,
    ipnumber *x,
    bool new_x,
    ipnumber *grad_f,
    UserDataPtr user_data)
{
    assert(n == 55);
    (void)n;

    (void)new_x;
    (void)user_data;

    grad_f[0] = (pow(x[9], 2) + pow(x[1], 2) + (2 * pow(x[2], 2)) + (2 * pow(x[3], 2)) + (2 * pow(x[4], 2)) + (2 * pow(x[5], 2)) + (2 * pow(x[6], 2)) + (2 * pow(x[7], 2)) + (2 * pow(x[8], 2))) / (16 * Ra_); // Derivatives of functions w.r.t all variables
    grad_f[1] = x[0] * x[1] / (8 * Ra_);
    grad_f[2] = x[0] * x[2] / (4 * Ra_);
    grad_f[3] = x[0] * x[3] / (4 * Ra_);
    grad_f[4] = x[0] * x[4] / (4 * Ra_);
    grad_f[5] = x[0] * x[5] / (4 * Ra_);
    grad_f[6] = x[0] * x[6] / (4 * Ra_);
    grad_f[7] = x[0] * x[7] / (4 * Ra_);
    grad_f[8] = x[0] * x[8] / (4 * Ra_);
    grad_f[9] = x[0] * x[9] / (8 * Ra_);

    for (int i = 10; i < n; i++)
    {
        grad_f[i] = 0;
    }

    return true;
} // End of eval_grad_f

/* Header files */
#include "IpStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

__uint64_t Ra_ = 1;
__uint64_t b_ = 2;
__uint64_t m_ = 3;
__uint64_t g_ = 4; // Gravity?
__uint64_t l0_ = 5;
__uint64_t k0_ = 6;
__uint64_t c_ = 7;
__uint8_t R_ = 8;
__uint64_t kt_ = 9;
__uint64_t J_ = 10;
__uint64_t La_ = 11;
__uint64_t kb_ = 12;

/* MyUserData contains our Problem */
struct MyUserData
{
    ipnumber g_offset[2]; /**< This is an offset for the constraints.  */
    IpoptProblem nlp;     /**< The problem to be solved. Required in intermediate_cb. */
};

/* Callback Functions */

/* ---------------bool eval_f---------------
 * Method to request the value of the objective function.
 */
static bool eval_f(
    ipindex n,   // Number of variables
    ipnumber *x, // Variables array
    bool new_x,
    ipnumber *obj_value, // Function
    UserDataPtr user_data)
{
    assert(n == 55);
    (void)n;

    (void)new_x;
    (void)user_data;
    *obj_value = x[0] *
                 (pow(x[9], 2) + pow(x[1], 2) + (2 * pow(x[2], 2)) + (2 * pow(x[3], 2)) + (2 * pow(x[4], 2)) +
                  (2 * pow(x[5], 2)) + (2 * pow(x[6], 2)) + (2 * pow(x[7], 2)) + (2 * pow(x[8], 2))) /
                 (16 * Ra_); // Function to be minimized

    return true;
} // End of eval_f

/* ---------------bool eval_g---------------
 * Method to request the constraint values.
 */
static bool eval_g(
    ipindex n, //
    ipnumber *x,
    bool new_x,
    ipindex m,
    ipnumber *g,
    UserDataPtr user_data)
{
    struct MyUserData *my_data = user_data;

    assert(n == 55);
    (void)n;
    assert(m == 40);
    (void)m;

    (void)new_x;

    /* Equality Constraints */
    for (int i = 0, j = 0; i < 8; i++, j += 5)
    {
        double A1 = 0 - (b_ * x[19 + i] / m) + (b_ * x[20 + i] / m_) + (x[10 + i] * pow(x[37 + i], 2)) - (x[11 + i] * pow(x[38 + i], 2)) - (g_ * sin(x[28 + i])) + (g_ * sin(x[29 + i])) - ((x[10 + i] - l0_) * k0_ / m_) + ((x[11 + i] - l0_) * k0_ / m_);
        double A2 = ((x[19 + i] + x[20 + i]) / 2) + ((x[0] / 64) * (A1));

        double B1 = ((0 - ((c_ * pow(R_, 2) * x[37 + i]) / (m * pow(x[10 + i], 2))) - ((2 * x[19 + i] * x[37 + i]) / (x[10 + i])) + ((kt_ * R_ * x[46 + i]) / (m_ * pow(x[10 + i], 2))) - ((g_ * cos(x[28 + i])) / (x[10 + i]))) / (1 + (pow(J_ * R_, 2) / (m_ * x[10 + i]))));
        double B2 = ((0 - ((c_ * pow(R_, 2) * x[38 + i]) / (m * pow(x[11 + i], 2))) - ((2 * x[20 + i] * x[38 + i]) / (x[11 + i])) + ((kt_ * R_ * x[47 + i]) / (m_ * pow(x[11 + i], 2))) - ((g_ * cos(x[29 + i])) / (x[11 + i]))) / (1 + (pow(J_ * R_, 2) / (m_ * x[11 + i]))));

        double C1 = ((x[10 + i] + x[11 + i]) / 2) + ((x[0] / 64) * (x[19 + i] - x[20 + i]));

        double D1 = ((x[46 + i] + x[47 + i]) / 2) + ((x[0] / 64) * ((x[1 + i] / La_) - (x[2 + i] / La_) - (kb_ * R_ * x[37 + i] / La_) + (kb_ * R_ * x[38 + i] / La_) - (Ra_ * x[46 + i] / La_) + (Ra_ * x[47 + i] / La_)));

        g[j] = x[10 + i] - x[11 + i] + ((x[0] / 48) * (x[19 + i] + x[20 + i] + 4 * A2));

        g[j + 1] = x[19 + i] - x[20 + i] + (x[0] / 48) * (A1 + 4 * ((((x[10 + i] + x[11 + i]) / 2) + ((x[0] / 64) * (x[19 + i] - x[20 + i]))) * pow((((x[37 + i] + x[38 + i]) / 2) + ((x[0] / 64) * ((B1) - (B2)))), 2) - (g_ * sin(((x[28 + i] + x[29 + i]) / 2) + ((x[0] / 64) * (x[37 + i] - x[38 + i])))) - ((0 - l0_ + ((x[10 + i] + x[11 + i]) / 2) + ((x[0] / 64) * (x[19 + i] - x[20 + i]))) * k0_ / m_) - (b_ * (A2) / m_)));

        g[j + 2] = x[28 + i] - x[29 + i] + ((x[0] / 48) * (x[37 + i] + x[38 + i] + (4 * (((x[37 + i] + x[38 + i]) / 2) + ((x[0] / 64) * (B1 - B2))))));

        g[j + 3] = x[37 + i] - x[38 + i] + ((x[0] / 48) * (B1 + B2 + (4 * (((kt_ * R_ * (D1)) / (m * pow(C1, 2))) - ((c_ * pow(R_, 2) * (((x[37 + i] + x[38 + i]) / 2) + ((x[0] / 64) * (B1 - B2)))) / (m * pow(C1, 2))) - ((g_ * cos(((x[28 + i] + x[29 + i]) / 2) + ((x[0] / 64) * (x[37 + i] - x[38 + i])))) / C1) - ((2 * (((x[37 + i] + x[38 + i]) / 2) + ((x[0] / 64) * (B1 - B2))) * A2) / C1)) / (1 + (pow(J_ * R_, 2) / (m_ * C1))))));

        g[j + 4] = x[46 + i] - x[47 + i] + ((x[0] / 48) * ((x[1 + i] / La_) + (x[2 + i] / La_) - (kb_ * R_ * (x[37 + i] + x[38 + i]) / La_) - (R_ * (x[46 + i] + x[47 + i]) / La_) + (4 * (((x[1 + i] + x[2 + i]) / (2 * La_)) - (Ra_ * D1 / La_) - ((kb_ * R_ * (((x[37 + i] + x[38 + i]) / 2) + ((x[0] / 64) * (B1 - B2)))) / La_)))));
    }

    return true;
} // End of eval_g

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

/* ---------------bool eval_jac_g---------------
 * Method to request either the sparsity structure or the values of the Jacobian of the constraints.
 *
 * The Jacobian is the matrix of derivatives where the derivative of constraint function g_i with respect
 * to variable x_j is placed in row i and column j.
 */
static bool eval_jac_g(
    ipindex n,
    ipnumber *x,
    bool new_x,
    ipindex m,
    ipindex nele_jac,
    ipindex *iRow,
    ipindex *jCol,
    ipnumber *values,
    UserDataPtr user_data)
{
    (void)n;
    (void)new_x;
    (void)m;
    (void)nele_jac;
    (void)user_data;

    if (values == NULL)
    {
        int idx = 0;
        /* return the structure of the jacobian */

        //----------- h[i+1]-----------
        for (int i = 0; i < 8; i++) {
            //     -----[[1]]-----
            iRow[idx] = (i * 5);    jCol[idx++] = 0;
            iRow[idx] = (i * 5);    jCol[idx++] = 10 + i;
            iRow[idx] = (i * 5);    jCol[idx++] = 11 + i;
            iRow[idx] = (i * 5);    jCol[idx++] = 19 + i;
            iRow[idx] = (i * 5);    jCol[idx++] = 20 + i;
            iRow[idx] = (i * 5);    jCol[idx++] = 28 + i;
            iRow[idx] = (i * 5);    jCol[idx++] = 29 + i;
            iRow[idx] = (i * 5);    jCol[idx++] = 37 + i;
            iRow[idx] = (i * 5);    jCol[idx++] = 38 + i;

            //     -----[[2]]-----
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 0;
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 10 + i;
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 11 + i;
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 19 + i;
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 20 + i;
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 28 + i;
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 29 + i;
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 37 + i;
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 38 + i;
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 46 + i;
            iRow[idx] = (i * 5) + 1;    jCol[idx++] = 47 + i;

            //     -----[[3]]-----
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 0;
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 10 + i;
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 11 + i;
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 19 + i;
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 20 + i;
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 28 + i;
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 29 + i;
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 37 + i;
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 38 + i;
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 46 + i;
            iRow[idx] = (i * 5) + 2;    jCol[idx++] = 47 + i;

            //     -----[[4]]-----
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 0;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 1 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 2 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 10 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 11 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 19 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 20 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 28 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 29 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 37 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 38 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 46 + i;
            iRow[idx] = (i * 5) + 3;    jCol[idx++] = 47 + i;

            //     -----[[5]]-----
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 0;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 1 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 2 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 10 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 11 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 19 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 20 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 28 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 29 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 37 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 38 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 46 + i;
            iRow[idx] = (i * 5) + 4;    jCol[idx++] = 47 + i;
        }
    }
    else
    {
        /* return the values of the jacobian of the constraints */
        int i;
        int idx = 0;

        for (i = 0; i < 8; i++)
        {
            /*----------- h[i+1]-----------*/
            double A1 = 0 - (b_ * x[19 + i] / m) + (b_ * x[20 + i] / m_) + (x[10 + i] * pow(x[37 + i], 2)) - (x[11 + i] * pow(x[38 + i], 2)) - (g_ * sin(x[28 + i])) + (g_ * sin(x[29 + i])) - ((x[10 + i] - l0_) * k0_ / m_) + ((x[11 + i] - l0_) * k0_ / m_);
            double A2 = ((x[19 + i] + x[20 + i]) / 2) + ((x[0] / 64) * (A1));
            double B1 = ((0 - ((c_ * pow(R_, 2) * x[37 + i]) / (m * pow(x[10 + i], 2))) - ((2 * x[19 + i] * x[37 + i]) / (x[10 + i])) + ((kt_ * R_ * x[46 + i]) / (m_ * pow(x[10 + i], 2))) - ((g_ * cos(x[28 + i])) / (x[10 + i]))) / (1 + (pow(J_ * R_, 2) / (m_ * x[10 + i]))));
            double B2 = ((0 - ((c_ * pow(R_, 2) * x[38 + i]) / (m * pow(x[11 + i], 2))) - ((2 * x[20 + i] * x[38 + i]) / (x[11 + i])) + ((kt_ * R_ * x[47 + i]) / (m_ * pow(x[11 + i], 2))) - ((g_ * cos(x[29 + i])) / (x[11 + i]))) / (1 + (pow(J_ * R_, 2) / (m_ * x[11 + i]))));
            double B3 = ((x[37 + i] + x[39 + i]) / 2) + (x[0] * (B1 - B2) / 64);
            double B4 = ((((2 * c_ * pow(R_, 2) * x[37 + i] / (m_ * pow(x[10 + i], 3))) + (2 * x[19 + i] * x[37 + i] / pow(x[10 + i], 2)) - (2 * kt_ * R_ * x[46 + i] / (m_ * pow(x[10 + i], 3))) + (g_ * cos(x[28 + i]) / pow(x[10 + i], 2))) / (1 + pow(J_ * R_, 2) / (m_ * x[10 + i]))) + (pow(J_ * R_, 2) * B1 / (m_ * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * pow(x[10 + i], 2))));
            double B5 = (0 - ((2 * c_ * pow(R_, 2) * x[38 + i] / (m_ * pow(x[11 + i], 3))) + (2 * x[20 + i] * x[39 + i] / pow(x[11 + i], 2)) - (2 * kt_ * R_ * x[48 + i] / (m_ * pow(x[11 + i], 3)) + (g_ * cos(x[29 + i] / pow(x[11 + i], 2)))) / (1 + pow(J_ * R_, 2) / (m_ * x[11 + i]))) - (pow(J_ * R_, 2) * B2 / (m_ * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * pow(x[11 + i], 2))));
            double B6 = -(c_ * pow(R_, 2) / (m_ * pow(x[10 + i], 2))) - (2 * x[19 + i] / x[10 + i]);
            double B7 = -(c_ * pow(R_, 2) / (m_ * pow(x[11 + i], 2))) - (2 * x[20 + i] / x[12 + i]);
            double C1 = ((x[10 + i] + x[11 + i]) / 2) + ((x[0] / 64) * (x[19 + i] - x[20 + i]));
            double D2 = (x[1 + i] / La_) - (x[2 + i] / La_) - (kb_ * R_ * x[37 + i] / La_) + (kb_ * R_ * x[38 + i] / La_) - (Ra_ * x[46 + i] / La_) + (Ra_ * x[47 + i] / La_);
            double D1 = ((x[46 + i] + x[47 + i]) / 2) + ((x[0] / 64) * (D2));

            /*     -----[[1]]-----     */
            values[idx++] = ((x[0] / 768) * A1) + ((x[19 + i] + x[20 + i] + (4 * A2)) / 48); // (row, 0)
            values[idx++] = 1 + (pow(x[0], 2) / 768) * (pow(x[37 + i], 2) - (k0_ / m));      // (row, 10)
            values[idx++] = -1 + (pow(x[0], 2) / 768) * ((k0_ / m) - pow(x[38 + i], 2));     // (row, 11)
            values[idx++] = (x[0] / 48) * (1 + 4 * (1 / 2 - (b_ * x[0] / (64 * m_))));       // (row, 19)
            values[idx++] = (x[0] / 48) * (1 + 4 * (1 / 2 + (b_ * x[0] / (64 * m_))));       // (row, 20)
            values[idx++] = -1 * g_ * pow(x[0], 2) * cos(x[28 + i]) / 768;                   // (row, 28)
            values[idx++] = g_ * pow(x[0], 2) * cos(x[29 + i]) / 768;                        // (row, 29)
            values[idx++] = pow(x[0], 2) * x[10 + i] * x[37 + i] / 384;                      // (row, 37)
            values[idx++] = -1 * pow(x[0], 2) * x[11 + i] * x[38 + i] / 384;                 // (row, 38)

            /*     -----[[2]]-----     */
            values[idx++] = (x[0] / 12) * ((C1 * (B1 - B2) * B3 / 32) + ((x[19 + i] - x[20 + i]) * pow(B3, 2) / 64) - (g_ * (x[37 + i] - x[38 + i]) * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / 64) - ((x[19 + i] - x[20 + i]) * k0_ / (64 * m_)) - (b_ * A1 / (64 * m_))) + ((A1 + 4 * ((C1 * pow(B3, 2)) - (g_ * sin(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64))) - ((C1 - l0_) * k0_ / m_) - (b_ * A2 / m_))) / 48); // (row+1, 0)
            values[idx++] = x[0] * (pow(x[37 + i], 2) - k0_ / m_ + 4 * ((x[0] * C1 * B4 * B3 / 32) + (pow(B3, 2) / 2) - (k0_ / (2 * m_)) - (b_ * x[0] * (pow(x[37 + i], 2) - k0_ / m_) / (64 * m_)))) / 48;                                                                                                                                                                                                                                                                           // (row+1, 10)
            values[idx++] = x[0] * (pow(x[38 + i], 2) - k0_ / m_ + 4 * ((x[0] * C1 * B5 * B3 / 32) + (pow(B3, 2) / 2) - (k0_ / (2 * m_)) - (b_ * x[0] * (k0_ / m_ - pow(x[38 + i], 2)) / (64 * m_)))) / 48;                                                                                                                                                                                                                                                                           // (row+1, 11)
            values[idx++] = 1 + x[0] * (-(b_ / m_) + 4 * (-(b_ * (1 / 2 - (b_ * x[0] / (64 * m_))) / m_) - (x[0] * C1 * x[37 + i] * B3 / (16 * (x[10 + i] + pow(J_ * R_, 2) / m_))) + (x[0] * pow(B3, 2) / 64) - (x[0] * k0_ / (64 * m_)))) / 48;                                                                                                                                                                                                                                     // (row+1, 19)
            values[idx++] = -1 + x[0] * (-(b_ / m_) + 4 * (-(b_ * (1 / 2 + (b_ * x[0] / (64 * m_))) / m_) + (x[0] * C1 * x[38 + i] * B3 / (16 * (x[11 + i] + pow(J_ * R_, 2) / m_))) - (x[0] * pow(B3, 2) / 64) + (x[0] * k0_ / (64 * m_)))) / 48;                                                                                                                                                                                                                                    // (row+1, 20)
            values[idx++] = x[0] * (-g_ * cos(x[28 + i]) + 4 * ((b_ * g_ * x[0] * cos(x[28 + i]) / (64 * m_)) - (g_ * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / 2) + (g_ * x[0] * C1 * B3 * sin(x[28 + i]) / (32 * (x[10 + i] + pow(J_ * R_, 2) / m_))))) / 48;                                                                                                                                                                                    // (row+1, 28)
            values[idx++] = x[0] * (-g_ * cos(x[29 + i]) + 4 * (-(b_ * g_ * x[0] * cos(x[29 + i]) / (64 * m_)) - (g_ * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / 2) - (g_ * x[0] * C1 * B3 * sin(x[29 + i]) / (32 * (x[11 + i] + pow(J_ * R_, 2) / m_))))) / 48;                                                                                                                                                                                   // (row+1, 29)
            values[idx++] = x[0] * ((2 * x[10 + i] * x[37 + i]) + 4 * (-(b_ * x[0] * x[10 + i] * x[37 + i] / (32 * m_)) + (2 * (1 / 2 + x[0] * (-(c_ * pow(R_, 2) / (m_ * pow(x[10 + i], 2))) - (2 * x[19 + i] / x[10 + i])) / (64 * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])))) * C1 * B3) - (g_ * x[0] * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / 64))) / 48;                                                                                    // (row+1, 37)
            values[idx++] = x[0] * ((2 * x[11 + i] * x[38 + i]) + 4 * (-(b_ * x[0] * x[11 + i] * x[38 + i] / (32 * m_)) + (2 * (1 / 2 - x[0] * (-(c_ * pow(R_, 2) / (m_ * pow(x[11 + i], 2))) - (2 * x[20 + i] / x[11 + i])) / (64 * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])))) * C1 * B3) - (g_ * x[0] * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / 64))) / 48;                                                                                    // (row+1, 38)
            values[idx++] = kt_ * R_ * pow(x[0], 2) * C1 * B3 / (384 * m_ * pow(x[10 + i], 2) * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])));                                                                                                                                                                                                                                                                                                                                            // (row+1, 46)
            values[idx++] = kt_ * R_ * pow(x[0], 2) * C1 * B3 / (384 * m_ * pow(x[11 + i], 2) * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])));                                                                                                                                                                                                                                                                                                                                            // (row+1, 47)

            /*     -----[[3]]-----     */
            values[idx++] = (x[0] * (B1 - B2) / 768) + (x[37 + i] + x[38 + i] + 4 * B3) / 48;                                                                                                         // (row+2, 0)
            values[idx++] = pow(x[0], 2) * B4 / 768;                                                                                                                                                  // (row+2, 10)
            values[idx++] = pow(x[0], 2) * B5 / 768;                                                                                                                                                  // (row+2, 11)
            values[idx++] = -(pow(x[0], 2) * x[37 + i]) / (384 * (x[10 + i] + pow(J_ * R_, 2) / m_));                                                                                                 // (row+1, 19)
            values[idx++] = (pow(x[0], 2) * x[38 + i]) / (384 * (x[11 + i] + pow(J_ * R_, 2) / m_));                                                                                                  // (row+1, 20)
            values[idx++] = -1 - (g_ * pow(x[0], 2) * sin(x[29 + i])) / (768 * (x[11 + i] + pow(J_ * R_, 2) / m_));                                                                                   // (row+2, 29)
            values[idx++] = x[0] * (1 + 4 * (1 / 2 + (x[0] * (-(c_ * pow(R_, 2) / (m_ * pow(x[10 + i], 2))) - (2 * x[19 + i] / x[10 + i]))) / (64 * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i]))))) / 48; // (row+2, 37)
            values[idx++] = x[0] * (1 + 4 * (1 / 2 - (x[0] * (-(c_ * pow(R_, 2) / (m_ * pow(x[11 + i], 2))) - (2 * x[20 + i] / x[11 + i]))) / (64 * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i]))))) / 48; // (row+2, 38)
            values[idx++] = kt_ * R_ * pow(x[0], 2) / (768 * m_ * (1 + (pow(J_ * R_, 2) / (m_ * x[10 + i]))) * pow(x[10 + i], 2));                                                                    // (row+2, 46)
            values[idx++] = kt_ * R_ * pow(x[0], 2) / (768 * m_ * (1 + (pow(J_ * R_, 2) / (m_ * x[11 + i]))) * pow(x[11 + i], 2));                                                                    // (row+2, 47)

            /*     -----[[4]]-----     */
            values[idx++] = ((B1 + B2 + (4 * ((kt_ * R_ * D1 / (m_ * pow(C1, 2))) - (c_ * pow(R_, 2) * B3 / (m_ * pow(C1, 2))) - (g_ * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / C1) - (2 * B3 * A2 / C1))) / (1 + pow(J_ * R_, 2) / (m_ * C1))) / 48) + (x[0] * (((pow(J_ * R_, 2) * (x[19 + i] - x[20 + i]) * ((kt_ * R_ * D1 / (m_ * pow(C1, 2))) - (c_ * pow(R_, 2) * B3 / (m_ * pow(C1, 2))) - (g_ * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / C1) - (2 * B3 * A2 / C1))) / (16 * m_ * pow((1 + pow(J_ * R_, 2) / (m_ * C1)), 2) * pow(C1, 2))) + (4 * ((kt_ * R_ * D2 / (64 * m_ * pow(C1, 2))) - (kt_ * R_ * (x[19 + i] - x[20 + i]) * D1 / (32 * m_ * pow(C1, 3))) - (c_ * pow(R_, 2) * (B1 - B2) / (64 * m_ * pow(C1, 2))) + (c_ * pow(R_, 2) * (x[19 + i] - x[20 + i]) * B3 / (32 * m_ * pow(C1, 3))) + (g_ * (x[19 + i] - x[20 + i]) * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / (64 * pow(C1, 2))) + (g_ * (x[37 + i] - x[38 + i]) * sin(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / (64 * C1)) - ((B3 * A1) / (32 * C1)) - (((B1 - B2) * A2) / (32 * C1)) + (((x[19 + i] - x[20 + i]) * B3 * A2) / (32 * pow(C1, 2)))) / (1 + pow(J_ * R_, 2) / (m_ * C1)))) / 48); // (row+3, 0)
            values[idx++] = (kt_ * R_ * pow(x[0], 2)) / (768 * La_ * m_ * (1 + pow(J_ * R_, 2) / (m_ * C1)) * pow(C1, 2));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             // (row+3, 1)
            values[idx++] = -(kt_ * R_ * pow(x[0], 2)) / (768 * La_ * m_ * (1 + pow(J_ * R_, 2) / (m_ * C1)) * pow(C1, 2));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            // (row+3, 2)
            values[idx++] = x[0] * (B4 + (4 * (-(kt_ * R_ * D1 / (m_ * pow(C1, 3))) - (c_ * pow(R_, 2) * x[0] * B4 / (64 * m_ * pow(C1, 2))) + (c_ * pow(R_, 2) * B3 / (m * pow(C1, 3))) + (g_ * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / (2 * pow(C1, 2))) - (x[0] * B3 * (pow(x[37 + i], 2) - k0_ / m_) / (32 * C1)) - (x[0] * B4 * A1 / (32 * C1)) + (B3 * A2 / pow(C1, 2))) / (1 + pow(J_ * R_, 2) / (m_ * C1))) + ((2 * pow(J_ * R_, 2) * ((kt_ * R_ * D1 / (m_ * pow(C1, 2))) - (c_ * pow(R_, 2) * B3 / (m_ * pow(C1, 2))) - (g_ * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / C1) - (2 * B3 * A2 / C1))) / (m_ * pow((1 + pow(J_ * R_, 2) / (m_ * C1)), 2) * pow(C1, 2)))) / 48;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           // (row+3, 10)
            values[idx++] = x[0] * (B5 + (4 * (-(kt_ * R_ * D1 / (m_ * pow(C1, 3))) - (c_ * pow(R_, 2) * x[0] * B5 / (64 * m_ * pow(C1, 2))) + (c_ * pow(R_, 2) * B3 / (m * pow(C1, 3))) + (g_ * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / (2 * pow(C1, 2))) - (x[0] * B3 * (-pow(x[37 + i], 2) + k0_ / m_) / (32 * C1)) - (x[0] * B5 * A1 / (32 * C1)) + (B3 * A2 / pow(C1, 2))) / (1 + pow(J_ * R_, 2) / (m_ * C1))) + ((2 * pow(J_ * R_, 2) * ((kt_ * R_ * D1 / (m_ * pow(C1, 2))) - (c_ * pow(R_, 2) * B3 / (m_ * pow(C1, 2))) - (g_ * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / C1) - (2 * B3 * A2 / C1))) / (m_ * pow((1 + pow(J_ * R_, 2) / (m_ * C1)), 2) * pow(C1, 2)))) / 48;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          // (row+3, 11)
            values[idx++] = x[0] * (-(2 * x[37 + i] / ((1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * x[10 + i])) + (4 * ((c_ * pow(R_, 2) * x[0] * x[37 + i] / (32 * m_ * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * x[10 + i] * pow(C1, 2))) - (kt_ * R_ * x[0] * D1 / (32 * m_ * pow(C1, 3))) + (c_ * pow(R_, 2) * x[0] * B3 / (32 * m_ * pow(C1, 3))) - (2 * (1 / 2 - b_ * x[0] / (64 * m_)) * B3 / C1) + (g_ * x[0] * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / (64 * pow(C1, 2))) + (x[0] * x[37 + i] * A2 / (16 * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * x[10 + i] * C1)) + ((x[0] * B3 * A2) / (32 * pow(C1, 2)))) / (1 + pow(J_ * R_, 2) / (m_ * C1))) + ((pow(J_ * R_, 2) * x[0] * ((kt_ * R_ * D1 / (m_ * pow(C1, 2))) - (c_ * pow(R_, 2) * B3 / (m_ * pow(C1, 2))) - (g_ * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / C1) - (2 * B3 * A2 / C1))) / (16 * m_ * pow((1 + pow(J_ * R_, 2) / (m_ * C1)), 2) * pow(C1, 2)))) / 48;                                                                                                                                                                                                                                                                                                         // (row+3, 19)
            values[idx++] = x[0] * (-(2 * x[38 + i] / ((1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * x[11 + i])) + (4 * (-(c_ * pow(R_, 2) * x[0] * x[38 + i] / (32 * m_ * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * x[11 + i] * pow(C1, 2))) + (kt_ * R_ * x[0] * D1 / (32 * m_ * pow(C1, 3))) - (c_ * pow(R_, 2) * x[0] * B3 / (32 * m_ * pow(C1, 3))) - (2 * (1 / 2 + b_ * x[0] / (64 * m_)) * B3 / C1) - (g_ * x[0] * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / (64 * pow(C1, 2))) - (x[0] * x[38 + i] * A2 / (16 * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * x[11 + i] * C1)) - ((x[0] * B3 * A2) / (32 * pow(C1, 2)))) / (1 + pow(J_ * R_, 2) / (m_ * C1))) - ((pow(J_ * R_, 2) * x[0] * ((kt_ * R_ * D1 / (m_ * pow(C1, 2))) - (c_ * pow(R_, 2) * B3 / (m_ * pow(C1, 2))) - (g_ * cos(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / C1) - (2 * B3 * A2 / C1))) / (16 * m_ * pow((1 + pow(J_ * R_, 2) / (m_ * C1)), 2) * pow(C1, 2)))) / 48;                                                                                                                                                                                                                                                                                                        // (row+3, 20)
            values[idx++] = x[0] * ((g_ * sin(x[28 + i]) / ((1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * x[10 + i])) + ((4 * ((g_ * x[0] * cos(x[28 + i] * B3 / (32 * C1))) - (c_ * g_ * pow(R_, 2) * x[0] * sin(x[28 + i]) / (64 * m_ * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * x[10 + i] * pow(C1, 2))) + (g_ * sin(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / (2 * C1)) - (g_ * x[0] * sin(x[28 + i]) * A2 / (32 * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * x[10 + i] * C1)))) / (1 + pow(J_ * R_, 2) / (m_ * C1)))) / 48;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         // (row+3, 28)
            values[idx++] = x[0] * ((g_ * sin(x[29 + i]) / ((1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * x[11 + i])) + ((4 * ((g_ * x[0] * cos(x[29 + i] * B3 / (32 * C1))) - (c_ * g_ * pow(R_, 2) * x[0] * sin(x[29 + i]) / (64 * m_ * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * x[11 + i] * pow(C1, 2))) + (g_ * sin(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / (2 * C1)) - (g_ * x[0] * sin(x[29 + i]) * A2 / (32 * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * x[11 + i] * C1)))) / (1 + pow(J_ * R_, 2) / (m_ * C1)))) / 48;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         // (row+3, 29)
            values[idx++] = 1 + x[0] * ((B6 / (1 + pow(J_ * R_, 2) / (m_ * x[10 + i]))) + (4 * (-(kb_ * kt_ * pow(R_, 2) * x[0] / (64 * La_ * m_ * pow(C1, 2))) - (c_ * pow(R_, 2) * (1 / 2 + x[0] * B6 / (64 * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])))) / (m_ * pow(C1, 2))) - (x[0] * x[10 + i] * x[37 + i] * B3 / (16 * C1)) + (g_ * x[0] * sin(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / (64 * C1)) - (2 * (1 / 2 + x[0] * B6 / (64 * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])))) * A2 / C1)) / (1 + pow(J_ * R_, 2) / (m_ * C1)))) / 48;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          // (row+3, 37)
            values[idx++] = -1 + x[0] * ((B7 / (1 + pow(J_ * R_, 2) / (m_ * x[11 + i]))) + (4 * ((kb_ * kt_ * pow(R_, 2) * x[0] / (64 * La_ * m_ * pow(C1, 2))) - (c_ * pow(R_, 2) * (1 / 2 - x[0] * B7 / (64 * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])))) / (m_ * pow(C1, 2))) - (x[0] * x[11 + i] * x[38 + i] * B3 / (16 * C1)) + (g_ * x[0] * sin(((x[28 + i] + x[29 + i]) / 2) + (x[0] * (x[37 + i] - x[38 + i]) / 64)) / (64 * C1)) - (2 * (1 / 2 - x[0] * B7 / (64 * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])))) * A2 / C1)) / (1 + pow(J_ * R_, 2) / (m_ * C1)))) / 48;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          // (row+3, 38)
            values[idx++] = x[0] * ((kt_ * R_ / (m_ * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * pow(x[10 + i], 2))) + (4 * ((kt_ * R_ * (1 / 2 - Ra_ * x[0] / (64 * La_)) / (m_ * pow(C1, 2))) - (c_ * kt_ * pow(R_, 3) * x[0] / (64 * pow(m_, 2) * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * pow(x[10 + i], 2) * pow(C1, 2))) - (kt_ * R_ * x[0] * A2 / (32 * m_ * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * pow(x[10 + i], 2) * C1))) / (1 + pow(J_ * R_, 2) / (m_ * C1)))) / 48;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 // (row+3, 46)
            values[idx++] = x[0] * ((kt_ * R_ / (m_ * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * pow(x[11 + i], 2))) + (4 * ((kt_ * R_ * (1 / 2 + Ra_ * x[0] / (64 * La_)) / (m_ * pow(C1, 2))) + (c_ * kt_ * pow(R_, 3) * x[0] / (64 * pow(m_, 2) * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * pow(x[11 + i], 2) * pow(C1, 2))) + (kt_ * R_ * x[0] * A2 / (32 * m_ * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * pow(x[11 + i], 2) * C1))) / (1 + pow(J_ * R_, 2) / (m_ * C1)))) / 48;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 // (row+3, 47)

            /*     -----[[5]]-----     */
            values[idx++] = (x[0] * (-(Ra_ * D2 / (64 * La_)) - (kb_ * R_ * (B1 - B2) / (64 * La_))) / 12) + ((x[1 + i] / La_) + (x[2 + i] / La_) - (kb_ * R_ * x[37 + i] / La_) - (kb_ * R_ * x[38 + i] / La_) - (Ra_ * x[46 + i] / La_) - (Ra_ * x[47 + i] / La_) + 4 * (((x[1 + i] + x[3 + i]) / (2 * La_)) - (Ra_ * D1 / La_) - (kb_ * R_ * B3 / La_))) / 48; // (row+4, 0)
            values[idx++] = x[0] * (1 / La_ + 4 * (1 / (2 * La_) - (Ra_ * x[0] / (64 * pow(La_, 2))))) / 48;                                                                                                                                                                                                                                                      // (row+4, 1)
            values[idx++] = x[0] * (1 / La_ + 4 * (1 / (2 * La_) - (Ra_ * x[0] / (64 * pow(La_, 2))))) / 48;                                                                                                                                                                                                                                                      // (row+4, 2)
            values[idx++] = -(kb_ * R_ * pow(x[0], 2) * B4) / (768 * La_);                                                                                                                                                                                                                                                                                        // (row+4, 10)
            values[idx++] = -(kb_ * R_ * pow(x[0], 2) * B5) / (768 * La_);                                                                                                                                                                                                                                                                                        // (row+4, 11)
            values[idx++] = -(g_ * kb_ * R_ * pow(x[0], 2) * sin(x[28 + i])) / (768 * La_ * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])) * x[10 + i]);                                                                                                                                                                                                                // (row+4, 28)
            values[idx++] = (g_ * kb_ * R_ * pow(x[0], 2) * sin(x[29 + i])) / (768 * La_ * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])) * x[11 + i]);                                                                                                                                                                                                                 // (row+4, 29)
            values[idx++] = x[0] * (-(kb_ * R_ / La_) + 4 * ((kb_ * R_ * Ra_ * x[0] / (64 * pow(La_, 2))) - (kb_ * R_ * ((1 / 2 + (x[0] * (-(c_ * pow(R_, 2) / (m_ * pow(x[10 + i], 2))) - (2 * x[19 + i] / x[10 + i]))) / (64 * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i]))))) / La_))) / 48;                                                                       // (row+4, 37)
            values[idx++] = x[0] * (-(kb_ * R_ / La_) + 4 * (-(kb_ * R_ * Ra_ * x[0] / (64 * pow(La_, 2))) - (kb_ * R_ * ((1 / 2 - (x[0] * (-(c_ * pow(R_, 2) / (m_ * pow(x[11 + i], 2))) - (2 * x[20 + i] / x[11 + i]))) / (64 * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i]))))) / La_))) / 48;                                                                      // (row+4, 38)
            values[idx++] = 1 + x[0] * (-Ra_ / La_ + 4 * (-(Ra_ * (1 / 2 - Ra_ * x[0] / (64 * La_)) / La_) - ((kb_ * kt_ * pow(R_, 2) * x[0]) / (64 * La_ * m_ * pow(x[10 + i], 2) * (1 + pow(J_ * R_, 2) / (m_ * x[10 + i])))))) / 48;                                                                                                                           // (row+4, 46)
            values[idx++] = -1 + x[0] * (-Ra_ / La_ + 4 * (-(Ra_ * (1 / 2 + Ra_ * x[0] / (64 * La_)) / La_) + ((kb_ * kt_ * pow(R_, 2) * x[0]) / (64 * La_ * m_ * pow(x[11 + i], 2) * (1 + pow(J_ * R_, 2) / (m_ * x[11 + i])))))) / 48;                                                                                                                          // (row+4, 47)
        }
    }

    return true;
} // end of eval_jac_g

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
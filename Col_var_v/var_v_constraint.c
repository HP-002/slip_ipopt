#include "var_header.c";

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

#include "var_header.c";

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
    *obj_value = (x[0]*(pow(x[9],2) + pow(x[1],2) + 2*pow(x[2],2) + 2*pow(x[3],2) + 2*pow(x[4],2) + 2*pow(x[5],2) + 2*pow(x[6],2) + 2*pow(x[7],2) + 2*pow(x[8],2)))/(16.*Ra_);
    
    return true;
} // End of eval_f
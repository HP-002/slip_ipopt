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

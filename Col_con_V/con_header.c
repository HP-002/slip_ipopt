/* Header files */
#include "IpStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>


/* MyUserData contains our Problem */
struct MyUserData {
   ipnumber g_offset[2]; /**< This is an offset for the constraints.  */
   IpoptProblem nlp;   /**< The problem to be solved. Required in intermediate_cb. */
};


/* Callback Functions */

// Start of eval_f
static bool eval_f(
   ipindex     n,             // Number of variables
   ipnumber*   x,             // Variables array 
   bool        new_x,
   ipnumber*   obj_value,     // Function
   UserDataPtr user_data )
{
   assert(n == 4);
   (void) n;

   (void) new_x;
   (void) user_data;

   *obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]; // Function to be minimized

   return true;
} // End of eval_f

// Start of eval_g
static bool eval_g(
   ipindex     n,             //
   ipnumber*   x,
   bool        new_x,
   ipindex     m,
   ipnumber*   g,
   UserDataPtr user_data
)
{
   struct MyUserData* my_data = user_data;

   assert(n == 4);
   (void) n;
   assert(m == 2);
   (void) m;

   (void) new_x;

   g[0] = x[0] * x[1] * x[2] * x[3] + my_data->g_offset[0];
   g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + my_data->g_offset[1];

   return true;
} // End of eval_g

// Start of eval_grad_f
static bool eval_grad_f(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber*   grad_f,
   UserDataPtr user_data
)
{
   assert(n == 4);
   (void) n;

   (void) new_x;
   (void) user_data;

   grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);      // Derivatives of functions w.r.t all variables
   grad_f[1] = x[0] * x[3];
   grad_f[2] = x[0] * x[3] + 1;
   grad_f[3] = x[0] * (x[0] + x[1] + x[2]);                    

   return true;
} // End of eval_grad_f

// eval_jac_g
static bool eval_jac_g(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipindex     m,
   ipindex     nele_jac,
   ipindex*    iRow,
   ipindex*    jCol,
   ipnumber*   values,
   UserDataPtr user_data
)
{
   (void) n;
   (void) new_x;
   (void) m;
   (void) nele_jac;
   (void) user_data;

   if( values == NULL )
   {
      /* return the structure of the jacobian */

      /* this particular jacobian is dense */
      iRow[0] = 0;
      jCol[0] = 0;
      iRow[1] = 0;
      jCol[1] = 1;
      iRow[2] = 0;
      jCol[2] = 2;
      iRow[3] = 0;
      jCol[3] = 3;
      iRow[4] = 1;
      jCol[4] = 0;
      iRow[5] = 1;
      jCol[5] = 1;
      iRow[6] = 1;
      jCol[6] = 2;
      iRow[7] = 1;
      jCol[7] = 3;
   }
   else
   {
      /* return the values of the jacobian of the constraints */

      values[0] = x[1] * x[2] * x[3]; /* 0,0 */
      values[1] = x[0] * x[2] * x[3]; /* 0,1 */
      values[2] = x[0] * x[1] * x[3]; /* 0,2 */
      values[3] = x[0] * x[1] * x[2]; /* 0,3 */

      values[4] = 2 * x[0]; /* 1,0 */
      values[5] = 2 * x[1]; /* 1,1 */
      values[6] = 2 * x[2]; /* 1,2 */
      values[7] = 2 * x[3]; /* 1,3 */
   }

   return true;
} // end of eval_jac_g

// Start of eval_h
static bool eval_h(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber    obj_factor,
   ipindex     m,
   ipnumber*   lambda,
   bool        new_lambda,
   ipindex     nele_hess,
   ipindex*    iRow,
   ipindex*    jCol,
   ipnumber*   values,
   UserDataPtr user_data
)
{
   (void) n;
   (void) new_x;
   (void) m;
   (void) new_lambda;
   (void) user_data;

   if( values == NULL )
   {
      ipindex idx; /* nonzero element counter */
      ipindex row; /* row counter for loop */
      ipindex col; /* col counter for loop */

      /* return the structure. This is a symmetric matrix, fill the lower left
       * triangle only. */

      /* the hessian for this problem is actually dense */
      idx = 0;
      for( row = 0; row < 4; row++ )
      {
         for( col = 0; col <= row; col++ )
         {
            iRow[idx] = row;
            jCol[idx] = col;
            idx++;
         }
      }

      assert(idx == nele_hess);
      (void) nele_hess;
   }
   else
   {
      /* return the values. This is a symmetric matrix, fill the lower left
       * triangle only */

      /* fill the objective portion */
      values[0] = obj_factor * (2 * x[3]); /* 0,0 */

      values[1] = obj_factor * (x[3]); /* 1,0 */
      values[2] = 0; /* 1,1 */

      values[3] = obj_factor * (x[3]); /* 2,0 */
      values[4] = 0; /* 2,1 */
      values[5] = 0; /* 2,2 */

      values[6] = obj_factor * (2 * x[0] + x[1] + x[2]); /* 3,0 */
      values[7] = obj_factor * (x[0]); /* 3,1 */
      values[8] = obj_factor * (x[0]); /* 3,2 */
      values[9] = 0; /* 3,3 */

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
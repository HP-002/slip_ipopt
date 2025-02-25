/* Header files */
#include "IpStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

__uint64_t Ra_ = 1;
__uint64_t b_ = 2;
__uint64_t m_ = 3;
__uint64_t g_ = 4;   // Gravity?
__uint64_t l0_ = 5;
__uint64_t k0_ = 6;
__uint64_t c_ = 7;
__uint8_t R_ = 8;
__uint64_t kt_ = 9;
__uint64_t J_ = 10;
__uint64_t La_ = 11;
__uint64_t kb_ = 12;

/* MyUserData contains our Problem */
struct MyUserData {
   ipnumber g_offset[2]; /**< This is an offset for the constraints.  */
   IpoptProblem nlp;   /**< The problem to be solved. Required in intermediate_cb. */
};


/* Callback Functions */

/* ---------------bool eval_f---------------
 * Method to request the value of the objective function. 
 */
static bool eval_f(
   ipindex     n,             // Number of variables
   ipnumber*   x,             // Variables array 
   bool        new_x,
   ipnumber*   obj_value,     // Function
   UserDataPtr user_data )
{
   assert(n == 55);
   (void) n;

   (void) new_x;
   (void) user_data;
   *obj_value = x[0] * 
               (pow(x[9], 2) + pow(x[1], 2) + (2*pow(x[2], 2)) + (2*pow(x[3], 2)) + (2*pow(x[4], 2)) + 
                  (2*pow(x[5], 2)) + (2*pow(x[6], 2)) + (2*pow(x[7], 2)) + (2*pow(x[8], 2))) / 
               (16 * Ra_); // Function to be minimized

   return true;
} // End of eval_f

/* ---------------bool eval_g---------------
 * Method to request the constraint values.
 */
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

   assert(n == 55);
   (void) n;
   assert(m == 2);   // ????
   (void) m;

   (void) new_x;

   /* Equality Constraints */
   /* ---------------h1--------------- */
   double A1 = 0 - (b_ * x[19] / m) + (b_ * x[20] / m_) + (x[10] * pow(x[37], 2)) - (x[11] * pow(x[38], 2)) - (g_ * sin(x[28])) +( g_ * sin(x[29])) - ((x[10]-l0_) * k0_ / m_) + ((x[11] - l0_) * k0_ / m_);
   double A2 = ((x[19] + x[20])/2) + ((x[0]/64)*(A1));

   double B1 = ((0 - ((c_ * pow(R_, 2) * x[37]) / (m * pow(x[10], 2))) - ((2 * x[19] * x[37]) / (x[10])) + ((kt_ * R_ * x[46]) / (m_ * pow(x[10], 2))) - ((g_ * cos(x[28])) / (x[10]))) / (1 + (pow(J_*R_, 2) / (m_ * x[10]))));
   double B2 = ((0 - ((c_ * pow(R_, 2) * x[38]) / (m * pow(x[11], 2))) - ((2 * x[20] * x[38]) / (x[11])) + ((kt_ * R_ * x[47]) / (m_ * pow(x[11], 2))) - ((g_ * cos(x[29])) / (x[11]))) / (1 + (pow(J_*R_, 2) / (m_ * x[11]))));

   double C1 = ((x[10]+x[11]) / 2) + ((x[0] / 64) * (x[19]-x[20]));

   double D1 =  ((x[46]+x[47]) / 2) + ((x[0] / 64) * ((x[1] / La_) - (x[2] / La_) - (kb_ * R_ * x[37] / La_) + (kb_ * R_ * x[38] / La_) - (Ra_ * x[46] / La_) + (Ra_ * x[47] / La_)));

   g[0] = x[10] - x[11] + ((x[0] / 48) * (x[19] + x[20] + 4 * A2));
   
   g[1] = x[19] - x[20] + (x[0] / 48) * (A1 + 4 * (( ((x[10]+x[11]) / 2) + ((x[0] / 64) * (x[19] - x[20])) ) * pow(( ((x[37] + x[38]) / 2) + ((x[0] / 64) * ( (B1) - (B2) )) ), 2) -
         ( g_ * sin(((x[28]+x[29])/2) + ((x[0]/64)*(x[37]-x[38]))) ) - ( (0 - l0_ + ( (x[10]+x[11]) / 2 ) + ( (x[0] / 64) * (x[19]-x[20]) ) ) * k0_ / m_ ) - ( b_ * (A2) / m_)
      )
   );

   g[2] = x[28] - x[29] + ((x[0] / 48) * (x[37] + x[38] + (4 * ( ((x[37]+x[38]) / 2) + ( (x[0] / 64) * (B1 - B2) )))));

   g[3] = x[37] - x[38] + ((x[0] / 48) * ( B1 + B2 + ( 4 * ( ((kt_ * R_ * (D1)) / (m * pow(C1, 2))) - 
      ( (c_ * pow(R_, 2) * (((x[37]+x[38]) / 2) + ((x[0] / 64) * (B1 - B2))) ) / (m * pow(C1, 2)) ) - ( (g_ * cos(((x[28]+x[29]) / 2) + ((x[0] / 64) * (x[37]-x[38])))) / C1 ) - 
      ( (2 * (((x[37]+x[38]) / 2) + ((x[0] / 64) * (B1 - B2))) * A2) / C1 ) ) / (1 + (pow(J_*R_, 2) / (m_ * C1))) ) )
   );

   g[4] = x[46] - x[47] + ( (x[0] / 48) * ( (x[1] / La_) + (x[2] / La_) - (kb_ * R_ * (x[37] + x[38]) / La_) - (R_ * (x[46] + x[47]) / La_) + (4 * ( ((x[1]+x[2]) / (2 * La_)) - (Ra_ * D1 / La_) - ( (kb_ * R_ * (((x[37]+x[38]) / 2) + ((x[0] / 64) * (B1 - B2)))) / La_) )) ) );



   /* ---------------h2--------------- */
   A1 = 0 - (b_ * x[20] / m) + (b_ * x[21] / m_) + (x[11] * pow(x[38], 2)) - (x[12] * pow(x[39], 2)) - (g_ * sin(x[29])) +( g_ * sin(x[30])) - ((x[11]-l0_) * k0_ / m_) + ((x[12] - l0_) * k0_ / m_);
   A2 = ((x[20] + x[21])/2) + ((x[0]/64)*(A1));

   B1 = ((0 - ((c_ * pow(R_, 2) * x[38]) / (m * pow(x[11], 2))) - ((2 * x[20] * x[38]) / (x[11])) + ((kt_ * R_ * x[47]) / (m_ * pow(x[11], 2))) - ((g_ * cos(x[29])) / (x[11]))) / (1 + (pow(J_*R_, 2) / (m_ * x[11]))));
   B2 = ((0 - ((c_ * pow(R_, 2) * x[39]) / (m * pow(x[12], 2))) - ((2 * x[21] * x[39]) / (x[12])) + ((kt_ * R_ * x[48]) / (m_ * pow(x[12], 2))) - ((g_ * cos(x[30])) / (x[12]))) / (1 + (pow(J_*R_, 2) / (m_ * x[12]))));

   C1 = ((x[11]+x[12]) / 2) + ((x[0] / 64) * (x[20]-x[21]));

   D1 =  ((x[47]+x[48]) / 2) + ((x[0] / 64) * ((x[2] / La_) - (x[3] / La_) - (kb_ * R_ * x[38] / La_) + (kb_ * R_ * x[39] / La_) - (Ra_ * x[47] / La_) + (Ra_ * x[48] / La_)));
   
   g[5] = x[11] - x[12] + ((x[0] / 48) * (x[20] + x[21] + 4 * A2));

   g[6] = x[20] - x[21] + (x[0] / 48) * (A1 + 4 * (( ((x[11]+x[12]) / 2) + ((x[0] / 64) * (x[20] - x[21])) ) * pow(( ((x[38] + x[39]) / 2) + ((x[0] / 64) * (B1 - B2)) ), 2) - 
         ( g_ * sin(((x[29]+x[30])/2) + ((x[0]/64)*(x[38]-x[39]))) ) - ( (0 - l0_ + ( (x[11]+x[12]) / 2 ) + ( (x[0] / 64) * (x[20]-x[21]) ) ) * k0_ / m_ ) - ( b_ * (A2) / m_)
      )
   );

   g[7] = x[29] - x[30] + ((x[0] / 48) * (x[38] + x[39] + (4 * ( ((x[38]+x[39]) / 2) + ( (x[0] / 64) * (B1 - B2) )))));

   g[8] = x[38] - x[39] + ((x[0] / 48) * ( B1 + B2 + ( 4 * ( ((kt_ * R_ * (D1)) / (m * pow(C1, 2))) - 
      ( (c_ * pow(R_, 2) * (((x[38]+x[39]) / 2) + ((x[0] / 64) * (B1 - B2))) ) / (m * pow(C1, 2)) ) - ( (g_ * cos(((x[29]+x[30]) / 2) + ((x[0] / 64) * (x[38]-x[39])))) / C1 ) - 
      ( (2 * (((x[38]+x[39]) / 2) + ((x[0] / 64) * (B1 - B2))) * A2) / C1 ) ) / (1 + (pow(J_*R_, 2) / (m_ * C1))) ) )
   );

   g[9] = x[47] - x[48] + ( (x[0] / 48) * ( (x[2] / La_) + (x[3] / La_) - (kb_ * R_ * (x[38] + x[39]) / La_) - (Ra_ * (x[47] +x[48]) / La_) + (4 * ( ((x[2]+x[3]) / (2 * La_)) - (Ra_ * D1 / La_) - ( (kb_ * R_ * (((x[38]+x[39]) / 2) + ((x[0] / 64) * (B1 - B2)))) / La_) )) ) );
   


   /* ---------------h3--------------- */
   A1 = 0 - (b_ * x[21] / m) + (b_ * x[22] / m_) + (x[12] * pow(x[39], 2)) - (x[13] * pow(x[40], 2)) - (g_ * sin(x[30])) +( g_ * sin(x[31])) - ((x[12]-l0_) * k0_ / m_) + ((x[13] - l0_) * k0_ / m_);
   A2 = ((x[21] + x[22])/2) + ((x[0]/64)*(A1));

   B1 = ((0 - ((c_ * pow(R_, 2) * x[39]) / (m * pow(x[12], 2))) - ((2 * x[21] * x[39]) / (x[12])) + ((kt_ * R_ * x[48]) / (m_ * pow(x[12], 2))) - ((g_ * cos(x[30])) / (x[12]))) / (1 + (pow(J_*R_, 2) / (m_ * x[12]))));
   B2 = ((0 - ((c_ * pow(R_, 2) * x[40]) / (m * pow(x[13], 2))) - ((2 * x[22] * x[40]) / (x[13])) + ((kt_ * R_ * x[49]) / (m_ * pow(x[13], 2))) - ((g_ * cos(x[31])) / (x[13]))) / (1 + (pow(J_*R_, 2) / (m_ * x[13]))));

   C1 = ((x[12]+x[13]) / 2) + ((x[0] / 64) * (x[21]-x[22]));

   D1 =  ((x[48]+x[49]) / 2) + ((x[0] / 64) * ((x[3] / La_) - (x[4] / La_) - (kb_ * R_ * x[39] / La_) + (kb_ * R_ * x[40] / La_) - (Ra_ * x[48] / La_) + (Ra_ * x[49] / La_)));

   g[10] = x[12] - x[13] + ((x[0] / 48) * (x[21] + x[22] + 4 * A2));

   g[11] = x[21] - x[22] + (x[0] / 48) * (A1 + 4 * (( ((x[12]+x[13]) / 2) + ((x[0] / 64) * (x[21] - x[22])) ) * pow(( ((x[39] + x[40]) / 2) + ((x[0] / 64) * (B1 - B2)) ), 2) - 
         ( g_ * sin(((x[30]+x[31])/2) + ((x[0]/64)*(x[39]-x[40]))) ) - ( (0 - l0_ + ( (x[12]+x[13]) / 2 ) + ( (x[0] / 64) * (x[21]-x[22]) ) ) * k0_ / m_ ) - ( b_ * (A2) / m_)
      )
   );

   g[12] = x[30] - x[31] + ((x[0] / 48) * (x[39] + x[40] + (4 * ( ((x[39]+x[40]) / 2) + ( (x[0] / 64) * (B1 - B2) )))));

   g[13] = x[39] - x[40] + ((x[0] / 48) * ( B1 + B2 + ( 4 * ( ((kt_ * R_ * (D1)) / (m * pow(C1, 2))) - 
      ( (c_ * pow(R_, 2) * (((x[39]+x[40]) / 2) + ((x[0] / 64) * (B1 - B2))) ) / (m * pow(C1, 2)) ) - ( (g_ * cos(((x[30]+x[31]) / 2) + ((x[0] / 64) * (x[39]-x[40])))) / C1 ) - 
      ( (2 * (((x[39]+x[40]) / 2) + ((x[0] / 64) * (B1 - B2))) * A2) / C1 ) ) / (1 + (pow(J_*R_, 2) / (m_ * C1))) ) )
   );

   g[14] = x[48] - x[49] + ( (x[0] / 48) * ( (x[3] / La_) + (x[4] / La_) - (kb_ * R_ * (x[39] + x[40]) / La_) - (Ra_ * (x[48] +x[49]) / La_) + (4 * ( ((x[3]+x[4]) / (2 * La_)) - (Ra_ * D1 / La_) - ( (kb_ * R_ * (((x[39]+x[40]) / 2) + ((x[0] / 64) * (B1 - B2)))) / La_) )) ) );
   


   /* ---------------h4--------------- */
   A1 = 0 - (b_ * x[22] / m) + (b_ * x[23] / m_) + (x[13] * pow(x[40], 2)) - (x[14] * pow(x[41], 2)) - (g_ * sin(x[31])) +( g_ * sin(x[32])) - ((x[13]-l0_) * k0_ / m_) + ((x[14] - l0_) * k0_ / m_);
   A2 = ((x[22] + x[23])/2) + ((x[0]/64)*(A1));

   B1 = ((0 - ((c_ * pow(R_, 2) * x[40]) / (m * pow(x[13], 2))) - ((2 * x[22] * x[40]) / (x[13])) + ((kt_ * R_ * x[49]) / (m_ * pow(x[13], 2))) - ((g_ * cos(x[31])) / (x[13]))) / (1 + (pow(J_*R_, 2) / (m_ * x[13]))));
   B2 = ((0 - ((c_ * pow(R_, 2) * x[41]) / (m * pow(x[14], 2))) - ((2 * x[23] * x[41]) / (x[14])) + ((kt_ * R_ * x[50]) / (m_ * pow(x[14], 2))) - ((g_ * cos(x[32])) / (x[14]))) / (1 + (pow(J_*R_, 2) / (m_ * x[14]))));

   C1 = ((x[13]+x[14]) / 2) + ((x[0] / 64) * (x[22]-x[23]));

   D1 =  ((x[49]+x[50]) / 2) + ((x[0] / 64) * ((x[4] / La_) - (x[5] / La_) - (kb_ * R_ * x[40] / La_) + (kb_ * R_ * x[41] / La_) - (Ra_ * x[49] / La_) + (Ra_ * x[50] / La_)));

   g[15] = x[13] - x[14] + ((x[0] / 48) * (x[22] + x[23] + 4 * A2));

   g[16] = x[22] - x[23] + (x[0] / 48) * (A1 + 4 * (( ((x[13]+x[14]) / 2) + ((x[0] / 64) * (x[22] - x[23])) ) * pow(( ((x[40] + x[41]) / 2) + ((x[0] / 64) * (B1 - B2)) ), 2) - 
         ( g_ * sin(((x[31]+x[32])/2) + ((x[0]/64)*(x[40]-x[41]))) ) - ( (0 - l0_ + ( (x[13]+x[14]) / 2 ) + ( (x[0] / 64) * (x[22]-x[23]) ) ) * k0_ / m_ ) - ( b_ * (A2) / m_)
      )
   );

   g[17] = x[31] - x[32] + ((x[0] / 48) * (x[40] + x[41] + (4 * ( ((x[40]+x[41]) / 2) + ( (x[0] / 64) * (B1 - B2) )))));

   g[18] = x[40] - x[41] + ((x[0] / 48) * ( B1 + B2 + ( 4 * ( ((kt_ * R_ * (D1)) / (m * pow(C1, 2))) - 
      ( (c_ * pow(R_, 2) * (((x[40]+x[41]) / 2) + ((x[0] / 64) * (B1 - B2))) ) / (m * pow(C1, 2)) ) - ( (g_ * cos(((x[31]+x[32]) / 2) + ((x[0] / 64) * (x[40]-x[41])))) / C1 ) - 
      ( (2 * (((x[40]+x[41]) / 2) + ((x[0] / 64) * (B1 - B2))) * A2) / C1 ) ) / (1 + (pow(J_*R_, 2) / (m_ * C1))) ) )
   );

   g[19] = x[49] - x[50] + ( (x[0] / 48) * ( (x[4] / La_) + (x[5] / La_) - (kb_ * R_ * (x[40] + x[41]) / La_) - (Ra_ * (x[49] + x[50]) / La_) + (4 * ( ((x[4]+x[5]) / (2 * La_)) - (Ra_ * D1 / La_) - ( (kb_ * R_ * (((x[40]+x[41]) / 2) + ((x[0] / 64) * (B1 - B2)))) / La_) )) ) );
   




   // g[0] = x[0] * x[1] * x[2] * x[3] + my_data->g_offset[0];
   // g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + my_data->g_offset[1];

   return true;
} // End of eval_g

/* ---------------bool eval_grad_f---------------
 * Method to request the gradient of the objective function.
 */
static bool eval_grad_f(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber*   grad_f,
   UserDataPtr user_data
)
{
   assert(n == 55);
   (void) n;

   (void) new_x;
   (void) user_data;

   grad_f[0] = (pow(x[9], 2) + pow(x[1], 2) + (2*pow(x[2], 2)) + (2*pow(x[3], 2)) + (2*pow(x[4], 2)) + (2*pow(x[5], 2)) + (2*pow(x[6], 2)) + (2*pow(x[7], 2)) + (2*pow(x[8], 2))) / (16 * Ra_);      // Derivatives of functions w.r.t all variables
   grad_f[1] = x[0] * x[1] / (8 * Ra_);                    
   grad_f[2] = x[0] * x[2] / (4 * Ra_);                    
   grad_f[3] = x[0] * x[3] / (4 * Ra_);                    
   grad_f[4] = x[0] * x[4] / (4 * Ra_);                    
   grad_f[5] = x[0] * x[5] / (4 * Ra_);                    
   grad_f[6] = x[0] * x[6] / (4 * Ra_);                    
   grad_f[7] = x[0] * x[7] / (4 * Ra_);                    
   grad_f[8] = x[0] * x[8] / (4 * Ra_);                    
   grad_f[9] = x[0] * x[9] / (8 * Ra_);                    

   for(int i = 10; i < n; i++) {
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

/* ---------------bool eval_h---------------
 * Method to request either the sparsity structure or the values of the Hessian of the Lagrangian.
 *
 * The Jacobian is the matrix of derivatives where the derivative of constraint function g_i with respect
 * to variable x_j is placed in row i and column j.
 */
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
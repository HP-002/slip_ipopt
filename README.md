# SLIP_IPopt

## Project Structure
```
src/
├── Col_con_v/
│   └── # C source files for the problem with constant voltage
│
├── Col_var_v/
│   └── # C source files for the problem with variable voltage
│
├── Equations/
│   ├── Col_con_v/
│   │   ├── parsed/
│   │   │   └── # Parsed text files containing equations formatted for C
│   │   └── # Raw equation text files exported from Mathematica
│   │
│   ├── Col_var_v/
│   │   ├── parsed/
│   │   │   └── # Parsed text files containing equations formatted for C
│   │   └── # Raw equation text files exported from Mathematica
│   │
│   └── Tests/
│       └── # Test files used to validate parser.py functionality
│
├── Mathematica/
│   ├── Col_con_v/
│   │   └── # Mathematica notebooks and scripts for constant voltage problem
│   │
│   └── Col_var_v/
│       └── # Mathematica notebooks and scripts for variable voltage problem
│
├── parser.py          # Python script to parse Mathematica equations into C-style equations
├── Makefile           # Build configuration file
└── README.md          # Project documentation

```

## Execution Steps
1) Print the equations from the Mathematica to text files in <code>src/Equations/</code> directory. This includes cost, cost gradient, equality constraints, inequality constraints, constraints jacobian, & constraint hessian.
2) Run the <code>parser.py</code> to parse the equations in text files in <code>src/Equations/</code> and output them in the same format in <code>src/Equations/parsed/</code>. The structure for each text file of equations is given below.
3) Copy the parsed equations in the C files for each problem in their respective function of IPopt. (TODO)
4) Execute IPopt using the Makefile (TODO)

## Equation Text Files Structure
### [Equations/Col_var_v](Equations/Col_var_v):
#### [cost_var_v_text.txt](Equations/Col_var_v/cost_var_v_text.txt):
Cost function in a single line of text. This is the objective function.

#### [costgrad_var_v_text.txt](Equations/Col_var_v/costgrad_var_v_text.txt):
Gradient of the cost function. Each line contains a single equation. Gradient is the derivative of the cost function w.r.t to all the variables. Thus, it is 55 x 1 Matrix.

#### [equality_var_v_text.txt](Equations/Col_var_v/equality_var_v_text.txt):
Contains equality constraints of the problem. There are 40 equality constraints. Each line contains a single constraint.


#### [inequality_var_v_text.txt](Equations/Col_var_v/inequality_var_v_text.txt):
Contains inequality constraint of the problem. There is only 1 inequality constraint present on a single line.

#### [eq_jacob_var_v_text.txt](Equations/Col_var_v/eq_jacob_var_v_text.txt):
Contains Jacobian of the equality constraints. The Jacobian is the matrix of derivatives where the derivative of constraint function g_i with respect to variable x_j is placed in row i and column j. There are 40 equality constraints and 55 variables, therefore, the jacobian is a 40 x 55 matrix. The file contains each equation on a single line.

The constraints are present in 8 batches (h1, h2,...,h8) of 5 (hi[[1]], hi[[2]],...hi[[5]]). In the text file, the jacobian is present as:
```
h_1
1
jacobian of h_1[[1]] w.r.t x[0] 
jacobian of h_1[[1]] w.r.t x[1] 
...
jacobian of h_1[[1]] w.r.t x[54]
2
jacobian of h_1[[2]] w.r.t x[0] 
jacobian of h_1[[2]] w.r.t x[1] 
...
jacobian of h_1[[2]] w.r.t x[54]
...
...
5
jacobian of h_1[[1]] w.r.t x[0] 
jacobian of h_1[[1]] w.r.t x[1] 
...
jacobian of h_1[[1]] w.r.t x[54]
h2
1
jacobian of h_1[[1]] w.r.t x[0] 
...
```
#### [Equations/Col_var_v/parsed](Equations/Col_var_v/parsed):
The parsed text files follow the same structure as above.
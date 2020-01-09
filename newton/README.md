# Newton

## newton.py

The file **newton.py** contains the implementation of the Newton algorithm to find the zero of a one-dimensional function. Thus the algorithm solves the equation `f(x)=0`.

For a description of the algorithm, see [wikipedia](https://en.wikipedia.org/wiki/Newton%27s_method). 

The algorithm here implemented does not require the number of iterations to be performed, because it stops when the function calculated at the n-th iteration is closer to zero than the precision parameter `eps` (set by default at `eps=1E-16`, i.e. machine precision). 

To call the algorithm, use the following line:

`xn = newton(f, Df, x0, eps=1e-16, max_n=50)`

As input parameters, the algortihm needs:
* f: the function of which the zero should be calculated;
* Df: the derivative of the function f;
* x0: a starting guess of the zero;
* eps: the precision of the algorithm (default value `eps=1e-16`);
* max_n: maximum iterations to be performed (default value `max_n=50`).

As output, the function return:
* xn: the zero of the function

Notice that an error arises if the derivative is zero in the initial guess point or in the iterated values. In this case, the function should print `Zero derivative` and return no value. 

If instead the number of iterations exceeds the number of maximum iterations `max_n`, then the function prints `Maximum number of iterations reached` and returns no value.
It is known from literature that the Newton algorithm has quadratic convergence: this means that the algorithm should find the zero at machine precision in 8/9 steps. 
Thus, if this error arises, be careful you implemented correctly the function and you have chosen an appropriate initial guess. 

### newton_test.py

The file **newton_test.py** contains the tests of the one-dimensional Newton algorithm implemented in **newton.py**. To perform the tests, go to the newton folder and digit the following command line in the iPython console:

`!pytest newton_test.py`

We here give a description of the tests performed.

* `test_one` applies the algorithm to the function `f(x)=x**2-4` with initial guesses `x0=3.0` and `x0=-1.5`; the test checks if the algorithm correctly finds as solution the zeros `x=2.0` and `x=-2.0`.
* `test_two` applies the algorithm to the function `f(x)=rx-x**2` for `r=3` and checks if the algorithm correctly finds, starting from `x0=2.0`, the zero `x=3`.
* `test_three` applies the algorithm to the function `f(x)=sin(x)` and checks if the algortihm, starting from initial guess `x0=0.5`, finds the zero `x=0`.

## newton2.py

The file **newton2_test.py** contains the two-dimensional Newton algorithm. Thus the algorithm solves the following system:
```
f1(x,y)=0
f2(x,y)=0
```
For a theoretical description of the two-dimensional Newton algortihm, see [mathfaculty.fullerton.edu](http://mathfaculty.fullerton.edu/mathews/n2003/FixPointNewtonMod.html).

To call the algorithm, digit the following line:

`pk = newton2(f,Jf,p0,eps=1e-8,max_iter=20)`

The input parameters needed for the function are:
* f: the list of two functions that form the system `f=[f1,f2]`;
* Jf: jacobian matrix, i.e. the array containing the derivatives of the functions with respect x and y;
* p0: two dimensional list containing the initial guess `p0=[x0,y0]`;
* eps: precision parameter (by default, `eps=1e-8`);
* max_iter: maximum number of iterations to be performed (by default 6`max_iter=20`).

As output, the function returns:
* pk: a two dimensional list containing the solution of the system.

Notice that if the maximum number of iterations is reached, the function prints `Maximum number of iterations reached`.

If the Jacobian matrix is singular in the initial guess point or in some iterated values, the function prints `Singular matrix: determinant 0`.

If one of the previous errors arise (maximum number of interations reached or singular matrix), then the output of the algorithm is a False boolean variable.

### newton2_test.py

The file **newton2_test.py** contains the tests of the bidimensional Newton algorithm implemented in **newton2.py**. To perform the test, go to the newton folder and digit the following instruction in the iPython console:

`!pytest newton2_test.py`

We give here a description of the tests performed.

* `test_one` applies the algorithm to the set of functions `f1(x,y)=1-4*x+2*x**2-2*y**3` and `f2(x,y)=-4+x**4+4*y+4*y**4` for which an analitical solution is known: with 6 decimals: `(0.06177,0.724491)`. The test checks if the algorithm finds the analitical solution (within 6 decimals).
* `test_two` applies the algorithm to the system composed by `f1(x,y)=x` and `f2(x,y)=y` and checks if the algorithm finds `(0,0)` as solution.
* `test_three` applies the algorithm to the system composed by `f1(x,y)=f2(x,y)=0` that should return a singular matrix. The test checks if the algorithm returns a `False` boolean value as it should do. 
* `test_four` applies the algorithm to the system `f1(x,y)=x**2+y**2+1` and `f2(x,y)=2*x` that should have no zeros, thus the algorithm should reach the maximum number of iterations. The test checks if the function returns a `False` boolean variable as it should do in this case.

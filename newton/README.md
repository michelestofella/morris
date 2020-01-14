# Newton

The present folder contains two scripts implementing the monodimensional (`newton.py`) and bidimensional (`newton2.py`) newton algorithm to find the zeros of a function (or a system of functions). 

Both the algorithm have been tested in the files newton_test.py (monodimensional) and newton2_test.py (bidimensional); the libraries pytest and hypothesis were used in order to develop testing strategies. 

A brief content of each file contained in the folder is given below. 

## newton.py

The file `newton.py` contains the implementation of the Newton algorithm to find the zero of a one-dimensional function. Thus the algorithm solves the equation `f(x)=0`.

For a description of the algorithm, see [wikipedia](https://en.wikipedia.org/wiki/Newton%27s_method). 

The algorithm here implemented does not require the number of iterations to be performed, because it stops when the function calculated at the n-th iteration is closer to zero than the precision parameter `eps` (set by default at `eps=1E-14`, i.e. machine precision). 

To call the algorithm, use the following line:

`xn = newton(f, Df, x0, eps=1e-14, max_n=100)`

As input parameters, the algortihm needs:
* `f`: the function of which the zero should be calculated;
* `Df`: the derivative of the function f;
* `x0`: a starting guess of the zero;
* `eps`: the precision of the algorithm (default value `eps=1e-14`);
* `max_n`: maximum iterations to be performed (default value `max_n=100`).

As output, the function return:
* `xn`: the zero of the function

An error (Exception) arises if one of the following conditions is reached:
* the derivative is zero 
* the maximum number of iterations is reached

### newton_test.py

The file `newton_test.py` contains the tests of the one-dimensional Newton algorithm implemented in `newton.py`. To perform the tests, go to the `newton` folder and digit the following command line:

`pytest newton_test.py`

We here give a description of the tests performed.

* `test_parabola_without_constant_terms` applies the algorithm to the function `f(x)=x**2` and checks if the algorithm correctly finds as solution the zero `x=0`. The strategy here used is to run the test using different initial conditions.
* `test_parabola_with_constant_terms` applies the algorithm to the function `f(x)=x**2-r` and checks if the algorithm correctly finds the zeros `x=sqrt(r)` or `x=-sqrt(r)`, depending on the initial guess. The strategy here implemented is to run the algorithm using different values of the initial guess: since the function is symmetric, the algorithm will find the positive zero `x=sqrt(r)` when the initial guess is positive and viceversa. Furthermore, the test is applied using different values of the parameter r that describes the function.
* `test_exception` tests if the algorithm raises an error when considering a function with no zeros. 

## newton2.py

The file `newton2.py` contains the two-dimensional Newton algorithm. Thus the algorithm solves the following system:
```
f1(x,y)=0
f2(x,y)=0
```
For a theoretical description of the two-dimensional Newton algortihm, see [mathfaculty.fullerton.edu](http://mathfaculty.fullerton.edu/mathews/n2003/FixPointNewtonMod.html).

To call the algorithm, digit the following line:

`pk = newton2(f,Jf,p0,eps=1e-8,max_iter=20)`

The input parameters needed for the function are:
* `f`: the list of two functions that form the system `f=[f1,f2]`;
* `Jf`: jacobian matrix, i.e. the array containing the derivatives of the functions with respect x and y;
* `p0`: two dimensional list containing the initial guess `p0=[x0,y0]`;
* `eps`: precision parameter (by default, `eps=1e-8`);
* `max_iter`: maximum number of iterations to be performed (by default 6`max_iter=20`).

As output, the function returns:
* `pk`: a two dimensional list containing the solution of the system.

An error (Exception) arises if one of the following conditions is reached:
* the determinant of the jacobian is zero
* the maximum number of iterations has been reached

### newton2_test.py

The file `newton2_test.py` contains the tests of the bidimensional Newton algorithm implemented in `newton2.py`. To perform the test, go to the `newton` folder and digit the following instruction:

`pytest newton2_test.py`

We give here a description of the tests performed.

* `test_unique_solution` considers a system where a unique solution exists, namely `f1(x,y)=x` and `f2(x,y)=y`, and checks if the correct solution is returned by the algorithm. The strategy here applied is to use different initial conditions.
* `test_two_possible_solution` considers a system where two zeros exist, namely `f1(x,y)=x-y` and `f2(x,y)=y**2-r`, and checks that both can be reached starting from a proper starting guess. The strategy is to use different initial conditions in order to reach both zeros. The parameter r that defines function f2 is also changed.
* `test_zero_determinant_exception` tests if the algorithm raises an exception when the determinant of the jacobian matrix is zero. The algorithm is applied to a particular set of function where each function is a constant that is varied within the testing strategy. 
* `test_max_iterations_exception` tests if the algorithm raises an exception when considering a system of function that has no solutions. 

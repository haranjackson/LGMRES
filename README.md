# LGMRES

The LGMRES algorithm is designed to avoid the slowing of convergence in restarted GMRES, due to alternating residual vectors.
Typically, it often outperforms GMRES(m) of comparable memory requirements by some measure, or at least is not much worse.

Another advantage in this algorithm is that you can supply it with ‘guess’ vectors in the `outer_v` argument that augment the Krylov subspace.
If the solution lies close to the span of these vectors, the algorithm converges faster.
This can be useful if several very similar matrices need to be inverted one after another, such as in Newton-Krylov iteration where the Jacobian matrix often changes little in the nonlinear steps.

This repository contains a C++ implementation of the LGMRES algorithm, with Python bindings. It matches SciPy's implementation to machine precision, but is around 30 times faster.

## Installation (Linux)

If you have Python 3 installed, you can simply add `LGMRES/bindings` to your path, and import the precompiled module with `import LGMRES`.
Alternatively, use the following commands to compile the Python bindings yourself (making sure you start in the top LGMRES directory):

```sh
$ sudo apt-get install python3-dev
$ cd bindings
$ bash compile_bindings.sh
```

## Python Usage

The LGMRES module has one function, solve:

`LGMRES.solve(A, b, x0=None, M=None, tol=1e-05, maxiter=1000, inner_m=30, outer_k=3, outer_v=[])`

### Parameters:	
* A : array
 The real or complex N-by-N matrix of the linear system.
* b : array
Right hand side of the linear system. Has shape (N,).
* x0 : array (optional)
Starting guess for the solution.
* M : array (optional)
Preconditioner for A. The preconditioner should approximate the inverse of A.
Effective preconditioning dramatically improves the rate of convergence, which implies that fewer iterations are needed to reach a given error tolerance.
* tol : float (optional)
Tolerance to achieve. The algorithm terminates when either the relative or the absolute residual is below tol.
* maxiter : int (optional)
Maximum number of iterations. Iteration will stop after maxiter steps even if the specified tolerance has not been achieved.
* inner_m : int (optional)
Number of inner GMRES iterations per each outer iteration.
* outer_k : int (optional)
Number of vectors to carry between inner GMRES iterations. Good values are in the range of 1...3.
However, note that if you want to use the additional vectors to accelerate solving multiple similar problems, larger values may be beneficial.
* outer_v : list of arrays (optional)
List containing vectors used to augment the Krylov subspace, and carried between inner GMRES iterations.
This parameter is modified in-place by lgmres, and can be used to pass “guess” vectors in and out of the algorithm when solving similar problems.

### Returns:	
* x : array
The converged solution.

### NOTE
Do not call `solve` with `x0=None` or `M=None` explicitly.
If no initial guess or preconditioner are required, simply leave out these arguments.

## C++ Usage

A class similar to the `Solution` class found in `lgmres.h` should be constructed with public functions `call` and `prec`
(corresponding to the actions of the linear operator required to be solved, and its preconditioner, respectively).
Note that although these functions must take and return a `Vec` (defined in `types.h`), they do not have to take the form of a matrix multiplication;
they may be any abstract linear operator.
Your custom class can then be passed to the `lgmres` function in the same manner as the `Solution` class.
Your code must be compiled against the eigen3 and pybind11 libraries, found in the `include` folder.
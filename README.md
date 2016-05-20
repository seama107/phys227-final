# PHYS227 Final

**Author:** Michael Seaman

[![Build Status](https://travis-ci.org/seama107/phys227-final.svg?branch=master)](https://travis-ci.org/seama107/phys227-final)

**Due date:** 2016/05/20

## Specification

Problem 1. Implement a python class with constructor Rossler(c, dt=0.001, T0=250, T=500) (10pt)
that implements a 4th-order Runge-Kutta integrator for the Eqs. (1). The class constructor should initialize
the following attributes: self.dt (float64), self.T (float64), self.t (numpy array), self.x (numpy array),
self.y (numpy array), self.z (numpy array). Just after initialization, self.t should store an evenly spaced
set of time points from 0 to T with spacing self.dt, while self.x, self.y, and self.z should store arrays
of zeros of the same length. This class should have a method run(self) that integrates Eqs. (1) from
the initial conditions ($ x_0 $, $ y_0 $, $ z_0 $) = (0, 0, 0) at t = 0 to the final time t = T and stores the results in the
corresponding attributes. This class should also have seven plotting methods that use matplotlib to visualize
the solution (10pt): plotx(self) (x vs. t), ploty(self) (y vs. t), plotz(self) (z vs. t), plotxy(self)
(y vs. x), plotyz(self) (z vs. y), plotxz(self) (z vs. x), and plotxyz(self) (3D x-y-z plot). The
plots vs. t should show the range $ t \in [0, T] $, while the 2D and 3D parametric plots should only include the
steady-state time range $ t \in [T_0, T] $. Plots should be labeled clearly, with the x and y ranges [−12, 12], and
the z range [0, 25]. Show that all functionality works in your notebook using c = 2.

Problem 2. Explore the onset of chaos (25pt). In your notebook, show a time plot x(t), a 2D parametric
plot y(x), and a 3D parametric plot z(x, y) for each of the following c values: 2, 3, 4, 4.15, 4.2, and 5.7. Show
other plots as desired. In your notebook, describe what is happening in each case. How do these results
compare to the logistic map from the midterm?

Problem 3. Examine structure of the maxima (25pt). Create a function findmaxima(x) that isolates local
maxima of a particular solution x(t) for t > T0 (to discard transient behavior). For each value c in the
range [2, 6] (with a small mesh spacing 0.001), solve Eqs. (1), use this function to isolate the set of maximal
points, then plot each maximal point (c, x) on the same scatter plot. After plotting all such points, you will
have a graph of (the multivalued function of) the asymptotic local maxima of x vs. c. The range of x in
the plot should be from [3, 12]. Comment on your findings. What happens if you replace x by y or z in this
procedure?



## Honor Pledge

I pledge that all the work in this repository is my own with only the following exceptions:

* Content of starter files supplied by the instructor;
* Code borrowed from another source, documented with correct attribution in the code and summarized here.

Signed,

Michael Seaman
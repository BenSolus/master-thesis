---
layout:    post
author:    Bennet Carstensen
title:     Let's get started - Motivation by a model problem
date:      March 30, 2017
permalink: /blog/intro
comments:  true
---

In this blog I will collect ideas for my master thesis with the working title
*"An efficient implementation of the Green cross approximation method on
GPGPUs"* based on the article
[Approximation of integral operators by Green quadrature and nested cross approximation](https://link.springer.com/article/10.1007/s00211-015-0757-y) [\[1\]]({{ site.baseurl }}/refs).

We consider integral equations e.g. of the form

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

$$
\int\limits_0^1 g(x, y)\ u(y)\ dy = f(x) \qquad \text{for all } x \in [0, 1].
$$

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

with

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

$$
g(x, y) :=
\begin{cases}
  -\log|x - y| & \text{if } x \neq y, \newline
  0            & \text{otherwise}
\end{cases}
\qquad \text{for all } x, y \in [0, 1]
$$

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

as an sample *kernel function* to construct the simple one-dimensional
[Fredholm integral equation](https://en.wikipedia.org/wiki/Fredholm_integral_equation).

<!--more-->

To compute an approximation of the solution, we consider a test space
V and look for the Galerkin approximation $$u \in V$$ satisfying the
variational equation

$$
\int\limits_0^1 v(x) \int\limits_0^1 g(x, y)\ u(y)\ dy\ dx =
\int\limits_0^1 v(x)\ f(x)\ dx \qquad \text{for all } v \in V.
$$

Since our first priority is the Galerkin discretization of the variational
formulation, we let $$n \in \mathbb{N}$$, introduce basis functions

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

$$
\phi_i(x) :=
\begin{cases}
  1 & \text{if } x \in \left [\frac{i - 1}{n}, \frac{i}{n} \right ]
  \newline
  0 & \text{otherwise}
\end{cases}
\qquad \text{for all } i \in [1 : n], x \in [0, 1]
$$

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

and replace the space $$V$$ by

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

$$
V_n := span\{\phi_i : i \in [1 : n]\}
$$

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

after the standard Galerkin approach and thus we get the following
finite-dimensional variational problem and look for an $$u_n \in V_n$$ such that

$$
\int\limits_0^1 v_n(x) \int\limits_0^1 g(x, y)\ u_n(y)\ dy\ dx =
\int\limits_0^1 v_n(x)\ f(x)\ dx \qquad \text{for all } v_n \in V_n.
$$

To compute $$u_n$$, we express it as a linear combination of our
basis with the coefficient vector $$z \in \mathbb{R}^n$$:

<!-- lint disable emphasis-marker-->

$$
u_n = \sum\limits_{j = 1}^n z_j\ \phi_j
$$

<!-- lint enable emphasis-marker-->

and thus we can describe an $$n$$-dimensional linear system of equations

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

$$
\sum\limits_{j = 1}^n z_j \int\limits_0^1 \phi_i(x) \int\limits_0^1
g(x, y)\ \phi_j(y)\ dy\ dx =
\int\limits_0^1 \phi_i(x)\ f(x)\ dx
\qquad \text{for all } i \in [1 : n].
$$

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

Defining the matrix $$G \in \mathbb{R}^{n \times n}$$ and the vector
$$b \in \mathbb{R}^n$$ by

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

$$
\begin{array}{rcl}
g_{ij} & := & \int\limits_0^1 \phi_i(x) \int\limits_0^1 g(x, y)\ \phi_j(y)\ dy\ dx, \newline
b_i & := & \int\limits_0^1 \phi_i(x)\ f(x)\ dx,
\end{array}
\qquad \text{for all } i, j \in [1 : n]
$$

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

we can write this linear system in the form

$$
Gz = b.
$$

Since $$G$$ is an $$n$$-dimensional square matrix, we have to store $$n^2$$
coefficients. To approximate $$u$$ reasonably accurate with $$u_n$$,
we have to choose $$n$$ rather large, thus, storing $$n^2$$ coefficients is
rather unattractive.

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

The matrix $$G$$ is typically non-sparse since the integral operator is
non-local, in other words we have $$g(x, y) \neq 0$$ for almost all $$x, y \in
[0, 1]$$, and thus $$g_{ij} \neq 0$$ for all $$i, j \in [1 : n]$$.

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

For those kind of matrices, we summarize common techniques proposed to handle
them into three categories: kernel-based-, matrix-based- and hybrid techniques.

The first category replaces the kernel function $$g$$ by a degenerated
approximation $$\tilde{g}$$ which the corresponding method handles efficiently.
The panel clustering method [\[2\]]({{ site.baseurl }}/refs) for example uses
Taylor expansion or interpolation [\[3\]]({{ site.baseurl }}/refs) to
approximate the kernel function.

Matrix-based techniques work directly with the matrix entries $$g_{ij}$$ to
construct an approximation. The cross-approximation
[\[4\]]({{ site.baseurl }}/refs) method in particular computes a small number of
*"crosses"*, each represented by one row and column of a submatrix of $$G$$,
which leads to a low-rank approximation. For an efficient implementation we
only partial evaluate the full matrix to retrieve the "crosses" thus we have
to specify a pivoting strategy and an error estimator to construct the
low-rank matrix which leads to the adaptive cross
approximation [\[5\]]({{ site.baseurl }}/refs).

Kernel- and matrix-based approximations both have advantages and disadvantages.
Kernel-based approximations are typically proven to converge at a certain
rate, and they do not depend on the choice of basis functions, but they are
frequently less efficient than matrix-based approximations. Matrix-based
approximations typically lead to high compression rates but finding a reliable
stopping criterion becomes challenging because of the partial evaluation and
efficient pivoting attempts [\[5\]]({{ site.baseurl }}/refs) presently rely on
further stability assumptions.

The third category tries to combine kernel and matrix-based techniques to
compensate the disadvantages while maintaining the advantages. The hybrid cross
approximation [\[6\]]({{ site.baseurl }}/refs) for example uses cross
approximation on a small submatrix resulting from interpolation
ergo avoiding error estimators.

Building up on the implementation of the newly developed Green cross
approximation (GCA) method [\[1\]]({{ site.baseurl }}/refs) for boundary
element matrices [\[7\]]({{ site.baseurl }}/refs) based on the
*H2Lib*[$$^1$$](#1), we are going to improve the performance of this method
by analyzing the algorithm and proposing techniques which divide the
upcoming tasks between CPU and GPU and take advantage of the architecture of
GPGPUs for computationally expensive operations.

--------------------------------------------------------------------------------
<b id="1">$$^1$$</b> Available at [http://www.h2lib.org](http://www.h2lib.org)

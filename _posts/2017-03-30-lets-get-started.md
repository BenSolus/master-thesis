---
layout:    post
author:    Bennet Carstensen
title:     Let’s get started - Motivation by a model problem
date:      March 30, 2017
permalink: /blog/intro
comments:  true
---

In this blog I will collect ideas for my master thesis with the working title
*"An efficient solver for the Green cross approximation method on GPUs"* based
on the article
[Approximation of integral operators by Green quadrature and nested cross approximation](https://link.springer.com/article/10.1007/s00211-015-0757-y) [\[1\]]({{ site.baseurl }}/refs).

We consider integral equations e.g. of the form

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

$$
\int\limits_{\Omega} g(x, y)\ u(y)\ dy = f(x) \qquad
\text{for all } x \in \Omega \subseteq \mathbb{R}^d.
$$

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

with $$d \in \{ 2, 3 \}$$ and

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

$$
g(x, y) :=
\begin{cases}
  -\frac{1}{2 \pi} \log{\left \Vert x - y \right \Vert_2} &
  \text{if } x \neq y, \ d = 2 \\
  \\
  \frac{1}{4 \pi \left \Vert x - y \right \Vert_2} &
  \text{if } x \neq y, \ d = 3 \\
  \\
  0 & \text{otherwise}
\end{cases}
\qquad \text{for all } x, y \in \Omega
$$

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

as an sample *kernel function* for a fundamental solution of the negative
Laplace operator $$-\Delta$$, e.g. occurring in the interpretation of
Poisson's equations in integral formulation.

<!--more-->

To compute an approximation of the solution, we consider a test space
V and look for the Galerkin approximation $$u \in V$$ satisfying the
variational equation

<!-- lint disable emphasis-marker-->

$$
\int\limits_{\Omega} v(x) \int\limits_{\Omega} g(x, y)\ u(y)\ dy\ dx =
\int\limits_{\Omega} v(x)\ f(x)\ dx \qquad \text{for all } v \in V.
$$

Since our first priority is the Galerkin discretization of the variational
formulation, we let $$n \in \mathbb{N}$$, introduce basis functions
$$(\phi_i)_{i \in \mathcal{I}}$$, where $$\mathcal{I}$$ is a suitable finite
index set, and replace the space $$V$$ by

<!-- lint enable emphasis-marker-->

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

$$
V_n := span\{\phi_i : i \in \mathcal{I} \}
$$

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

after the standard Galerkin approach. Furthermore, we get the following
finite-dimensional variational problem and are now searching for an $$u_n \in
V_n$$ such that

<!-- lint disable emphasis-marker-->

$$
\int\limits_{\Omega} v_n(x) \int\limits_{\Omega} g(x, y)\ u_n(y)\ dy\ dx =
\int\limits_{\Omega} v_n(x)\ f(x)\ dx \qquad \text{for all } v_n \in V_n.
$$

<!-- lint enable emphasis-marker-->

To compute $$u_n$$, we express it as a linear combination of our
basis with the coefficient vector $$z \in \mathbb{R}^{\mathcal{I}}$$:

<!-- lint disable emphasis-marker-->

$$
u_n = \sum\limits_{i \in \mathcal{I}} z_j\ \phi_j
$$

<!-- lint enable emphasis-marker-->

and thus we can describe an $$n$$-dimensional linear system of equations

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

<!-- lint disable emphasis-marker-->

$$
\sum\limits_{i \in \mathcal{I}} z_j
\int\limits_{\Omega} \phi_i(x) \int\limits_{\Omega} g(x, y)\ \phi_j(y)\ dy\ dx =
\int\limits_{\Omega} \phi_i(x)\ f(x)\ dx
\qquad \text{for all } i \in [1 : n].
$$

<!-- lint enable emphasis-marker-->

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

Defining the matrix $$G \in \mathbb{R}^{\mathcal{I} \times \mathcal{I}}$$ as
well as the vector $$b \in \mathbb{R}^{\mathcal{I}}$$ by

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

<!-- lint disable emphasis-marker-->

$$
\begin{array}{rcl}
g_{ij} & := & \int\limits_{\Omega}
\phi_i(x) \int\limits_{\Omega} g(x, y)\ \phi_j(y)\ dy\ dx, \newline
b_i & := & \int\limits_{\Omega} \phi_i(x)\ f(x)\ dx,
\end{array}
\qquad \text{for all } i, j \in \mathcal{I}
$$

<!-- lint enable emphasis-marker-->

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
\Omega$$, and thus $$g_{ij} \neq 0$$ for all $$i, j \in \mathcal{I}$$.

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

For this kind of matrices, we summarize common techniques proposed to handle
them into three categories: kernel-based, matrix-based and hybrid techniques.

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
only partially evaluate the full matrix to retrieve the "crosses" thus we have
to specify a pivoting strategy and an error estimator to construct the
low-rank matrix which leads to the adaptive cross
approximation [\[5\]]({{ site.baseurl }}/refs).

Kernel- and matrix-based approximations both have advantages as well as
disadvantages. Kernel-based approximations are typically proven to converge at
a certain rate, and they do not depend on the choice of basis functions, but
they are frequently less efficient than matrix-based approximations.
Matrix-based approximations typically lead to high compression rates.
Finding a reliable stopping criterion, however, becomes challenging because of
the partial evaluation, and furthermore, efficient pivoting strategies
[\[5\]]({{ site.baseurl }}/refs) presently rely on further stability
assumptions.

The third category tries to combine kernel and matrix-based techniques to
compensate the disadvantages while maintaining the advantages. The hybrid cross
approximation [\[6\]]({{ site.baseurl }}/refs) for example uses cross
approximation on a small submatrix resulting from interpolation,
ergo avoiding error estimators.

--------------------------------------------------------------------------------
<b id="1">$$^1$$</b> Available at [http://www.h2lib.org](http://www.h2lib.org)

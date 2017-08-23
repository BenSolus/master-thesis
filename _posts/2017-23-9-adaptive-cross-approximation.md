---
layout:    post
author:    Bennet Carstensen
title:     Adaptive cross approximation
date:      August 23, 2017
permalink: /blog/aca
comments:  true
---

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

The low-rank approximation we prepared in
[this]({{ site.baseurl }}/blog/gca/1) section still offers room
for improvment [\[1\]]({{ site.baseurl }}/refs): the matrices
$$A_{t, +}$$ and $$A_{t, -}$$ describe the influence of Neumann
and Dirichlet values of $$g(\cdot, y)$$ on $$\partial \omega_{t}$$,
respectively. Using the Pioncaré-Steklov operator, we can construct the Neumann
from the Dirichlet values. In other words we can gain $$A_{t, +}$$ from
$$A_{t, -}$$ and thus the rank of $$A_{t} = \left ( A_{t, +} A_{t, -} \right )
\in \mathbb{R}^{\hat{t} \times 2K}$$ is $$K$$.

Though it's possible to discard $$A_{t, +}$$ explicitly by approximating the
Poincaré-Steklov operator, this approach needs to solve an integral equation on
the boundary $$\partial \omega_{t}$$ and a high accuraccy for the sake of
preserving the exponential convergence of the quadrature approximation
[\[1\]]({{ site.baseurl }}/refs).

Because of their implactions on performance, we choose an implicit approach: we
assume that the properties mentioned above holds true, e.g. the rank of
$$A_{t}$$ is lower than the number of columns, and we can use an algebraic
procedure to approximate $$A_{t}$$ by a lower-rank matrix. The
adaptive cross approximation [\[1, 5\]]({{ site.baseurl }}/refs) is an
interesting method for this situation since it allows us to construct an
algebraic interpolation operator we can use to approximate matrix blocks
based on a small number of their entries.

<!--more-->

Assuming, our index sets $$\mathcal{I}$$ and $$\mathcal{J}$$ have the form of
$$\mathcal{I} = \{ 1, \cdots, n \}$$ and $$\mathcal{J} = \{ 1, \cdots, m \}$$,
$$n, m \in \mathbb{N}$$, we interpret a matrix
$$X \in \mathbb{R}^{\mathcal{I} \times \mathcal{J}}$$ by standard definitions
of matrices with $$X \in \mathbb{R}^{n \times m}$$. To find an adaptive cross
approximation of $$X$$, we first consider an LU factorization

$$
  X = LU
$$

with a lower triangular matrix $$L \in \mathbb{R}^{n \times l}$$
and an upper triangular matrix $$U \in \mathbb{R}^{l \times m}$$,
$$l = \min{\{n, m\}}$$.

With $$k \in \{ 1, \cdots, l \}$$, we can split the factors into

$$
  L = \begin{pmatrix} L_{kk} & \\ L_{*k} & L_{**} \end{pmatrix} \qquad
  L_{kk} \in \mathbb{R}^{k \times k} \qquad
  L_{*k} \in \mathbb{R}^{(n - k) \times k} \qquad
  L_{**} \in \mathbb{R}^{(n - k) \times (l - k)} \\
  U = \begin{pmatrix} U_{kk} & U_{k*} \\ & U_{**} \end{pmatrix} \qquad
  U_{kk} \in \mathbb{R}^{k \times k} \qquad
  U_{k*} \in \mathbb{R}^{k \times (m - k)} \qquad
  L_{**} \in \mathbb{R}^{(l - k) \times (m - k)}
$$

We construct row pivot indices
$$\{ i_{1}, \cdots, i_{k} \} \subseteq \mathcal{I}$$ and column pivot indices
$$\{ j_{1}, \cdots, j_{k} \} \subseteq \mathcal{J}$$ and apply the standard
LU factorization to the matrix $$\hat{X} \in \mathbb{R}^{k \times k}$$ with
$$\hat{x}_{\nu \mu} = x_{i_{\nu} j_{\mu}}, 1 \leq \nu, \mu \leq k$$. We
construct this pivot indices parallel to the LU factorization to avoid a zero
on the diagonal.

Starting with indices $$i_{1} \in \mathcal{I}$$ and $$j_{1} \in \mathcal{J}$$
such that $$x_{i_{1}, j_{1}} \neq 0$$, the first row and column of L and U are
given by

$$
l_{i_{1}, 1} = 1, \qquad u_{1, j_{1}} = x_{i_{1} j_{1}}, \qquad
l_{i, 1} = \frac{x_{i, j_{1}}}{u_{i, j_{1}}}, \qquad
u_{1, j} = x_{i_{1}, j} \qquad \text{for all } i \in \mathcal{I}, j \in
\mathcal{J}.
$$

We proceed with constructing the remainder

$$
  X^{(1)} = X - L_{\mathcal{I} \times \{1\}} U_{\{1\} \times \mathcal{J}}
$$

and look for its LU factorization. The resulting algorithm looks as follows:  
$$ \textbf{procedure } \text{aca}(X, \textbf{var } L, R, \epsilon) \\
\quad \tau_{0} \leftarrow \emptyset; \quad \sigma_{0} \leftarrow \emptyset;
\quad X^{(0)} \leftarrow X; \quad k \leftarrow 0 \\
\quad \textbf{while} \left\lVert X^{(k)}\right\rVert \geq \epsilon
\textbf{ do begin} \\
\qquad \text{Choose } i_{k + 1} \in \mathcal{I} \setminus \tau_{k}
\text{ and } j_{k + 1} \in \mathcal{J} \setminus \sigma_{k} \text{ with }
x^{(k)}_{i_{k + 1}, j_{k + 1}} \neq 0; \\
\qquad L|_{\mathcal{I} \times \{k + 1\}} \leftarrow
X^{(k)}|_{\mathcal{I} \times \{j_{k + 1}\}} / x^{(k)}_{i_{k + 1}, j_{k + 1}}; \\
\qquad U|_{\{k + 1\} \times \mathcal{J}} \leftarrow
X^{(k)}|_{\{i_{k + 1}\} \times \mathcal{J}}; \\
\qquad X^{(k + 1)} \leftarrow X^{(k)} - L|_{\mathcal{I} \times \{k + 1\}}
U|_{\{k + 1\} \times \mathcal{J}}; \\
\qquad \tau_{k + 1} \leftarrow \tau_{k} \cup \{i_{k + 1}\}; \quad
\sigma_{k + 1} \leftarrow \sigma_{k} \cup \{j_{k + 1}\}; \\
\qquad k \leftarrow k + 1; \\
\quad \textbf{end} \\
\textbf{end}
$$

This constructs matrices
$$L \in \mathbb{U}^{\mathcal{I} \times \{ 1, \cdots, k \}}$$ and
$$R \in \mathbb{R}^{\{1, \cdots, k\} \times \mathcal{J}}$$ such that

$$
  X = LU + X^{(k)}
$$

holds and the error

$$
  \left\lVert X^{(k)} \right\rVert = \left\lVert X - LU \right\rVert  < \epsilon
  \in \mathbb{R}
$$

is sufficiently small. We can use

$$
  \widetilde{X} := LU
$$

as a rank-k approximation of $$X$$.

We can interpret this approximation, like mentioned before, as an algebraic
procedure: we introduce the matrix
$$P^{(\gamma)} \in \mathbb{R}^{\gamma \times \mathcal{I}}, \gamma \in
\{ 1, \cdots, k \},$$ by

$$
  P^{(\gamma)}z := \begin{pmatrix} z_{i_{1}} \\ \vdots \\ z_{i_{\gamma}}
  \end{pmatrix} \quad \text{for all } z \in \mathbb{R}^{\mathcal{I}},
$$

mapping a vector to the first $$\gamma$$-th selected pivot elements. Interpeted
as an algebraic interpolation, the pivot elements play the role of interpolation
points.

<!-- For an equivalent of Lagrange polynomials, we first need to realize that we
can interprete $$L$$ and $$U$$ as an array of column and row vectors,
respectively:

$$
  \begin{align*}
    L = \left ( l^{(1)} \cdots l^{(k)} \right )
      &= \left (
    \begin{pmatrix} l_{1, j_{1}} \\ \vdots \\ \vdots \\ l_{k, j_{1}} \end{pmatrix}
    \begin{pmatrix} 0 \\ l_{2, j_{2}} \\ \vdots \\ l_{k, j_{2}} \end{pmatrix}
    \cdots
    \begin{pmatrix} 0 \\ \vdots \\ 0 \\ l_{k, j_{k}} \end{pmatrix}
    \right ) \\


    U = \begin{pmatrix} u^{(1)} \\ \vdots \\ u^{(k)} \end{pmatrix}
      &= \begin{pmatrix}
          ( & u_{i_{1}, 1} & \cdots & \cdots & u_{i_{1}, k} & ) \\
          ( & 0 & u_{i_{2}, 2} & \cdots & u_{i_{2}, k} & ) \\
          ( & 0 & \cdots & 0 & u_{i_{k}, k} & )
        \end{pmatrix}.
  \end{align*}
$$ -->

For an equivalent of Lagrange polynomials, we first need to realize that the
algorithm effectively sets all entries of the $$i$$-th row and
$$j$$-th column, $$i \in \{ i_1, \cdots, i_{\gamma} \}, j \in
\{ j_{1}, \cdots, j_{\gamma} \}$$, to zero in the $$l$$-th remainder matrix

$$
  X^{(\gamma)} = X - \tilde{X}^{(1)} - \cdots - \tilde{X}^{\gamma}.
$$

Applying the permutations on this matrix results in a block matrix

$$
  P^{(l)}X^{(l)} = \begin{pmatrix}
    0_{\gamma \times \gamma} & 0_{\gamma \times (m - \gamma)} \\
    0_{(n - \gamma) \times \gamma} & \hat{X}
  \end{pmatrix} \qquad \hat{X} \in
  \mathbb{R}^{(n - \gamma) \times (m - \gamma)}  
$$

and thus with

$$
  l_{i, \gamma} = l_{i}^{(\gamma)} = 0 \qquad \text{for all } i \in
  \{ i_{1}, \cdots, i_{\gamma-1} \},
$$

we can describe the entries of $$PL \in \mathbb{R}^{k \times k}$$ by

$$
  (PL)_{\nu \mu} = l_{i_{\nu}}^{(\mu)} \qquad \text{for all } \nu, \mu \in
  {1, \cdots, k}.
$$

This matrix is lower triangular and due to our choice of scaling, the diagonal
entries are equal to one, ergo the matrix is also invertible, resulting in

$$
  V := L(PL)^{-1} \in \mathbb{\mathcal{I} \times k}
$$

being well-defined. The columns are now acting as the Lagrange polynomials, and
the algebraic operator given by

$$
  \mathfrak{I} := VP
$$

satisfies the projection property of

$$
  \mathfrak{I}L = VPL = L(PL)^{-1}PL = L.
$$

Since the rows $$i_{1}, \cdots, i_{k}$$ of $$X^{(k)}$$ are set to zero,
we have

$$
  0 = PX^{(k)} = P(X - LU),
$$

which leads to

$$
  \mathfrak{I}X = (\mathfrak{I}L)U = LU,
$$

the low-rank approximation LU is resulting from algebraic interpolation.

---
layout:    post
author:    Bennet Carstensen
title:     \(\mathcal{H}^2\)-Matrices
date:      May 22, 2017
permalink: /blog/h2matrices
comments:  true
---

<!-- lint disable emphasis-marker-->

As we already [mentioned]({{ site.baseurl }}/blog/intro), matrices resulting
from integral equations are typically non-sparse, thus we consider
*data-sparsity*:
We are interested in an alternative representation or approximation of the
matrix requiring just a small number of coefficients. Hierarchical matrix
methods [\[8\]]({{ site.baseurl }}/refs), [\[9\]]({{ site.baseurl }}/refs)
achieve this by dividing the matrix into submatrices approximating those who
hold an *admissibility condition* or storing them in standard two-dimensional
arrays if they are small enough.

<!-- lint enable emphasis-marker-->

<!--more-->

**Definition (Tree)** [\[10\]]({{ site.baseurl }}/refs) *Let* $$N$$ *be
a finite set, let* $$r \in N$$ *and* $$E \subseteq N \times N$$*. Define*
$$\mathcal{T} := (N, r, E)$$*. We call* $$\mathcal{T}$$ *a tree if for each*
$$n \in N$$ *there is precisely one sequence* $$n_0, n_1, \cdots, n_l \in N$$
*with* $$l \in \mathbb{N_0}$$ *such that*

<!-- lint disable emphasis-marker-->

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

$$ v_0 = r, \qquad v_l = v, \qquad (v_{i-1}, v_i) \in E \qquad \textit{for all }
i \in [1:l]. $$

<!-- lint enable no-shortcut-reference-link no-undefined-references-->

<!-- lint disable emphasis-marker-->

*We call* $$N$$ *the nodes,* $$r$$ *the root and* $$E$$ *the edges of*
$$\mathcal{T}$$. *Furthermore, we call*

$$
  \text{sons}(n) := \{m \in N : (n, m) \in E \}
$$

*the set of sons of* $$n$$ *and* $$n$$ *the father of m. The father is
unique for each* $$n \in N$$ *except the root, which has no father.
Hereinafter we are using* $$\mathcal{T}$$ *as a synonym for* $$N$$*, e.g., we
write* $$n \in \mathcal{T}$$.

One might interpret a tree as a directed graph where each node represents a
subset of a index set, e.g. those resulting from the discretization of the
Galerkin approximation. To divide the set and thus the matrix into subsets we
consider a clustering method.

**Definition (Cluster tree)** *Let* $$\mathcal{I}$$ *be a finite set, let*
$$\mathcal{T} = (N, r, E)$$ *be a tree and let be a subset*
$$\hat{t} \subseteq \mathcal{I}$$ *given for each* $$t \in \mathcal{T}$$.  
*We call*
$$\mathcal{T}_{\mathcal{I}} := (N, r, E, \mathcal{I}, (\hat{t}_{t \in V}))$$ *a
cluster tree for $$\mathcal{I}$$ if the tree fullfills the follwing
requirements:*

*   *The index set corresponding to the root is $$\mathcal{I}$$, i.e.,*

$$
\hat{r} = \mathcal{I}
$$

*   *Given a non-leaf node, the index set consists of the union of the index
    sets of its sons, i.e.,*

$$
\hat{t} = \bigcup_{t' \in \text{sons}(t)} \hat{t}' \qquad \textit{for all } t
\in \mathcal{T}_{\mathcal{I}}, \text{sons}(t) \neq \emptyset.
$$

*   *Given two different sons of a node, their index sets are disjoint, i.e.,*

$$
\hat{t} \cap \hat{t}_2 \neq \emptyset \Rightarrow t_1 = t_2 \qquad
\textit{for all } t \in \mathcal{T}_{\mathcal{I}}, t_1, t_2 \in \text{sons}(t).
$$

This way, we can represent submatrices of a matrix
$$G \in \mathbb{R}^{\mathcal{I} \times \mathcal{J}}$$, with index sets
$$\mathcal{I}$$ and $$\mathcal{J}$$, by pairs of cluster chosen
from two cluster trees $$\mathcal{T}_{\mathcal{I}}$$ and
$$\mathcal{T}_{\mathcal{J}}$$. To find suitable submatrices to approximate, we
introduce another tree structure.

**Definition (Block tree)** *Let* $$\mathcal{T}_{\mathcal{I}}$$ *and*
$$\mathcal{T}_{\mathcal{J}}$$ *be cluster trees for index sets* $$\mathcal{I}$$
and $$\mathcal{J}$$*. We call* $$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}
= (N, r, E)$$ *a block tree for* $$\mathcal{T}_{\mathcal{I}}$$ *and*
$$\mathcal{T}_{\mathcal{J}}$$ *if the tree fullfills the follwing
requirements:*

*   *The nodes are pairs of clusters, i.e.,*

$$
V \subseteq \mathcal{T}_{\mathcal{I}} \times \mathcal{T}_{\mathcal{J}}.
$$

*   *The root consists of the roots of $$\mathcal{T}_{\mathcal{I}}$$ and
    $$\mathcal{T}_{\mathcal{J}}$$, e.g.,*

$$
r = (\text{root}(\mathcal{T}_{\mathcal{I}}),
    \text{root}(\mathcal{T}_{\mathcal{J}})).
$$

*   *If* $$b = (t, s) \in \mathcal{T}$$ *has sons, one can form them by the
    sons of* $$t$$ *and* $$s$$*, i.e.,*

$$
\text{sons}(b) =
\begin{cases}
  \text{sons}(t) \times \text{sons}(s) &
  if \text{ sons}(t) \neq \emptyset, \text{sons}(s) \neq \emptyset \\
  \{t\} \times \text{sons}(s) &
  if \text{ sons}(t) = \emptyset, \text{sons}(s) \neq \emptyset \\
  \text{sons}(t) \times s &
  if \text{ sons}(t) = \emptyset, \text{sons}(s) = \emptyset \\
\end{cases}
$$

*We call the nodes of a block tree*
$$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$
*for cluster trees* $$\mathcal{T}_{\mathcal{I}}$$ *and*
$$\mathcal{T}_{\mathcal{J}}$$ *blocks and denote the set of leaves by*
$$\mathcal{L}_{\mathcal{I} \times \mathcal{J}}$$.

A block tree $$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$ represents a
cluster tree for the Cartesian product index set
$$\mathcal{I} \times \mathcal{J}$$.

For any cluster tree, the labels of the leaves form a disjoint
partition of the index set, as one can show by induction. The leaves
of a block tree $$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$ in particular
form the disjoint partition

$$
\{\hat{t} \times \hat{s} : b = (t, s) \in
\mathcal{L}_{\mathcal{I} \times \mathcal{J}}\}
$$

of the index set $$\mathcal{I} \times \mathcal{J}$$, i.e, resulting from the
discretization of the Galerkin approximation.

Since we won't be able to approximate each block, we define an admissiblity
condition

$$
\text{adm}: \mathcal{T}_{\mathcal{I}} \times \mathcal{T}_{\mathcal{J}}
\rightarrow \{\text{true}, \text{false}\}
$$

which indicates whether we approximate a submatrix
$$G|_{\hat{t} \times \hat{s}}$$. We call a block $$(t, s)$$ admissible if
adm$$(t, s) = true holds.

**Definition (Admissibility condition)** *Let* $$\mathcal{T}_{\mathcal{I}}$$
*and* $$\mathcal{T}_{\mathcal{J}}$$ *be cluster trees for index sets*
$$\mathcal{I}$$ and $$\mathcal{J}$$. *We call a mapping*

$$
\text{adm}: \mathcal{T}_{\mathcal{I}} \times \mathcal{T}_{\mathcal{J}}
\rightarrow \{\text{true}, \text{false}\}
$$

*an* admissibility condition *for* $$\mathcal{T}_{\mathcal{I}}$$ *and*
$$\mathcal{T}_{\mathcal{J}}$$ *if*

$$
\begin{align*}
  \text{adm}(t, s) \Rightarrow \text{adm}(t', s) \qquad \textit{for all } t \in
  \mathcal{T}_{\mathcal{I}}, \ s \in \mathcal{T}_{\mathcal{J}}, \ t' \in
  \text{sons} (t), \\
  \text{adm}(t, s) \Rightarrow \text{adm}(t, s') \qquad \textit{for all } t \in
  \mathcal{T}_{\mathcal{I}}, \ s \in \mathcal{T}_{\mathcal{J}}, \ s' \in
  \text{sons} (s). \\
\end{align*}
$$

Since being a disjoint partition of the index set, the leaves of a block tree
$$\mathcal{T}_{\mathcal{I}} \times \mathcal{T}_{\mathcal{J}}$$ give rise to a
decomposition of a matrix $$G$$ into submatrices: We expect the leaves to
correspond to those submatrices we approximate, i.e. the ones who are
admissible, or those we represent as an array of coefficients since they are
small enough.

**Definition (Admissible and inadmissible leaves)** *Let*
$$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$ *be a block tree
for cluster trees* $$\mathcal{T}_{\mathcal{I}}$$ *and*
$$\mathcal{T}_{\mathcal{J}}$$*and let* adm *be an admissiblity condition. We
split the set of leaves*
$$\mathcal{L}_{\mathcal{I} \times \mathcal{J}}$$ *of*
$$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$ *into the sets*

$$
  \mathcal{L}^{+}_{\mathcal{I} \times \mathcal{J}} := \{(t, s)
  \in \mathcal{L}_{\mathcal{I} \times \mathcal{J}} \ : \ \text{adm}(t, s) =
  \text{true}\},
$$

*of* admissible *and*

$$
  \mathcal{L}^{-}_{\mathcal{I} \times \mathcal{J}} := \{(t, s)
  \in \mathcal{L}_{\mathcal{I} \times \mathcal{J}} \ : \ \text{adm}(t, s) =
  \text{false} \} = \mathcal{L}_{\mathcal{I} \times \mathcal{J}} \setminus
  \mathcal{L}^{+}_{\mathcal{I} \times \mathcal{J}},
$$

*inadmissible leaves.*

**Definition (Admissible block tree)** *Let*
$$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$ *be a block
tree, let* $$\mathcal{L}^{+}_{\mathcal{I} \times \mathcal{J}} := \{(t, s)$$
*and* $$\mathcal{L}^{-}_{\mathcal{I} \times \mathcal{J}}$$ *denote the
admissible and inadmissible leaves of
$$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$.

*We call a block tree*
$$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$ admissible *if*

$$
(t, s) \in \mathcal{L}^{-}_{\mathcal{I} \times \mathcal{J}} \Rightarrow
(t \in \mathcal{L}_{\mathcal{I}} \lor s \in
\mathcal{L}_{\mathcal{J}}) \qquad \textit{holds for all }
(t, s) \in \mathcal{T}_{\mathcal{I} \times \mathcal{J}},
$$

*and we call it* strictly admissible if

$$
(t, s) \in \mathcal{L}^{-}_{\mathcal{I} \times \mathcal{J}} \Rightarrow
(t \in \mathcal{L}_{\mathcal{I}} \land s \in
\mathcal{L}_{\mathcal{J}}) \qquad \textit{holds for all }
(t, s) \in \mathcal{T}_{\mathcal{I} \times \mathcal{J}}.
$$

With cluster trees $$\mathcal{T}_{\mathcal{I}}$$, $$\mathcal{T}_{\mathcal{J}}$$
and an admissibility condition given, we can construct a (strictly) admissible
block tree by checking the admissibility condition for the root
$$(r_{\mathcal{I}}, r_{\mathcal{I}})$$. We finished our construction if the
condition holds true. Otherwise, we check the sons and further descendants.

Constructing an approximation of $$G$$ now means defining approximations for
all $$G_{\hat{t} \times \hat{s}}$$ with $$b = (t, s) \in
\mathcal{L}^{+}_{\mathcal{I} \times \mathcal{J}}$$.

**Definition (Hierarchical matrix)** *Let*
$$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$ *be a block tree for cluster
trees* $$\mathcal{T}_{\mathcal{I}}$$, $$\mathcal{T}_{\mathcal{J}}$$*, and let*
$$\mathcal{L}^{+}_{\mathcal{I} \times \mathcal{J}}$$ *denote the admissible
leaves of* $$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$.

*We call a matrix* $$G \in \mathbb{R}^{\mathcal{I} \times \mathcal{J}}$$ *a
hierarchical matrix* (or $$\mathcal{H}$$-matrix) *of local rank*
$$k \in \mathbb{N}$$ of $$\mathcal{T}_{\mathcal{I} \times \mathcal{J}}$$*, if
for each* $$b = (t, s) \in \mathcal{L}^{+}_{\mathcal{I} \times \mathcal{J}}$$
*we find* $$A_b \in \mathbb{R}^{\hat{t} \times k}$$ *and*
$$B_b \in \mathbb{R}^{\hat{s} \times k}$$ *such that*

$$
G|_{\hat{t} \times \hat{s}} = A_b B_b^*,
$$

<!-- lint disable no-inline-padding-->

*where* $$B_b^* \in \mathbb{R}^{k \times \hat{s}}$$ *denotes the transposed of*
$$B_b$$.

<!-- lint enable no-inline-padding-->

Approximating admissible submatrices leads in typical applications to reduced
storage requierements of a hierarchical matrix to $$\mathcal{O}(nk \log{(n)})$$,
where $$n := \max{\{ \# \mathcal{I}, \# \mathcal{J} \}}$$.

We can avoid the logarithmic factor by refining the representation: we choose
sets of basis vectors for all clusters $$t \in \mathcal{T}_{\mathcal{I}}$$,
$$s \in \mathcal{T}_{\mathcal{J}}$$ and represent the admissible blocks by
these basis vectors.

**Definition (Cluster index)** *We call a family*
$$K := \left ( K_t \right )_{t \in \mathcal{T}_{\mathcal{I}}}$$ *of an finite
index set* $$\mathcal{I}$$ *a cluster index for* $$\mathcal{T}_{\mathcal{I}}$$.

**Definition (Cluster basis)** *Let*
$$K = \left ( K_t \right )_{t \in \mathcal{T}_{\mathcal{I}}}$$ *be a cluster
index for* $$\mathcal{T}_{\mathcal{I}}$$*. We call a family*
$$V = \left ( V_t \right )_{t \in \mathcal{T}_{\mathcal{I}}}$$ *of matrices
a cluster basis of rank* $$K$$ *if*

*   $$V_t \in \mathbb{R}^{\hat{t} \times K_t}$$ *for all*
    $$t \in \mathcal{T}_{\mathcal{I}}$$ *holds and*
*   *there is a matrix* $$E_{t} \in \mathbb{R}^{K_{t'} \times K_{t}}$$
    *for all* $$t \in \mathcal{T}_{\mathcal{I}}$$ *and all*
    $$t' \in \text{sons}(t)$$ *satisfying*

$$
V_t|_{\hat{t}' \times K_t} = V_{t'}E_{t'}. \qquad (*)
$$

*In this case, we call the matrices $$V_t$$ basis matrices and the matrices
$$E_t$$ transfer matrices.*

**Definition (Nested representation)** *Let*
$$V = \left ( V_t \right )_{t \in \mathcal{T}_{\mathcal{I}}}$$ *be a cluster
basis of rank* $$K$$ *and let*
$$\left ( E_{t} \right )_{t \in \mathcal{T}_{\mathcal{I}}}$$ *be a family of
transfer matrices satisfying* $$(*)$$.

*We call $$\left ( \left ( V_t \right )_{t \in \mathcal{L}_{\mathcal{I}}},  
\left ( E_{t} \right )_{t \in \mathcal{T}_{\mathcal{I}}} \right )$$ a nested
representation of the cluster basis $$V$$.

In a nested representation, we give the matrices
$$V_t$$ of $$t \in \mathcal{L}_{\mathcal{I}}$$ explicitly and otherwise
implicitly by $$(*)$$.

**Definition ($$\mathcal{H}^2$$-matrix)** *Let*
$$G \in \mathbb{R}^{\mathcal{I} \times \mathcal{J}}$$ *be a hierarchical matrix
for nested cluster bases*
$$\left ( V_t \right )_{t \in \mathcal{T}_{\mathcal{I}}}$$ *of rank*
$$K = \left ( K_t \right )_{t \in \mathcal{T}_{\mathcal{I}}}$$ *and*
$$\left ( V_s \right )_{s \in \mathcal{T}_{\mathcal{J}}}$$ *of rank*
$$L = \left ( L_s \right )_{s \in \mathcal{T}_{\mathcal{J}}}$$*, if for each*
$$b = (t, s) \in \mathcal{L}^{+}_{\mathcal{I} \times \mathcal{J}}$$ *we can
find* $$S_b \in \mathbb{R}^{K_t \times L_s}$$ *such that*

$$
G|_{\hat{t} \times \hat{s}} = V_t S_b W_s^*,
$$

*we call* $$G$$ *a* $$\mathcal{H}^2$$-matrix*. We call the matrices* $$S_b$$
coupling matrices.

Representing all admissible submatrices and the cluster basis by the coupling
and transfer matrices respectively, we can reduces the storage requirements of
an $$\mathcal{H}^2$$-matrix* to $$\mathcal{O}(n k)$$, where
$$k := \max{\{ \# K_t, \# L_s : t \in \mathcal{T}_{\mathcal{I}}, s \in \mathcal{T}_{\mathcal{J}}\} }$$.

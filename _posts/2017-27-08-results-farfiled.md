---
layout:    post
author:    Bennet Carstensen
title:     First results - Substituting the farfield
date:      August 27, 2017
permalink: /blog/results/farfield
comments:  true
---

<!-- lint disable no-shortcut-reference-link
                  no-shortcut-reference-image
                  no-undefined-references-->

After some time now, I have the first results to present. First the bad news: We
are (currently) don't have a faster solver compared to previous methods, which
is rather untypical for utilizing GPUs. On the bright side, we can save 20-26%
of memory used by a $$\mathcal{H}^2$$-matrices at the expense of 20-55% of
increased computation time in the $$\mathcal{H}^2$$-matrix vector
multiplication, decreasing antiproportional with the problem size. This in term
is heavily utilized in, e.g., Krylov subspace methods to solve system of
equations like the [discretized integrals]({{ site.baseurl }}/blog/intro) we
consider.

<!--more-->

First lets construct $$\mathcal{H}^2$$-matrix. We construct the matrix based
on the unit sphere as an example surface consisting of $$n$$ triangles. $$n$$
will also describe the dimension of the matrix. We set a maximum leaf size of
64 in our cluster and a maximum error of the adaptive cross approximation of
$$10^{-4}$$. Furthermore, we apply quadrature of the order of 2 to regular
integrals, 4 to singular integrals and 8 to approximate Green's second identity,
resulting in the following memory requirements for $$\mathcal{H}^2$$-matrix of
differrent dimensions:

{:.center}
![h2-memory]({{ site.baseurl }}/assets/graphs/h2-memory.png)

{:.center}
**Figure 1**: Memory requirements of $$\mathcal{H}^2$$-matrices and their
proportion of neafield and farfield coupling matrices.

For our current algorithm, we concentrate on the farfield part of
$$\mathcal{H}^2$$-matrix which occupies the 20-26% of memory we want save.
The details of our algorithm will follow. Regardless, we are basically computing
the all farfield coupling matrices on-the-fly during the matrix vector
multiplication phase instead of conventionaly saving the matrices and
multiplying them with the input vector.

{:.center}
![farfield-memory]({{ site.baseurl }}/assets/graphs/farfield-memory.png)

{:.center}
**Figure 2**: Memory requirements for saving the coupling matrices compared to
              the informations needed to compute those matrices.

The information needed to compute the coupling matrices and to perform the
follwing matrix vector multiplication, including the
*   geometry,
*   quadrature points and weights and
*   bookkeeping informations of the
    *   dimensions of each farfield matrix and
    *   their row and column offset relative to the whole
        $$\mathcal{H}^2$$-matrix,

occupies only a fraction of the memory used to save the matrices directly.

While this a lot of memory one could potentially save, for 3D problems, the
tasks for computing those matrices involves calculating 4D integrals for each
entrie of each farfield matrix.

Although this is an utterly computational heavy undertaking, evaluating
quadrature rules can be efficiently vectorized thus being well suited for GPUs.
We performed the following tests on a
[Intel Xeon E5-4640](https://ark.intel.com/de/products/64603/Intel-Xeon-Processor-E5-4640-20M-Cache-2_40-GHz-8_00-GTs-Intel-QPI)
CPU and a [NVIDIA Tesla K20X](http://international.download.nvidia.com/tesla/pdf/tesla-k20x-board-spec.pdf)
GPU accelerator.

{:.center}
![farfield-performance]({{ site.baseurl }}/assets/graphs/farfield-performance.png)

{:.center}
**Figure 2**: Time needed to perform the standard $$\mathcal{H}^2$$-matrix
              vector multiplication compared to the same task while computing
              the farfield matrices and their part of the product on the GPU.

As already mentioned, computing the coupling matrices on-the-fly, even when done
on the GPU, performes slightly worse than solely performing the matrix vector
multiplication. Still, saving the relative amount of memory observed here, this
algorithm might be a real consideration where memory is sparse compared to
computing power. One might also consider distributing the computations accross
multiple GPUs thus speeding up the performance.

Furthermore, consulting Figure 1, another 70-75% of the memory required for a
$$\mathcal{H}^2$$-matrix occupies the neafield. Since this represents the
singular integrals, additional work needs to be done to compute it on-the-fly
and integrate it in the current algorithm while keeping work done by threads
and GPUs balanced. The result would be saving almost 90% of the global
memory requirements of $$\mathcal{H}^2$$-matrices.

---
layout:    post
author:    Bennet Carstensen
title:     \(\mathcal{H}^2\)-Matrix-vector multiplication
date:      May 22, 2017
permalink: /blog/h2mvm
comments:  true
---

<!-- lint disable no-shortcut-reference-link no-undefined-references-->

Alongside the linear complexity of the storage requirements is the time
complexity of the $$\mathcal{H}^2$$-Matirx-vector multiplication strictly linear
to the problem size [\[9\]]({{ site.baseurl }}/refs). Yet, because
of the nested representation of a cluster basis, for a given admissible leaf
$$b = (t, s) \in \mathcal{L}^{+}_{\mathcal{I} \times \mathcal{J}}$$, evaluating

$$
  G|_{\hat{t} \times \hat{s}} = V_{t} S_{b} W_{s}^{*}
$$

requires special algorithms for computing $$W_{s}^{*}x$$ and $$V_{t} y_{t}$$
for given vectors $$x \in \mathbb{F}^{\hat{s}}$$ and
$$y_{t} \in \mathbb{F}^{K_{t}}$$.

<!--more-->

Let $$W = (W_{s})_{s \in \mathcal{T}_{\mathcal{J}}}$$ be a cluster basis of rank
$$L = (L_{s})_{s \in \mathcal{T}_{\mathcal{J}}}$$ and let
$$F = (F_{s})_{ \mathcal{T}_{\mathcal{J}}}$$ be a family of corresponding transfer
matrices.

For $$s \in  \mathcal{T}_{\mathcal{J}}$$ with sons$$(s) = \emptyset$$  we can
access $$W_{s}$$ in the nested representation of $$W$$ and compute

$$
x_{s} := W_{s}x
$$

directly for all $$x \in \mathbb{F}^{\hat{s}}$$.

If sons$$(s) \neq \emptyset$$, we define $$n := \#\text{sons}(s)$$

---
layout:    post
author:    Bennet Carstensen
title:     Let's get started - Constructing a basic model problem
date:      March 10, 2017
permalink: /blog/intro
comments:  true
---

In this blog I will collect ideas for my master thesis with the working title
*"An efficient implementation of the Green cross approximation method on
GPGPUs"* based on the article
[Approximation of integral operators by Green quadrature and nested cross approximation](https://link.springer.com/article/10.1007/s00211-015-0757-y) [\[1\]]({{ site.baseurl }}/refs).

We will start with considering a model problem to motivate the project and to
show some results later on. Some examples could be solving some kind of problem
for the [Laplace equation](https://en.wikipedia.org/wiki/Laplace%27s_equation).
Since I am not quite sure which problem might be suitable we start of with the
most basic model problem from which a more complex problem can inherit from.

<!--more-->

The most fundamental model problem might be just representing some points in a
multidimensional space. This can implemented in C with a simple struct:

```c
   1: /** @brief Representation of objects in d-dimensional space. */
   2: typedef struct _objects objects;
   3:
   4: /** @brief Pointer to @ref objects. */
   5: typedef objects* pobjects;
   6:
   7: /** @brief Pointer to constant @ref objects.*/
   8: typedef const objects* pcobjects;
   9
  10: struct _objects {
  11:   /** @brief Location of the objects. */
  12:   real **x;
  13:
  14:   /** @brief Number of objects. */
  15:   uint n;
  16:
  17:   /** @brief Dimension of the space containing the objects. */
  18:   uint dim;
  19:
  20:   /** @brief Internal field containing all informations about the objects */
  21:   real *mem;
  22: };
```

This struct represents `n` objects in a `dim`-dimensional space. Here, `x` is
an array of length `n`, each entry pointing to a disjunct section of `mem` of
length `dim` representing the coordinates of an object.

This objects can easily be translated into an
[cluster geometry](http://www.h2lib.org/doc/d4/dbd/group__clustergeometry.html)
object provided by the [H2Lib](http://www.h2lib.org):

```c
   1: pclustergeometry
   2: build_from_objects_clustergeometry(pobjects objs)
   3: {
   4:   const uint n   = objs->n;
   5:   const uint dim = objs->dim;
   6:
   7:   pclustergeometry cg = new_clustergeometry(dim, n);
   8:
   9:   for(uint i = 0; i < n; ++i)
  10:     for(uint j = 0; j < dim; ++j)
  11:       cg->x[i][j] = cg->smin[i][j] = cg->smax[i][j] = objs->x[i][j];
  12:
  13:   return cg;
  14: }
```

With the help of an index set representing the index of each object

```c
   1: uint *idx = (uint *) allocmem((size_t) objs->n * sizeof(uint));
   2:
   3: for(uint i = 0; i < objs->n; ++i)
   4:   idx[i] = i;
```
and the cluster geometry object we can construct a
[cluster](http://www.h2lib.org/doc/d7/de3/group__cluster.html) tree of our
objects with `build_adaptive_cluster()`. Each cluster, represented by a
bounding box, contains a number of objects whose coordinates are inside the
bounding box.

A first test is to check if for each cluster each object is actually inside the
bounding box of the cluster. This way we can ensure that this step of the
construction of a
[$$\mathcal{H}^2$$-matrix](http://www.h2lib.org/doc/d7/ddd/group__h2matrix.html)
later on won't cause any trouble. Furthermore, we can continue from this point
by furhter specifing the model problem and by continuing the construction of a
$$\mathcal{H}^2$$-matrix.

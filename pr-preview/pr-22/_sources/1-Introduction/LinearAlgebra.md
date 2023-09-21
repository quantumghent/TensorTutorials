---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Julia
  language: julia
  name: julia-1.9
---

# (Multi-) Linear Algebra

```{contents} Contents
:depth: 2
```

## Overview

This lecture covers some basic linear algebra concepts and operations, that will serve as
foundation for most of what follows. The goal is to provide some intuitive understanding of
the concepts, without insisting on too much mathematical rigour. The most important goal is
to introduce and define the concept of a tensor, without resorting to the usual mathematical
definition, which is not very intuitive.

Simultaneously, the lecture also showcases some of the features of
[TensorKit.jl](https://github.com/Jutho/TensorKit.jl), a Julia package that is extremely
well-suited for the demonstration of the concepts that are discussed.

```{code-cell} julia
using TensorKit
```

Before discussing tensor networks, it is necessary to understand what tensors are.
Furthermore, before really understanding tensors, it is instructive to reiterate some basic
concepts of linear algebra for the case of vectors and matrices, which are nothing but
specific cases of tensors. In fact, many of the concepts and ideas that are introduced and
discussed are defined in terms of thinking of tensors as vectors or matrices.

In what follows, vectors and matrices will be thought of from the viewpoint of computers,
where they are represented using regular one- and two-dimensional arrays of either real or
complex numbers. Nevertheless, much of the discussion can be readily generalized to
arbitrary vector spaces and linear maps.

## Vectors and Matrices

In general, a vector is an object in a vector space, which can be described by a list of
numbers that correspond to the components of the vector in some basis. For example, a vector
in a two-dimensional space is in its most general form described by
$\vec{v} = \left[v_1, v_2\right]^T$.

As a reminder, the defining properties of vector spaces make sure that the following
operations are well-defined:

* Vectors can be added together, i.e. $\vec{v} + \vec{w}$ is a vector.
* Vectors can be multiplied by scalars, i.e. $\alpha \vec{v}$ is a vector.
* These operations behave as expected, i.e. there is a notion of associativity, commutativity, and distributivity.

Given two such vector spaces (not necessarily distinct) it is possible to define a linear
map between them, which is just a function that preserves the vector space structure. In
other words, a linear map $A \colon V \rightarrow W$ maps vectors from one vector space $V$
to another vector space $W$. Because of the structure of vector spaces, and the requirement
of linearity, such a map is completely determined by its action on the basis vectors of $V$.
This leads in a very natural way to the notion of a matrix by considering the following
construction, where $v_i$ are the components of $\vec{v}$ and $w_i$ are the components of
$\vec{w}$:

```{math}
:label: eq:linear_map
\begin{array}{rcl}
A & : & V \rightarrow W\\
  &   & \vec{v} ↦ A(\vec{v}) \equiv \sum_j A_{ij} v_j = w_i \equiv \vec{w}
\end{array}
```

where $A_{ij}$ are the components of the matrix $A$ in these bases. In other words, the
abstract notion of a linear map between vector spaces can be represented by a concrete
matrix, and the action of the map is the usual matrix product.

In particular, it is instructive to think of the columns of the matrix $A$ as labelling the
components of the input vector space, also called _domain_, while the rows label the
component of the output vector space, or _codomain_.

In the context of Julia, we can create vector spaces, vectors and matrices through a syntax
that follows this very closely:

```{code-cell} julia
V = ℂ^2             # type as \bbC<TAB> 
W = ComplexSpace(3) # equivalent to ℂ^3

A = TensorMap(rand, Float64, W ← V) # ← as \leftarrow<TAB>
v = Tensor(rand, Float64, V)
w = A * v

w[1] ≈ A[1,1] * v[1] + A[1,2] * v[2]
```

````{note}
For linear maps, both notations $V \rightarrow W$ and $W \leftarrow V$ are used to denote
their codomain and domain. The choice of notation is mostly a matter of taste, as left to
right might seem more conventional for a language that reads from left to right, while right
to left is more natural when considering the mathematical usage, where matrices typically
act on vectors from left to right. In TensorKit, both notations are supported through the
`→` and `←` operators, and a Unicode-less version is also available, which defaults to `←`.
Thus, the following are all equivalent:

```{code-block} julia
A = TensorMap(rand, Float64, V → W)
A = TensorMap(rand, Float64, W ← V)
A = TensorMap(rand, Float64, W, V)
```
````

## Tensors and Tensor Products

Using the same logic as above, it is possible to generalize the notion of a linear map by
making use of the [tensor product](https://en.wikipedia.org/wiki/Tensor_product), which is
nothing but an operation that can combine two vector spaces $V$ and $W$ into a new vector
space $V \otimes W$. The tensor product is defined in such a way that the combination of
vectors from the original vector spaces preserves a natural notion of linearity, i.e. the
following equality holds for all vectors $v \in V$, $w \in W$, and scalars $\lambda$:

```{math}
:label: eq:tensor_product
(\lambda v) \otimes w = v \otimes (\lambda w) = \lambda (v \otimes w)
```

```{code-cell} julia
λ = rand()
(λ * v) ⊗ w ≈ v ⊗ (λ * w) ≈ λ * (v ⊗ w)
```

This new vector space can be equipped with a canonical basis, which is constructed by taking
the tensor product of the basis vectors of the original vector spaces. For example, if $V$
and $W$ are two-dimensional vector spaces with basis vectors $v_i$ and $w_j$, respectively,
then the basis vectors of $V \otimes W$ are given by $v_i \otimes w_j$. In other words, the
vectors in $V \otimes W$ are linear combinations of all combinations of the basis vectors of
$V$ and $W$.

When considering how to represent a vector in this new vector space, it can be written as a
list of numbers that correspond to the components of the vector in that basis. For example,
a vector in $V \otimes W$ is described by:

```{math}
:label: eq:tensor_basis
t = \sum_{i_1,i_2} t_{i_1i_2} (v_{i_1} \otimes w_{i_2})
```

```{code-cell} julia
t = Tensor(rand, Float64, V ⊗ W)
t[] # shorthand for extracting the multi-dimensional array of components
```

Here, the tentative name $t$ was used to denote that this is in fact a tensor, where
$t_{i_1i_2}$ are the components of that tensor $t$ in the basis $v_{i_1} \otimes w_{i_2}$.
Because of the induced structure of the tensor product, it is more natural and very common
to express this object not just as a list of numbers, but by reshaping that list into a
matrix. In this case, the components of the $i_1$-th row correspond to basis vectors that
are built from $v_{i_1}$, and similarly the $i_2$-th column corresponds to basis vectors
that are built from $w_{i_2}$.

As the tensor product can be generalized to more than two vector spaces, this finally leads
to the general definition of a tensor as an element of the vector space that is built up
from the tensor product of an arbitrary number of vector spaces. Additionally, the
components of these objects are then naturally laid out in a multi-dimensional array, which
is then by a slight misuse of terminology also called a tensor.

```{note}
The reshaping operation of components from a list of numbers into a multi-dimensional array
is a mapping between linear indices $I$ and Cartesian indices $i_1, i_2, \cdots,
i_N$. This is a very common and useful trick which allows reinterpreting tensors as vectors,
or vice versa.
```

```{code-cell} julia
LinearIndices((1:2, 1:3))
```

```{code-cell} julia
collect(CartesianIndices((1:2, 1:3))) # collect to force printing
```

## Tensors and Multi-Linear Maps

Due to the fact that the tensor product of vector spaces is a vector space in of itself, it
is again possible to define linear maps between such vector spaces. Keeping in mind the
definition of a linear map from {eq}`eq:linear_map`, the columns now label components of the
input vector space, while the rows label components of the output vector space. Now however,
the components of the input and output vector spaces are themselves comprised of a
combination of basis vectors from the original vector spaces. If a linear order of these
combinations can be established, the linear map can again be represented by a matrix:

```{math}
:label: eq:multilinear_map
\begin{array}{rcl}
A & : & W_1 \otimes W_2 \otimes \cdots \otimes W_M \leftarrow 
        V_1 \otimes V_2 \otimes \cdots \otimes V_N \\
  &   & v_1 \otimes v_2 \otimes \cdots \otimes v_N \mapsto 
        A(v_1 \otimes v_2 \otimes \cdots \otimes v_N) \\
  &   & = \sum_{j_1, j_2, \cdots, j_N} A_{i_1, i_2, \cdots, i_M; j_1, j_2, \cdots, j_N}
          v_{1, j_1} \otimes v_{2, j_2} \otimes \cdots \otimes v_{N, j_N} \\
  &   & = \sum_{J} A_{I;J} v_J \\
  &   & = w_1 \otimes w_2 \otimes \cdots \otimes w_M \\
\end{array}
```

```{code-cell} julia
V1 = ℂ^2
V2 = ℂ^2
W1 = ℂ^2
W2 = ℂ^2

A = TensorMap(rand, Float64, W1 ⊗ W2 ← V1 ⊗ V2)
v = Tensor(rand, Float64, V1 ⊗ V2)
w = A * v
w[] ≈ reshape(reshape(A[], 4, 4) * reshape(v[], 4), 2, 2)
```

The attentive reader might have already noted that the definition of a linear map as a
matrix strongly resembles the definition of a vector in a tensor product vector space. This
is not a coincidence, and in fact the two can easily be identified by considering the
following identification (isomorphism):

```{math}
:label: eq:tensor_isomorphism
(W \leftarrow V) \cong (W \otimes V^*) 
```

```{code-cell} julia
A = TensorMap(rand, Float64, W ← V)
B = Tensor(rand, Float64, W ⊗ V')
space(A, 2) == space(B, 2)
```

```{note}
For finite-dimensional real or complex vector spaces without additional structure, this
isomorphism is *trivial* and is just the reshaping operation of the components of a vector
into a matrix. However, note that this is a choice, which is not unique, and already differs
for
[row- and column-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order). In
a more general setting, the identification between $V \otimes W^*$ and $V \leftarrow W$ is
not an equivalence but an isomorphism. This means that it is still possible to relate one
object to the other, but the operation is not necessarily trivial.
```

## Conclusion

The entire discussion can be summarized and leads to the following equivalent definitions of
a tensor:

* A tensor is an element of a tensor product of vector spaces, which can be represented as a multi-dimensional array of numbers that indicate the components along the constituent basis vectors. Thus, a tensor is _vector-like_.
* A tensor is a multi-linear map between vector spaces, which can be represented as a matrix that represents the action of the map on the basis vectors of the input vector space. Thus, a tensor is _matrix-like_.

The equivalence of these two definitions leads to the lifting of many important facets of
linear algebra to the multi-linear setting.

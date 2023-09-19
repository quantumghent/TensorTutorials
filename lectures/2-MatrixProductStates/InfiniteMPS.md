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

# Infinite Matrix Product States

This section discusses matrix product states (MPS) in the thermodynamic limit and their
properties. Our discussion is mostly based on the excellent review
{cite}`vanderstraeten2019tangentspace`, which provides a thorough technical overview of
tangent-space methods for uniform MPS. The formal exposition is supplemented with some very
basic code examples on working with infinite MPS using
[MPSKit.jl](https://github.com/maartenvd/MPSKit.jl) at the end of this section. For more
details on the numerical implementation of routines for uniform MPS we refer to the Julia
version of the [tutorials on uniform MPS](https://github.com/leburgel/uniformMpsTutorial),
which is again based on {cite}`vanderstraeten2019tangentspace`.

```{contents} Contents
:depth: 3
```

## Matrix Product States in the Thermodynamic Limit

### Representation

The finite MPS representation introduced in the previous previous section can be readily
extended to the thermodynamic limit. Consider an infinite one-dimensional chain with a local
physical Hilbert space $\mathbb{C}^d$ of dimension $d$ at every site in the chain. An
infinite matrix product state representing the quantum state of such a system has the form

```{math}
:label: inf_mps_formula

\left | \Psi({A}) \right \rangle = \sum_{\{s\}} \boldsymbol{v}_L^\dagger \left[ \prod_{m\in\mathbb{Z}} A[m]^{s_m} \right] \boldsymbol{v}_R \left | \{s\} \right \rangle.
```

<!-- TODO: uniformize notation, conventions and discussion with finite MPS section -->

Here, each $A[m]^s$ is a $D \times D$ matrix for every value of the physical index $s$. As
before, we can alternatively view $A[m]$ as a tensor with three indices of dimensions $D
\times d \times D$, where $D$ is the the so-called *bond dimension*, which we assume to be
same at every site. As in the case of finite systems, this bond dimension controls the
entanglement of the corresponding state, where in the limit $D\to\infty$ the MPS
{eq}`inf_mps_formula` can represent any quantum state up to arbitrary accuracy. For certain
classes of quantum states such low-energy states of gapped systems however we expect that we
can accurately approximate these states using a much smaller bond dimension. Note that while
in Eq. {eq}`inf_mps_formula` we have also introduced two boundary vectors $\boldsymbol{v}_L$
and $\boldsymbol{v}_R$, but, as we work on an infinite system, the boundary conditions will
never have any physical meaning. These can therefore safely be ignored in the following, and
all bulk properties of the states are faithfully captured by the MPS tensors $A[m]$.

For infinite systems which are invariant under translations, it is natural to also impose
transation-invariance on the MPS {eq}`inf_mps_formula`. This leads to a *uniform* MPS which
has the same tensor $A[m] := A$ at every site. In the diagramatic notation introduced in the
[first section](tensor_networks), a uniform MPS can be represented as

```{image} /_static/figures/imps/umps.svg
:name: mps_state
:align: center
```

````{note}
In some cases, instead of assuming an MPS has the same tensor at each site it is more
natural to use a state with a non-trivial repeating unit cell. A uniform MPS with a unit
cell of size three would for example correspond to the state

```{image} /_static/figures/imps/umps3.svg
:name: mps_state3
:align: center
```

While we will restrict our discussion to MPS with a single-site unit cell, most concepts and
techniques apply just as well to the multi-site unit cell case.
````

One of the central objects in any unform MPS calculation is the transfer operator or
*transfer matrix*, defined in our case as

```{image} /_static/figures/imps/tm.svg
:name: transfer_matrix
:align: center
```

The transfer matrix corresponds to an operator acting on the space of $D\times D$ matrices,
and can be interpreted as a 4-leg tensor $\mathbb C^D \otimes \mathbb C^D \leftarrow \mathbb
C^D \otimes \mathbb C^D$. The transfer matrix can be shown to be a completely positive map,
such that its leading eigenvalue is a positive number. The eigenvalues of the transfer
matrix characterize the normalization and correlation length of a uniform MPS, while its
eigenvectors can be used to evaluate expectation values of local observables.


### Normalization

The norm of a uniform MPS corresponds to a contraction of the form

```{image} /_static/figures/imps/mpsNorm.svg
:name: mps_norm
:align: center
```

Clearly, this norm is nothing more than an infinite product of MPS transfer matrices defined
above. Consider the spectral decomposition of the $n$th power $\mathbb E^n$,

```{image} /_static/figures/imps/tmDecomp.svg
:name: tm_decomp
:align: center
```

where $l$ and $r$ are the left and right fixed points which correspond to the largest
magnitude eigenvalue $\lambda_0$ of $\mathbb E$,

```{image} /_static/figures/imps/fixedPoints.svg
:name: fixed_points
:align: center
```

and the $\lambda_i$ represent the remaining eigenvalues of smaller mangitude. Taking the
limit of this expression, it follows that the infinite product of transfer matrices reduces
to a projector onto the fixed points corresponding to the leading eigenvalue $\lambda_0$,

```{image} /_static/figures/imps/tmPower.svg
:name: tm_power
:align: center
```

To ensure a properly normalized state we should therefore rescale the leading eigenvalue
$\lambda_0$ to one by rescaling the MPS tensor as $A \leftarrow A / \sqrt{\lambda_0}$, as
well as normalize the fixed points $l$ and $r$ by requiring their trace to be equal to one:

```{image} /_static/figures/imps/traceNorm.svg
:name: trace_norm
:align: center
```

With these properties in place, the norm of an MPS reduces to the overlap between the
boundary vectors and the fixed points. Since there is no effect of the boundary vectors on
the bulk properties of the MPS, we can always choose these such that MPS is properly
normalized as $ \left \langle \Psi(\bar{A})\middle | \Psi(A) \right \rangle = 1$.


### Expectation Values

The fixed points of the transfer matrix can for example be to compute expectation values of
operators. Suppose we wish to evaluate expectation values of an extensive operator,

```{math}
O = \frac{1}{\mathbb{Z}} \sum_{n \in \mathbb{Z}} O_n.
```

If we assume that each $O_n$ acts on a single site and we are working with a properly
normalized MPS, translation invariance dictates that the expectation value of $O$ is given
by the contraction

```{image} /_static/figures/imps/expVal.svg
:name: exp_val
:align: center
```

In the uniform gauge, we can use the fixed points of the transfer matrix to contract
everything to the left and to the right of the operator, such that we are left with the
contraction

```{image} /_static/figures/imps/expVal2.svg
:name: exp_val2
:align: center
```

### Correlation Functions

Correlation functions are computed similarly. Let us look at

```{math}
c^{\alpha\beta}(m,n) = \bra{\Psi(\bar A)} (O^\beta_m)^\dagger O^\alpha_n \ket{\Psi(A)},
```

where $m$ and $n$ are abritrary locations in the chain, and, because of translation
invariance, the correlation function only depends on the difference $m-n$. Again, we
contract everything to the left and right of the operators by inserting the fixed points $l$
and $r$, so that

```{image} /_static/figures/imps/corrFunc.svg
:name: corr_func
:align: center
```

From this expression, we learn that it is the transfer matrix that determines the
correlations in the ground state. Indeed, if we again use the spectral decomposition of the
transfer matrix, recalling that now $\lambda_0 = 0$, we can see that the correlation
function reduces to

```{image} /_static/figures/imps/corrFunc2.svg
:name: corr_func2
:align: center
```

The first part is just the product of the expectation values of $O^\alpha$ and $O^\beta$,
called the disconnected part of the correlation function, and the rest is an exponentially
decaying part. This expression implies that connected correlation functions of an MPS
*always* decay exponentially, which is one of the reasons why MPS generally have a harder
time dealing with critical states. The correlation length $\xi$ is determined by the second
largest eigenvalue of the transfer matrix $\lambda_1$ as

```{math}
\xi = -\frac{1}{\log|\lambda_\mathrm{max}|}.
```

```{note}
The subleading eigenvalues of the transfer matrix typically also have a physical meaning,
because they correspond to subleading correlations in the system. For example, by focussing
on eigenvalues in a specific symmetry sector one can target the correlations associated to
exitations corresponding to that particular symmetry. The subleading eigenvalues also play a
crucial role in the powerful technique of *finite entanglement scaling* for infinite MPS
{cite}`rams2018precise`. Using this framework we can accurately capture critical phenomena
using MPS, despite the ansatz inherently having exponentially decaying correlations.
```


## Gauging Revisited

### Gauging in the Thermodynamic Limit

<!-- TODO: uniformize with finite MPS section -->

While a given MPS tensor $A$ corresponds to a unique state $\left | \Psi(A) \right \rangle$,
the converse is not true, as different tensors may give rise to the same state. This is
easily seen by noting that the gauge transform

```{image} /_static/figures/imps/gaugeTransform.svg
:name: gauge_transform
:align: center
```

leaves the physical state invariant. We may use this freedom in parametrization to impose
canonical forms on the MPS tensor $A$.

We start by considering the *left-orthonormal form* of an MPS, which is defined in terms of
a tensor $A_L$ that satisfies the condition

```{image} /_static/figures/imps/leftOrth.svg
:name: left_orthonormal
:align: center
```

We can find the gauge transform $L$ that brings $A$ into this form

```{image} /_static/figures/imps/leftGauge.svg
:name: left_gauge
:align: center
```

using an iterative procedure based on the QR docomposition, where starting from some initial
guess $L^0$ we repeatedly perform the QR-based update

```{image} /_static/figures/imps/qrStep.svg
:name: qr_step
:align: center
```

This iterative procedure is bound to converge to a fixed point for which
$L^{(i+1)}=L^{(i)}=L$ and $A_L$ is left orthonormal by construction:

```{image} /_static/figures/imps/qrConv.svg
:name: qr_convergence
:align: center
```

Note that this left gauge choice still leaves room for unitary gauge transformations

```{image} /_static/figures/imps/unitaryGauge.svg
:name: unitary_gauge
:align: center
```

which can be used to bring the right fixed point $r$ into diagonal form. Similarly, we can
find the gauge transform that brings $A$ into *right-orthonormal form*

```{image} /_static/figures/imps/rightGauge.svg
:name: right_gauge
:align: center
```

such that

```{image} /_static/figures/imps/rightOrth.svg
:name: right_orthonormal
:align: center
```

and the left fixed point $l$ is diagonal. A right-orthonormal tensor $A_R$ and a matrix $R$
such that $A R = R A_R$ can be found using a similar iterative procedure.

Finally, we can define a *mixed gauge* for the uniform MPS by choosing one site, the 'center
site', and bringing all tensors to the left of it in the left-orthonormal form and all the
tensors to the right of it in the right-orthonormal form. Defining a new tensor $A_C$ on the
center site, we obtain the form

```{image} /_static/figures/imps/mixedGauge.svg
:name: mixed_gauge
:align: center
```

By contrast, the original representation using the same tensor at every site is commonly
referred to as the *uniform gauge*. The mixed gauge has an intuitive interpretation.
Defining $C = LR$, this tensor then implements the gauge transform that maps the
left-orthonormal tensor to the right-orthonormal one, thereby defining the center-site
tensor $A_C$:

```{image} /_static/figures/imps/mixedGauge2.svg
:name: mixed_gauge2
:align: center
```

This relation is called the mixed gauge condition and allows us to freely move the center
tensor $A_C$ through the MPS, linking the left- and right orthonormal tensors.

Finally we may bring $C$ into diagonal form by performing a singular value decomposition $C
= USV^\dagger$ and absorbing $U$ and $V^\dagger$ into the definition of $A_L$ and $A_R$
using the residual unitary gauge freedom

```{image} /_static/figures/imps/diagC.svg
:name: mixed_gauge3
:align: center
```

````{note}
When working in the mixed gauge, the normalization of the MPS is entirely determined by that
of the center tensors $A_C$ and $C$. Indeed, it is easily seen that requiring that an MPS is
normalized now reduces to
```{image} /_static/figures/imps/normAC.svg
:name: norm_mixed
:align: center
```
or alternatively to ${\rm tr}(C^\dagger C) = 1$.
````

### Expectation Values Revisited

In the mixed gauge, we can locate the center site where the operator is acting, and then
contract everything to the left and right to the identity to arrive at the particularly
simple expression for the expectation value

```{image} /_static/figures/imps/expVal3.svg
:name: exp_val3
:align: center
```

### Entanglement Entropy

The mixed canonical form with a diagonal $C$ now allows to straightforwardly write down a
Schmidt decomposition of the state across an arbitrary bond in the chain

```{math}
\left | \Psi(A) \right \rangle = \sum_{i=1}^{D} C_i \left | \Psi^i_L(A_L) \right \rangle \otimes \left | \Psi^i_R(A_R) \right \rangle,
```

where the states $\left | \Psi^i_L(A_L) \right \rangle$ and $\left | \Psi^i_R(A_R) \right
\rangle$ are orthogonal states on half the lattice. The diagonal elements $C_i$ are exactly
the Schmidt numbers of any bipartition of the MPS, and as such determine its bipartite
entanglement entropy

```{math}
S = -\sum_i C_i^2 \log(C_i^2) .
```

### Truncation

The mixed canonical form also enables efficient truncatation of an MPS. The sum in the above
Schmidt decomposition can be truncated, giving rise to a new MPS that has a reduced bond
dimension for that bond. This truncation is optimal in the sense that the norm between the
original and the truncated MPS is maximized. To arrive at a translation invariant truncated
MPS, we can truncate the columns of the absorbed isometries $U$ and $V^\dagger$
correspondingly, thereby transforming *every* tensor $A_L$ or $A_R$. The truncated MPS in
the mixed gauge is then given by

```{image} /_static/figures/imps/truncMPS.svg
:name: truncate_mps
:align: center
```

We note that the resulting state based on this local truncation is not guaranteed to
correspond to the MPS with a lower bond dimension that is globally optimal. This would
require a variational optimization of the cost function.

```{math}
\left | \left | ~\left | \Psi(A) \right \rangle - \left | \Psi(\tilde{A}) \right \rangle ~\right | \right |^2.
```


### Code Example: `MPSKit.InfiniteMPS`

The Julia package [MPSKit.jl](https://github.com/maartenvd/MPSKit.jl) provides many tools
for working with infinite MPS. Without going into much detail, we can already check some
aspects of our discussion above with this numerical implementation.

We can construct an
[`MPSKit.InfiniteMPS`](https://maartenvd.github.io/MPSKit.jl/dev/lib/lib/#MPSKit.InfiniteMPS)
by specifying the physical and virtual vector spaces of the MPS. We will use standard
complex vector spaces as specified by a
[`TensorKit.ComplexSpace`](https://jutho.github.io/TensorKit.jl/latest/lib/spaces/#TensorKit.ComplexSpace),
and choose a physical dimension $d = 3$ and bond dimension $D = 5$.

```{code-cell} julia
using MPSKit, TensorKit

d = 3 # physical dimension
D = 5 # bond dimension
mps = InfiniteMPS(ℂ^d, ℂ^D)
```

The infinite MPS is automatically stored in the mixed canonical form introduced above. For
example, we can check that its normalization is indeed characterized by the center gauge
tensors $A_C$ and $C$.

```{code-cell} julia
using LinearAlgebra

@show norm(mps)
@show norm(mps.AC[1])
@show norm(mps.CR[1]);
```

We can also explicitly verify the mixed gauge conditions on $A_L$, $A_R$, $A_C$ and $C$ by
evaluating the corresponding tensor network diagrams using the
[`TensorOperations.@tensor` macro](https://jutho.github.io/TensorOperations.jl/stable/man/indexnotation/#The-@tensor-macro).

```{code-cell} julia
using TensorOperations

@tensor AL_id[-1; -2] := mps.AL[1][1 2; -2] * conj(mps.AL[1][1 2; -1])
@tensor AR_id[-1; -2] := mps.AR[1][-1 1; 2] * conj(mps.AR[1][-2 1; 2])

@assert AL_id ≈ id(space(mps.AL[1], 3)') "AL not in left-orthonormal form!"
@assert AR_id ≈ id(space(mps.AR[1], 1)) "Ar not in right-orthonormal form!"

@tensor LHS[-1 -2; -3] := mps.AL[1][-1 -2; 1] * mps.CR[1][1; -3]
@tensor RHS[-1 -2; -3] := mps.CR[1][-1; 1] * mps.AR[1][1 -2; -3]

@assert LHS ≈ RHS && RHS ≈ mps.AC[1] "Center gauge MPS tensor not consistent!"
```

We can also easily evaluate the expectation value of local operators

```{code-cell} julia
O = TensorMap(randn, ℂ^d ← ℂ^d)
expectation_value(mps, O)
```

as well as compute the correlation length encoded in the MPS.

```{code-cell} julia
correlation_length(mps)
```

MPSKit.jl exports a variety of infinite MPS algorithms, some of which will be discussed in
the next section.

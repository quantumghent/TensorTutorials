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

(symmetries)=
# Symmetries in Quantum Many-Body Physics

The goal of this section is to give a very gentle introduction to the concept of symmetries
in quantum many-body physics, and the notion of symmetric tensors. The general mathematical
framework of symmetries in physics (or at least the framework we will restrict to) is that
of group - and representation theory. Our goal is not to take this framework as a given and
illustrate it, but rather to first discuss a couple of important applications of symmetries
in the context of some concrete models and gradually build up to the more general framework.
We will finish our discussion with an outlook to generalizations of the framework presented
here. It goes without saying that we will only scratch the surface of this vast topic. The
interested reader is referred to the immense literature on this topic, or to a more
specialized course.

## Examples and Applications

### Symmetry Breaking, Order Parameters and Phases

Recall the one-dimensional transverse field Ising model defined above. Its degrees of
freedom are qubits ordered on a one-dimensional lattice, and its Hamiltonian reads

```{math}
H = -\sum_{i} \sigma^z_i\sigma^z_{i+1} -h_x\sum_i\sigma^x_i.
```

Let us simply consider periodic boundary conditions. Besides the obvious translation
symmetry, which we will discuss below, this model is also invariant under flipping all spins
simultaneously in the Z-direction, i.e. in the Pauli Z basis:
$\ket{\uparrow}\leftrightarrow\ket{\downarrow}$. That this operation constitutes a symmetry
is clear from the Hamiltonian as the energy of the first term only depends on the
neighbouring spins being (anti-)aligned, which is clearly spin flip-invariant. The second
spin is trivially invariant as this models an external magnetic field which is orthogonal to
the Z-direction.

This spin flip is "implemented", or more correctly "represented", by the unitary operator
$P=\bigotimes_i \sigma^x_i$. Notice that $P^2=1$ in accordance with our intuition that
flipping all the spins twice is equivalent with leaving all spins untouched. The fact that
this operator represents a symmetry of the model then translates to $[H,P]=0$, or
equivalently $P^\dagger HP=H$. Notice that the identity operator is also trivially a
symmetry (of every model) and thus the set $\{1,P\}$ is closed under taking the product.

Even though the Hamiltonian has the symmetry regardless of the value of the parameter $h_x$,
you might know from a previous course that the ground state or ground state subspace are not
necessarily invariant under the symmetry, a phenomenon known as spontaneous symmetry
breaking (SSB) or symmetry breaking for short. Let us investigate the ground state subspace
of the transverse field Ising model in the extremal case of vanishing and infinite magnetic
field.

- $h_x\rightarrow \infty$ In this case the model effectly reduces to a paramagnet. The
  unique ground state is the product state $\ket{\Psi_+}=\ket{+}^{\otimes N}$ where
  $\ket{+}=\frac{1}{\sqrt{2}}(\ket{\uparrow}+\ket{\downarrow})$ is the unique eigenvalue 1
  eigenvector of $\sigma^x$. Notice that this state is invariant under the symmetry operator
  $P$, $P\ket{\Psi_+}=\ket{\Psi_+}$. In other words, the ground state in this case is
  symmetric. For reasons mentioned below this state is also considered to be disordered.
- $h_x=0$ In this case the energy is minimized by aligning all the spins and the model
  behaves as a classical ferromagnet. Obviously, two distinct ground states are
  $\ket{\Psi_\uparrow}=\ket{\uparrow\uparrow...\uparrow}$ and
  $\ket{\Psi_\downarrow}=\ket{\downarrow\downarrow...\downarrow}$. Contrary to the previous
  case they span a two-dimensional ground state subspace, and these states are not
  symmetric. In fact, under the action of $P$ they get mapped onto the other:
  $P\ket{\Psi_\uparrow}=\ket{\Psi_\downarrow}$ and vice versa. The ground state in this case
  is thus symmetry broken, or ordered.

Since the ground state degeneracy is necessarily an integer, it is clear that it can not
change smoothly from two to one when the magnetic field is slowly turned on from
$h_x = 0 \rightarrow \infty$. Therefore the Ising model for small $h_x$ and large $h_x$ are
said to belong to different phases, and for some finite value of $h_x$ a phase transition
where the ground state degeneracy changes abruptly is expected to take place. As it turns
out, this change happens for $h_x = 1$, at which point the Ising model becomes critical.

Inspired by the credo of symmetry we can introduce a local operator which probes the phase
and can witness the phase transition. In the case of the Ising model this order parameter is
the local magnetisation on every site: $O=\sum_i\sigma^z_i$. It is clear that this order
parameter anticommutes with the symmetry, $P^\dagger OP=-O$, from which it follows that in
the symmetric phase the expectation value of the order parameter vanishes,
$\braket{\Psi_+|O|\Psi_+}=0$, while in the ferromagnetic phase
$\braket{\Psi_\uparrow|O|\Psi_\uparrow}>0$ and
$\braket{\Psi_\downarrow|O|\Psi_\downarrow}<0$. Notice however that for the latter we could
also have chosen the ground state $\ket{\Psi_\uparrow}+\ket{\Psi_\downarrow}$ in which case
the expectation value of $O$ becomes 0. So it seems that the expectation value of the order
parameter is ill-defined is this phase. This can be remedied by first adding a small
symmetry breaking term $\lambda\sum_i\sigma^z_i$ to the Hamiltonian which, depending on the
sign of $\lambda$, selects one of the ground states $\ket{\Psi_{\uparrow/\downarrow}}$ after
which the limit $\lambda\rightarrow 0$ is taken.

The synopsis of this example is thus the following. Symmetries in quantum many-body physics
(but also in single-particle quantum mechanics) are represented by unitary operators which
are closed under multiplication. Depending on the parameters in the Hamiltonian, part of
these symmetries can be broken by the ground state subspace, and this pattern of symmetry
breaking is a hallmark feature of different phases of the model. Different phases can be
probed by a local order parameter which does not commute with the symmetries. This paradigm
of classifying phases based on symmetry principles was first put forward by Landau
{cite}`Landau:1937obd`, and since then bears his name.

### Noether and Conserved Quantitites

You might remember Noether's theorem from a course on field theory. It states that every
continuous symmetry of a system (in field theory most often defined via its Lagrangian)
gives rise to a conserved current. In the context of quantum physics Noether's theorem
becomes almost trivial and states that the expectation value of every operator that commutes
with the Hamiltonian has a conserved expectation value:

```{math}
[H, O] = 0 \implies \frac{d}{dt}\braket{\Psi(t)|O|\Psi(t)} = 0.
```

The proof is almost trivial and is left as a simple exercise.

- The simplest example of this principle is obviously the Hamiltonain that trivially
  commutes with itself. The consequence is that the expectation value of the total energy is
  conserved.

- Another example is that of translation symmetry. Translation symmetry is implemented by
  the operator $T$ that acts on local operators $O_i$ via $T^\dagger O_iT=O_{i+1}$. Since
  for a system with $N$ sites we obviously have the identity $T^N=1$, and $T$ is unitary,
  the eigenvalues of $T$ are phases $\exp(2\pi ip/N)$ where the quantum number
  $p=0,1,...,N-1$ is the momentum. By virtue of Noether, translation invariance is
  understood to give rise to conservation of momentum, and thus momentum acts as a good
  quantum number for the eigenstates of translationally invariant models.

Let us consider another non-trivial example to illustrate the implications of this theorem.
Recall the spins $s$ XXZ Heisenberg model whose Hamiltonian reads

```{math}
H = -J\sum_i \left(S^x_iS^x_{i+1}+S^y_iS^y_{i+1}+\Delta S^z_iS^z_{i+1}\right).
```

The spin operators are $2s + 1$-dimensional and satisfy the $\mathfrak{su}(2)$ commutation
relations

```{math}
[\sigma^a_i,\sigma^b_j]=i\delta_{i,j}\sum_c \varepsilon_{abc}S^c_i 
```

Let us define the total spin

```{math}
S^a = \sum_i S^a_i.
```

From a direct computation it follows that in the case where $\Delta=1$, and the model thus
reduces to the Heisenberg XXX model, $H$ commutes with all $S^a$, $[H,S^a]=0$, $a=x,y,z$.
However, when $\Delta\neq 1$ only the Z component $S^z$ commutes with $H$, $[H, S^z]=0$.
Notice the difference with the Ising model where the same symmetry was present for all
values of $h_x$.

This means that in the $\Delta=1$ case the Hamiltonian is symmetric under the full $SU(2)$
(half integer s) or $SO(3)$ (integer s) symmetry (see below), whereas when $\Delta\neq 1$
only an $SO(2)\simeq U(1)$ symmetry generated by $S^z$ is retained. If $H$ commutes with
$S^z$ it follows that it automatically also commutes with $\exp(i\theta S^z)$,
$\theta\in[0,2\pi)$. This operator has an interpretation as a rotation around the Z-axis
with an angle $\theta$.

According to Noether the Heisenberg model thus has conserved quantities associated with
these operators. Regardless of $\Delta$ the Z component of the total spin is conserved, and
for $\Delta=1$ all components of the total spin are conserved. In particular this means that
the eigenvalue $M_z$ of $S^z$ and $S(S+1)$ of $\vec{S}\cdot\vec{S}$ are good quantum numbers
to label the eigenstates of the Heisenberg Hamiltonian.

## Group and Representation Theory

Motivated by the examples from above, we will gently introduce some notions of group - and
representation theory that form the backbone of a general theory of symmetries.

### Group Theory

Roughly speaking a group $G$ is a set of symmetry operators and a multiplication rule on how
to compose them. Let us motivate the definition one step a time.

First of all notice that a model can have a finite or infinite (discrete or continuous)
number of symmetries. Clearly, the spin flip symmetry of the Ising model consists of only
one non-trivial symmetry operation, namely flipping all spins. The operator carrying out
this transformation is $P=\bigotimes_i\sigma^x_i$. The XXZ model however has a continuous
symmetry, namely rotations around the Z-axis, that is implemented via $\exp(i\theta S^z)$,
$\theta\in[0,2\pi)$, where we should really think about every value of $\theta$ as labeling
a different symmetry operation.

These symmetries can be composed or multiplied to form a new symmetry operation. Take for
example flipping all the spins. Flipping all spins twice results in not flipping any spins
at all, which is trivially also a symmetry of the Hamiltonian. Next, consider also the
$U(1)$ symmetry of the XXZ model. First rotating over $\theta_2$ and then over $\theta_1$
gives a new rotation over $\theta_1+\theta_2$: $\exp(i\theta_1 S^z)\exp(i\theta_2
S^z)=\exp(i(\theta_1+\theta_2) S^z)$. This leads to the first part of the definition of what
a group is.

1. A group $G$ is a set $G=\{g_1,g_2,...\}$ endowed with a multiplication $G\times
   G\rightarrow G$. There exists an identity $1\in G$ for the multiplication such that
   $1g=g1=g, \forall g\in G$.

Note that this multiplication is not necessarily abelian. A simple example is the full
$SU(2)$ symmetry of the XXX model defined above. However, the composition of symmetries is
still associative:

2. For all group elements $g,h,k$ we have that $g(hk)=(gh)k$.

A property we also would like to formalize is the fact that every symmetry transformation
can be undone. Take for example a $U(1)$ rotation $\exp(i\theta S^z)$, if we compose it with
the opposite rotation $\exp(i(2\pi-\theta) S^z)$ we get the identity. Hence:

3. Every group element $g$ has a unique inverse $g^{-1}$: $gg^{-1}=g^{-1}g=1$.

Together 1. 2. and 3. constitute the definition of a group. Before mentioning some examples
let us also introduce the concept of a subgroup. As the name suggests, a subgroup is a
subset of a group which itself constitutes a group. Note for example that a rotation over
$\pi$, $\exp(i\pi S^z)$, together with the identity, generates a subgroup of $\{\exp(2\pi
i\theta S^z|\theta\in[0,2\pi)\}$ with two elements.

The concept of subgroups lies at the heart of symmetry breaking. Recall that in the
ferromagnetic phase, the Ising model breaks the spin flip symmetry. In Landau's paradigm we
say that the pattern of symmetry breaking is $\mathbb{Z}_2\rightarrow \{1\}$ (see below for
an explanation of the notation). In other words, the full symmetry group ($\mathbb{Z}_2$) is
broken in the ferromagnetic phase to a subgroup (the trivial group). More generally, a
theory with a $G$ symmetry can undergo a pattern of symmetry breaking $G\rightarrow H$ where
$H$ is a subgroup of $G$. The meaning of this symbolic expression is that the ground states
keep an H symmetry, and the ground state degeneracy is $|G|/|H|$. 

#### Examples

- The trivial group is a group with only one element that is than automatically also the
  identity, and a trivial multiplication law. Above, it was denoted by $\{1\}$.
- $\mathbb{Z}_N$ is the additive group of integers modulo $N$. The group elements are the
  integers $\{0,1,...,N-1\}$ and the group multiplication is addition modulo $N$. Hence it
  is clearly a finite group. In particular, the spin flip symmetry from above corresponds to
  the group $\mathbb{Z}_2$. Notice that for all $N$ $\mathbb{Z}_N$ is abelian.
- Another abelian group is $U(1)$. This group is defined as
  $U(1)=\left\{z\in\mathbb{C}:|z|^2 = 1\right\}$, with group multiplication the
  multiplication of complex numbers. Note we encountered this group in the XXZ model as
  being the rotations around the Z axis: $\{\exp(2\pi i\theta S^z|\theta\in[0,2\pi)\}$.
- $SU(2)$ is the group of unimodular unitary $2\times 2$ matrices:
  ```{math}
  SU(2) := \left\{U \in \mathbb{C}^{2\times 2} | 
    \det U = 1, UU^\dagger = U^\dagger U = \mathbb{I}\right\}.
  ```
  The group multiplication is given by group multiplication. Similarly, one defines
  $SU(N),N\geq 2$. Note that none of these groups are abelian.
  
- The 3D rotation group or special orthogonal group $SO(3)$ is the group of real $3\times 3$
  orthogonal matrices with unit determinant:
  ```{math}
  SO(3) := \left\{M\in\mathbb{R}^{3\times 3}|MM^T=M^TM=\mathbb{I},\det M=1\right\}.
  ```
  Similarly, one defines $SO(N),N\geq 2$. Note that only $SO(2)$ is abelian.

(representation_theory)=
### Representation Theory

In the above examples, we were dealing with the question which symmetry transformations
leave the Hamiltonian (and in the absence of symmetry breaking also the ground states)
invariant. These symmetry representations were implemented (represented) by invertible
linear operators, non-singular matrices, that form a closed set under multiplication. This
multiplication structure is what we identified as a group. What we could now do, is to take
a group as given, and wonder which linear transformations we can come up with that multiply
according to these multiplication rules. This is exactly the underlying idea of
representation theory. Representation theory deals with the question how groups can linearly
act on vector spaces.

This immediately raises a plethora of questions such as if we can classify all
representations (up to some kind of equivalence), if there exists 'minimal' representations
and how we can construct new representations of known ones. A minimal answer to these
questions is the goal of this section.

#### Definition

For the sake of these notes, a representation of a group $G$ is thus a set of matrices
indexed by the group elements, $\{X_g|g\in G\}$ that multiply according to the
multiplication rule of $G$:

```{math}
X_gX_h = X_{gh}
```

Note that the identity is always mapped to the identity matrix!

We call the dimension of the matrices $X_g$ the dimension of the representation.

##### Examples

- Every group can be trivially represented by mapping every group element to the 'matrix'
  (1). Obviously, this representation is one-dimensional and is called the trivial
  representation.
- Probably the simplest non-trivial representation, is the representation of $\mathbb{Z}_2$
  that maps the non-trivial element to -1. Concretely, $X_0=1, X_1=-1$, and indeed
  $X_1X_1=(-1)^2=X_0$. This representation is called the sign representation.
- Let us construct a two-dimensional representation of $\mathbb{Z}_2$. Since the Pauli
  matrix $\sigma^x$ (as any other Pauli matrix) squares to the identity, $\sigma^x$ together
  with the two-dimensional identity matrix constitutes a two-dimensional representation of
  $\mathbb{Z}_2$. In the notation from above, $X_0=\mathbb{I}_2$, $X_1=\sigma^x$. This
  representation is called the regular representation of $\mathbb{Z}_2$.

#### Complex Conjugate Representation, Tensor Product and Direct Sum Representation

Given a representation $\{X_g|g\in G\}$, the complex conjugate representation $\bar X$ is
defined as $\bar X=\{\bar X_g|g\in G\}$ which satisfies the defining property of
representations via $\bar X_g\bar X_h= \overline{X_gX_h}=\bar X_{gh}$.

Given two representations of $G$, $X\equiv\{X_g|g\in G\}$ and $Y\equiv\{Y_g|g\in G\}$, there
are two obvious ways to construct a new representation.

The first one is the tensor product representation defined via the Kronecker product of
matrices:

```{math}
\{X_g\otimes Y_g|g\in G\}.
```

You should check that these still satisfy the defining property of a representation. The
dimension of the tensor product is the product of the dimensions of the two representations
$X$ and $Y$.

The other one is the direct sum:

```{math}
\{X_g\oplus Y_g|g\in G\}.
```

Its dimension is that the sum of the dimensions of $X$ and $Y$.

#### Irreducible Representations

It is clear that physical observables should not depend on any choice of basis. Therefore
two representations are (unitarily) equivalent when there is a unitary basis transformation
$U$ such that $X_g' =UX_gU^\dagger$. Note that $U$ is independent of $g$.

##### Example

Consider again the two-dimensional regular representation of $\mathbb{Z}_2$ from above. The
basis transformation

```{math}
H=\frac{1}{\sqrt 2}
\begin{pmatrix}
    1 & 1\\
    1 & -1
\end{pmatrix}
```

shows that this representation is equivalent to one where the non-trivial element of
$\mathbb{Z}_2$ is represented by $H\sigma^x H^\dagger=\sigma^z$. This illustrates that the
regular representation is equivalent to the direct sum of the trivial representation and the
sign representation!

The crux of this example is the following. Some representations can, by an appropriate
choice of basis, be brought in a form where all $X_g$ are simultaneously block-diagonal:

```{math}
X_g'=UX_gU^\dagger=
\begin{pmatrix}
    \fbox{$X^1_g$} & 0 &\cdots\\
    0& \fbox{$X^2_g$} & \cdots\\
    \vdots & \vdots & \ddots
\end{pmatrix}.
```

These blocks correspond to invariant subspaces of the representation, i.e. subspaces that
transform amongst themselves under the action of the group.

An irreducible representation, irrep for short, can then be defined as a representation that
can not be brought in a (non-trivial) block-diagonal form by any change of basis.

It can be shown that every finite group has a finite number of irreps. The sum of the
dimensions squared is equal to the number of elements in the group: $\sum_\alpha
d_\alpha^2=|G|$ where the sum is over all irreps labeled by $\alpha$ and $d_\alpha$ denote
their respective dimensions.

One of the key questions of representation theory is what the irreps of a given group are
and how the tensor product of irreps (which is in general not an irrep!) decomposes in a
direct sum of irreps. The latter are sometimes known as the fusion rules. The basis
transformation that reduce a given representation in a direct sum of irreps is sometimes
called the Clebsch-Gordan coefficients, and are for some groups known explicitly. Before
discussing the example of $SU(2)$, let us first state the most important result in
representation theory which is due to Schur.

[Schur's lemma] If a matrix $Y$ commutes with all representation matrices of an irreducible
representation of a group G, $X_gY=YX_g$ $\forall g\in G$, then $Y$ is proportional to the
identity matrix.

(su2_irreps)=
#### Example

The answer to the questions posed above is very well understood for the case of $SU(2)$. You
probably know the answer from a previous course on quantum mechanics.

The irreps of $SU(2)$ can be labeled by its spin, let us call it $s$, that takes values
$s=0,1/2,1,3/2,...$. The dimension of the spin $s$ representation is equal to $2s+1$, so
there is exactly one irrep of every dimension. The spin $s=0$ irrep corresponds to the
trivial representation.

The fusion rules can be summarized as

```{math}
s_1\otimes s_2 \simeq \bigoplus_{s=|s_1-s_2|}^{s_1+s_2}s.
```

For example: $\frac{1}{2}\otimes\frac{1}{2}\simeq 0\oplus 1$. The Clebsch-Gordan
coefficients for $SU(2)$ have been computed analytically, and for low-dimensional irreps
have been tabulated for example
[here](https://pdg.lbl.gov/2018/reviews/rpp2018-rev-clebsch-gordan-coefs.pdf).

(symmetric_tensors)=
## Symmetric Tensors

In physics we are often dealing with tensors that transform according to the tensor product
representation of a given group $G$. A symmetric tensor can then be understood as a tensor
that transforms trivially under the action of $G$, or more concretely under the tensor
product representation $X\otimes\bar Y\otimes\bar Z$:

```{figure} ../_static/SymmetricTensors/symmtens.svg
:scale: 12%
:name: symmtens
:align: center
```

This has strong implications for the structure of the tensor $T$. Notice that we didn't
assume the representations $X,Y$ and $Z$ to be irreducible. As we argued above, an
appropriate change of basis can bring the representations $X,Y$ and $Z$ in block-diagonal
form where every block corresponds to an irrep of the group and every block can appear
multiple times, which we call the multiplicity of an irrep in the representation. Schur's
lemma then implies that in this basis, the tensor becomes block-diagonal. In an appropriate
matricization of $T$ we can thus write $T=\bigoplus_c B_c\otimes\mathbb{I}_c$ where the
direct sum over $c$ represents the decomposition of $X\otimes\bar Y\otimes\bar Z$ in irreps
$c$ that can appear multiple times. In other words, the generic symmetric tensor $T$ can be
stored much more efficiently by only keeping track of the different blocks $B_c$.

TensorKit is particularly well suited for dealing with symmetric tensors. What TensorKit
does is exactly what was described in the previous paragraph, it keeps track of the block
structure of the symmetric tensor, hereby drastically reducing the amount of memory it takes
to store these objects, and is able to efficiently manipulate them by exploiting its
structure to the maximum.

As a simple exercise, let us construct a rank 3 $SU(2)$ symmetric tensor as above. For
example the spin $1/2$ and spin $1$ representation can be called via respectively

```{code-cell} julia
:tags: ["hide-output"]
using TensorKit

s = SU₂Space(1/2 => 1)
l = SU₂Space(1 => 1)
```

Here, `` => 1`` essentially means that we consider only one copy (direct summand) of these
representations. If we would want to consider the direct sum $\frac{1}{2}\oplus\frac{1}{2}$
we would write

```{code-cell} julia
:tags: ["hide-output"]
ss = SU₂Space(1/2 => 2)
```

A symmetric tensor can now be constructed as

```{code-cell} julia
A = TensorMap(l ← s ⊗ s)
```
This tensor then has, by construction, the symmetry property that it transforms trivially
under $1\otimes\bar{\frac{1}{2}}\otimes\bar{\frac{1}{2}}$. The blocks can then be inspected
by calling ``blocks`` on the tensor, and we can also check that the dimensions of the domain
and codomain are as expected:

```{code-cell} julia
@assert dim(domain(A)) == 4
@assert dim(codomain(A)) == 3
blocks(A)
```

We see that this tensor has one block that we can fill up with some data of our liking. Let
us consider another example

```{code-cell} julia
B = TensorMap(s ← s ⊗ s)
blocks(B)
```

This tensor does not have any blocks! This is compatible with the fact that two spin 1/2's
cannot fuse to a third spin 1/2. Finally let us consider a tensor with with more blocks:

```{code-cell} julia
C = TensorMap(ss ← ss)
blocks(C)
```

This tensor has four non-trivial entries.

## Outlook and generalizations

Let us conclude with an outlook and some generalizations.

- Besides the "global" symmetries we considered here, you might also be familiar with gauge
  symmetries from another course. Gauge theories are ubiquitous in physics and describe a
  plethora of interesting physical phenomena. Gauge symmetries should however not be thought
  of as actual symmetries transforming physically different states into each other, but
  rather describe a redundancy in the description of the system. Nevertheless, group theory
  also lies at the heart of these theories.
- In this brief overview we mostly neglected spatial symmetries. Spatial symmetries can be
  understood as transformations that translate, rotate or reflect the lattice. These kind of
  symmetries thus don't act "on site" anymore. The full classification of spatial symmetry
  groups is notoriously rich and beautiful, especially in higher dimensions, and exploiting
  them in algorithms can result in tremendous speedup and stability. We already encountered
  the example of translation symmetry. One of the benefits of exploiting this symmetry in
  tensor networks is e.g. that if the ground state of an infinite one-dimensional model does
  not break translation invariance, this ground state can be well modelled by a uniform
  matrix product state, a matrix product state consisting of one tensor repeated
  indefinitely.
- Inspired by the discovery of topological phases of matter and their anyonic excitations,
  there has been a growing fascination with the exploration of non-invertible, or
  categorical symmetries. These symmetries are beyond the scope of these notes. These
  categorical symmetries are not described by groups but by more general and intricate
  algebraic structures called fusion categories, of which (finite) groups and their
  representations are specific examples. For an example of how spin chains with categorical
  symmetries can be constructed, see for example {cite}`feiguin2007interacting`. TensorKit
  allows for an efficient construction and storage of tensors which are symmetric with
  respect to these more general kind of symmetries.
  

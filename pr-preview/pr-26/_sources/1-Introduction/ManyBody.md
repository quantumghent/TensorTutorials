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

(many_body)=
# The Hilbert Space of Many-Body Physics

All of the previous axioms remain valid for a composite system consisting of several quantum
degrees of freedom. However, we need to know how to describe the state of the system, and 
thus more specifically, how to define the Hilbert space associated to such a system. It
turns out that quantum mechanics forces us to distinguish two cases.

## Distinguisable Particles and Tensor Products

Consider a quantum system composed out of two subsystems, which we call $A$ and $B$,
sometimes referred to as Alice and Bob in quantum information context. These can themselves
already be many-body systems. Suppose we know the Hilbert space $\mathbb{H}^A$ in which to
describe states of subsystem $A$ when considered as an isolated system on itself, and
analoguously for $\mathbb{H}^B$. Now consider both systems together, but where they do not
interact, so that we can still treat them independently. In particular, we can prepare
subsystem $A$ in a state $\ket{\psi^A}$ and subsystem $B$ in a state $\ket{\psi^B}$. We
should also be able to describe these two independent subsystems jointly, so that there must
exist a map from the two arguments $(\ket{\psi^A}, \ket{\varphi^B}) \in \mathbb{H}^A \times
\mathbb{H}^B$ to a single state which we denote as $\ket{ \psi^A} \otimes \ket{\varphi^B}$
and that lives in a joint Hilbert space $\mathbb{H}^{AB}$ that we have yet to determine.

Now, it makes sense that, if we build superpositions in one of the two subsystems, while
keeping the other fixed, this also correspond to a superposition in the joint description of
both systems together. This leads to

```{math}
\left (a_1 \ket{\psi^A_1} + a_2 \ket{\psi^A_2}\right) \otimes \ket{\varphi^B} = a_1
\ket{\psi^A_1 }\otimes \ket{\varphi^B} + a_2 \ket{\psi^A_1 }\otimes \ket{\varphi^B}
```

and similarly

```{math}
\ket{\psi^A} \otimes \left(b_1 \ket{\varphi^B_1} + b_2 \ket{\varphi^A_2}\right) = b_1
\ket{\psi^A }\otimes \ket{\varphi^B_1} + b_2 \ket{\psi^A }\otimes \ket{\varphi^B_2}.
```

Hence, the Hilbert space $\mathbb{H}^{AB}$ that we are trying to construct must contain all
states $\ket{\psi^A} \otimes \ket{\varphi^B}$ for all $\ket{\psi^A} \in \mathbb{H}^A$ and
all $\ket{\varphi^B} \in \mathbb{H}^B$, all possible linear combinations thereof (in order
to be a vector space), but in such a way that the above equalities hold. This construction,
which can be made mathematically precise, is known as the tensor product of vector spaces
$\mathbb{H}^{AB} = \mathbb{H}^A \otimes \mathbb{H}^B$.

We have also denoted the output of the map from two states $(\ket{\psi^A}, \ket{\varphi^B})
\in \mathbb{H}^A \times \mathbb{H}^B$ to $\mathbb{H}^A \otimes \mathbb{H}^B$ using the same
tensor product symbol, and refer to such a state as a (tensor) product state $\ket{ \psi^A}
\otimes \ket{\varphi^B}$. Importantly, however, the tensor product space $\mathbb{H}^A
\otimes \mathbb{H}^B$ certainly contains vectors which are not product states, such as

```{math}
a_1 \ket{ \psi_1^A} \otimes \ket{\varphi_1^B} + a_2 \ket{ \psi_2^A} \otimes
\ket{\varphi_2^B}.
```

This forms the basis for quantum correlations and the concept of (quantum) entanglement,
which will be a fundamental property of quantum many-body systems. That the Hilbert space of
a composite system is given by the tensor product of the individual Hilbert spaces is often
introduced as a separate axiom. The deductive (but informal) argument just given can however
be turned into a proof that depends only on the axioms given above (in fact only on the
first two).

As expected (and required), it can be shown that the tensor product of two Hilbert spaces is
again a Hilbert space, if we define its inner product in the following way. We first define
the inner product for product states as

```{math}
\braket{\psi_1^A \otimes \varphi_1^B | \psi_2^A \otimes \varphi_2^B } = \braket{\psi_1^A |
\psi_2^A} \braket{\varphi_1^B | \varphi_2^B}
```

and then extend this definition by linearity (in the second argument and antilinearity in
the first argument).

In practice, given two finite-dimensional Hilbert spaces $\mathbb{H}^A \cong
\mathbb{C}^{d^A}$ and $\mathbb{H}^B \cong \mathbb{C}^{d^B}$ with a basis $\{ \ket{j},
j=1,\dots, d^A\}$ and $\{\ket{k}, k=1,\dots, d^B\}$, the tensor product space is spanned by
a basis composed of all products

```{math}
\{ \ket{j,k} = \ket{j} \otimes \ket{k}, j=1,\dots, d^A, k=1,\dots, d^B\}
```

and thus has dimension $d^A \cdot d^B$. A general state $\ket{\Psi} \in
\mathbb{H}^{A}\otimes \mathbb{H}^B$ can then be expanded as

```{math}
\ket{\Psi} = \sum_{j=1}^{d^A }\sum_{k=1}^{d^B} \Psi_{jk} \ket{j,k}
```

The expansion coefficients $\Psi_{jk}$ thus have two indices, and it is often useful to
think of them as a matrix. Note that we will almost always use this product basis, also
referred to as the *computational basis*, for working with tensor product spaces. However,
one can certainly also use more complicated basis choices, where the basis vectors are not
simple product states. One well known choice that you might remember from your quantum
mechanics course is in the case of two spin-1/2 systems. If we denote the basis for a single
spin-1/2 system as $\{\ket{\uparrow},\ket{\downarrow}\}$, then the product basis for a
system consisting of two spin-1/2 systems is given by $\{\ket{\uparrow,\uparrow},
\ket{\downarrow,\uparrow}, \ket{\uparrow,\downarrow}, \ket{\downarrow,\downarrow}, \}$.
However, in the context of spin coupling (see Section on Symmetries), one also uses the
coupled basis

```{math}
\ket{0,0} &= \frac{1}{\sqrt{2}} \left(\ket{\uparrow,\downarrow} - \ket{\downarrow,\uparrow}\right)\\
\ket{1,+1} &= \ket{\uparrow,\uparrow}\\
\ket{1,0} &= \frac{1}{\sqrt{2}} \left(\ket{\uparrow,\downarrow} + \ket{\downarrow,\uparrow}\right)\\
\ket{1,-1} &= \ket{\downarrow,\downarrow}
```

Note that we also use the same tensor product notation as an operation to map operators from
the subsystems into operators acting on the full tensor product Hilbert space. In
particular, the process of measuring operator $\hat{A}$ in subsystem $A$ and simultaneously
operator $\hat{B}$ in subsystem $B$ is associated with an operator $\hat{A}\otimes \hat{B}$
acting on $\mathbb{H}^A \otimes \mathbb{H}^B$, the action of which is first defined on the
product states as

```{math}
\left(\hat{A} \otimes \hat{B}\right) \left(\ket{\psi^A}\otimes \ket{\varphi^B}\right) = \left(\hat{A}\ket{\psi^A}\right) \otimes \left(\hat{B}\ket{\varphi^B}\right)
```

and then extended by linearity. It furthermore holds that

```{math}
(\hat{A}_1 \otimes \hat{B}_1) (\hat{A}_2 \otimes \hat{B}_2) = (\hat{A}_1 \hat{A}_2) \otimes (\hat{B}_1 \hat{B}_2).
```

With respect to a product basis, the matrix representation of $\left(\hat{A} \otimes
\hat{B}\right)$ is given by the
[Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product).

When we are only interested in an operator $\hat{O}$ acting on subsystem $A$ without doing
anything on subsystem $B$, we should create the operator $\hat{O} \otimes \hat{1}_B$, with
$\hat{1}_B$ the identity operator of the Hilbert space $\mathbb{H}^B$. Often, we will omit
this explicit tensor product with the identity operator, and simply use some notation which
indicates that an operator acts on a certain subsystem, such as $\hat{O}^{(A)} = \hat{O}
\otimes \hat{1}_B$. This also makes it explicit that operators defined on different
subsystems, when lifted to act on the full Hilbert space, commute, i.e.

```{math}
\left[\hat{O}_1^{(A)} , \hat{O}_2^{(B)}\right] =
\left[ \hat{O}_1 \otimes \hat{1}_B, \hat{1}_A \otimes \hat{O}_2\right] = 0.
```

The tensor product construction extends readily to systems with multiple subsystems.
Consider for example a system consisting of qubits, where every individual qubit has an
associated Hilbert space $\mathbb{C}^2$ with basis denoted as $\{\ket{0},\ket{1}\}$. The
Hilbert space $\mathbb{H}^N$ of $N$ qubits is then spanned by a computational basis which we
can denote as

```{math}
\{\ket{s_1, s_2, \ldots, s_N} = \ket{s_1} \otimes \ket{s_2} \otimes \cdots \otimes
\ket{s_N}; s_1 =0,1; s_2 =0,1; \ldots; s_n =0,1\}.
```

Hence, the Hilbert space thus has dimension $2^N$, and a general state $\ket{\Psi}$ has
expansion coefficients

```{math}
\Psi_{s_1,s_2, \ldots, s_N}
```

which can be interpreted as a single vector of length $2^N$, or as a $N$-dimensional tensor,
where every tensor index ranges over the two values 0 and 1. This exponential increase of
the Hilbert space dimension with the number of particles is exactly why the quantum
many-body problem is so difficult, but also essential for providing a quantum computer with
its speed-up. It is exactly these type of quantum states living in a many-body Hilbert
space, which is thus composed of many tensor product factors, that we will represent as a
tensor network.

Finally, we also have to specify the Hamiltonian of a many-body system. It typically takes
the form of a sum of terms, where every individual term acts nontrivially on only a few
subsystems. One important example that will reappear throughout these tutorials is the
"Quantum Ising Model with transverse magnetic field", which acts on a system composed of
qubits or spin-1/2 particles, and is defined as

```{math}
\hat{H} = - J \sum_{\langle i, j \rangle} \sigma^z_i \otimes \sigma^z_j - h \sum_i
\sigma^x_i
```

Here, the summation variables $i$ and $j$ correspond to the sites of a lattice. The notation
$\sum_{\langle i,j \rangle}$ denotes a sum over pairs of neighbouring lattice sites $i$ and
$j$. The second sum contains terms $\sigma^x_i$ which act nontrivially only on the site $i$,
and as the identity operator elsewhere. If, for example, we enumerate the sites from $1$ to
$N$, it would act as

```{math}
\sigma^x_i = \underbrace{1 \otimes 1 \otimes \ldots \otimes 1}_{\text{$i-1$ factors}} \otimes \sigma^x \otimes \underbrace{1 \otimes \ldots \otimes 1}_{\text{$N-i-1$ factors}}
```

with $\sigma^x = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}$ the Pauli x matrix, and $1$
the $2 x 2$ unit matrix. The first set of terms in $\hat{H}$ acts nontrivially on two sites,
and is defined analoguously, using the Pauli z matrices $\sigma^z = \begin{bmatrix} 1 & 0 \\
0 & -1 \end{bmatrix}$.

## Identical Particles and Pauli's Exclusion Principle

The tensor product construction needs to be revised when discussing the Hilbert space of a
system composed of identical particles. Consider for example a system made out of $N$
identical particles. To every individual particle we can associate a particular Hilbert
space, which we denote as $\mathbb{H}^{(1)}$, for example $\mathbb{H}^{(1)} =
L^2(\mathbb{R})$ for a particle moving on the real line, or $\mathbb{H}^{(1)} =
\mathbb{C}^L$ for a particle living on the sites of a chain of length $L$.

If we temporarily assign each of the $N$ particles a label $n=1, \dots, N$, then the Hilbert
space of the composite system would be given by the $N$-fold tensor product
$\widetilde{\mathbb{H}}^{(N)} = \left(\mathbb{H}^{(1)}\right)^{\otimes N}$. However, for
identical particles, our labeling is completely arbitrary. For the case of $N=2$ particles
on a chain of $L$ sites, we cannot distinguish between the state $\ket{j_1, j_2}$ where
particle $1$ is on site $j_1$ and particle $2$ is on site $j_2$ versus the state $\ket{j_2,
j_1}$ where site $j_1$ is occupied by the particle that we gave label $2$ and site $j_2$ is
occupied by the particle with label $1$. A general redefinition of the particle labels
amounts to a permutation, and we have to require that no physical measurement can
distinguish between such permutations. Hence, this permutation invariance does not behave
like a regular symmetry (like e.g. rotation symmetry, one can still construct observables
along preferred directions such that they can detect rotations).

We are forced to restrict our tensor product Hilbert space
$\left(\mathbb{H}^{(1)}\right)^{\otimes N}$ to the subspace $\mathbb{H}^{(N)}$ of physical
states which are not affected by acting with such permutations. Note that, due to the fact
that quantum states actually correspond to rays of vectors, it is still allowed that the
vectors in $\mathbb{H}^{(N)}$ pick up a phase factor when applying certain permutations. It
is a result in the representation theory of the permutation group that there are only two
possibilities. Either the phase factor is always absent (or thus 1), or the phase factor is
(-1) for odd permutations and (+1) for even permutations, i.e. the phase factor equals the
sign(ature) of the permutation. Identical particles for which the phase factor is always one
are known as *bosons*, whereas those with the nontrival phase factor choice correspond to
*fermions*. Indeed, the nontrivial phase factor automatically gives rise to *Pauli's
exclusion principle*: two fermions cannot be in the same quantum state, since $P_{12}
\ket{j_1,j_2} = \ket{j_2,j_1} = -\ket{j_1,j_2}$ and for $j_1=j_2$ we would thus find
$\ket{j,j} = -\ket{j,j}$.

Bosons are thus described by states which are symmetric under permutations, whereas fermions
are described by states which are called antisymmetric. We can define an operator on
$\tilde{\mathbb{H}}^{(N)} = \left(\mathbb{H}^{(1)}\right)^{\otimes N}$ that maps any given
state onto such a (anti)symmeric state, namely by first defining its action on product
states as

```{math}
\hat{S}^{\pm} \ket{\psi_1} \otimes \ket{\psi_2} \otimes \cdots \otimes \ket{\psi_N} =
\frac{1}{\sqrt{N!}} \sum_{\sigma \in S_N} \epsilon_\sigma \ket{\psi_{\sigma(1)}} \otimes
\ket{\psi_{\sigma(2)}} \otimes \cdots \otimes \ket{\psi_{\sigma(N)}}
```

and then extending it by linearity. Here, $S_N$ is the symmetric group containing all
permutations $\sigma$ of $N$ elements, where the permutation $\sigma$ is a bijective map
from integers $j \in \{1,\dots,N\}$ to a new number $\sigma(j) \in \{1,\dots,N\}$. The
sign(ature) $\epsilon_\sigma$ of the permutation takes the value $+1$ or $-1$, depending on
whether the permutation $\sigma$ can be obtained by composing an even or odd number of
elementary transpositions. An elementary transposition $\tau_{i,j}$ is a permutation which
only interchanges the two numbers $i$ and $j \neq i$:

```{math}
\tau_{i,j}(i) =j, \tau_{i,j}(j) = i, \tau_{i,j}(k) =k, \forall k\neq i \land k \neq j
```

Note that $\hat{S}^{\pm}$ does not necessarily yield a normalised state, and can indeed even
map a state to zero, in order to give rise to Pauli's exclusion principle:
$\hat{S}^-\ket{j,j} = 0$. The image of $\hat{S}^{\pm}$ contains all states with the proper
behaviour under relabeling permutations, and thus correspond to the physical Hilbert space
for bosons or fermions:

```{math}
\mathbb{H}^{(N)} = \hat{S}^{\pm} \widetilde{\mathbb{H}}^{(N)} = \hat{S}^{\pm} \left(\mathbb{H}^{(1)}\right)^{\otimes N}
```

Note that in this case, the physical Hilbert space is not a tensor product. However, we can
think of it as a subspace of an auxiliary Hilbert space, $ \widetilde{\mathbb{H}}^{(N)}$,
which is a tensor product. The restriction to this subspace can thus be thought of as a
constraint, and the same scenario happens in other constrained quantum systems. The most
notable example is that of quantum gauge theories, where there is an extensive set of
constraints, namely that physical quantum states need to be gauge invariant.

Now consider a single particle Hilbert space $\mathbb{H}^{(1)}$ with an orthonormal basis
$\{\ket{j}, j=1,\ldots,L\}$, for example where $\ket{j}$ corresponds to the particle being
positioned on site $j$ of a lattice with $L$ sites. We also refer to these single particle
states as modes. To construct a basis for $\mathbb{H}^{(N)}$, we can start from the tensor
product basis of $\widetilde{\mathbb{H}}^{(N)}$ and apply $\hat{S}^{\pm}$ to each of its
$L^N$ elements. Let us henceforth denote these states as

```{math}
\ket{j_1,j_2,\ldots ,j_N} = \hat{S}^{\pm} \left(\ket{j_1} \otimes \ket{j_2} \otimes \cdots
\otimes \ket{j_N}\right)
```

The application of $\hat{S}^{\pm}$ will create certain linear dependences. In particular,
states $ \ket{j_1,j_2, \ldots, j_N}$ that contain the same set of modes $j_k$, i.e. for
which the $j_k$'s are related by a permutation, are equal (up to a sign in the case of
$\hat{S}^-$). We can thus select a single state by ordering the $j_k$ arguments.
Furthermore, in the case of $\hat{S}^{-}$, the state is mapped to zero as soon as two $j_k$
values coincide, so we can eliminate such states. If we thus restrict the set to states
$\ket{j_1,j_2,\ldots ,j_N}$ which are such that the modes are ordered as $j_1 < j_2 < \ldots
< j_N$ (for fermions) or $j_1 \leq j_2 \leq \ldots \leq j_N$ (for bosons), then we have a
linearly independent set of states. For fermions, this implies in particular that we need to
have $N \leq L$, there cannot be more fermions in the system then there are linearly
independent modes (single particle states).

Finally, one can wonder about the normalisation of these states. For fermions, the
superposition created by $\hat{S}^-$ contains $N!$ terms, which are mutually orthogonal, so
that the resulting state is normalised, because of the $1/\sqrt{N!}$ prefactor in the
definition of $\hat{S}^{-}$. More generally, one then finds

```{math}
\braket{i_1 < i_2 < \ldots < i_N | j_1 < j_2 < \ldots < j_N} = \delta_{i_1,j_1}
\delta_{i_2,j_2} \cdots \delta_{i_N,j_N}
```

For bosons, the situation is more complicated in the case that some $j_k$ values coincide.
Some of the $N!$ terms created by $\hat{S}^+$ are then equal and contribute differently to
the norm. If we denote with $n_1, n_2, \ldots, n_L$ the number of $j$ values that equal the
value $1, 2, \ldots, L$, i.e. the number of particles in mode $1, 2, \ldots, L$, then we
find

```{math}
\braket{i_1 \leq i_2 \leq \ldots \leq i_N | j_1 \leq j_2 \leq \ldots \leq j_N} = (n_1! n_2!
\cdots n_L!) \delta_{i_1,j_1} \delta_{i_2,j_2} \cdots \delta_{i_N,j_N}
```

This more general exprression is also valid for fermions, where every $n_j$ is restricted to
be zero or one. In fact, the values $n_j$ for $j=1,\ldots,L$ completely characterise the
state, and can thus be used to relabel the basis. Instead of specifying the mode $j_k$ that
each particle $k=1,\ldots,N$ occupies (where the labeling of the particles is arbitrary
because they are identical), we can move to a mode-based description and thus specify the
number of particles in each mode, also known as the mode occupation number. We can then
refer to the basis vectors as

```{math}
\ket{n_1, n_2, \ldots, n_L}
```

where $n_j = 0, 1$ (fermions) or $n_j = 0,1,2, \ldots $ (bosons) and furthermore
$\sum_{j=1}^{L} n_j = N$. Furthermore, we define these states to be normalised to 1, i.e. we
absorb a suitable normalisation factor when defining $\ket{n_1, n_2, \ldots, n_L}$ in terms
of the construction above.

This way of labelling the basis states now is again reminiscent of a tensor product
structure, i.e. we could think of $\ket{n_1, n_2, \ldots, n_L}$ as the tensor product of
states $\ket{n_j}$ associated to every mode, and where the Hilbert space associated with
such a mode is two-dimensional in the case of fermions, or infinite-dimensional in the case
of bosons. However, there is still a global constraint $\sum_{j=1}^{L} n_j = N$ so that we
cannot let the different $n_j$ values vary completely independently from each other.
Furthermore, some caution is now needed as to what it means to have operators acting on
these different "mode Hilbert spaces". The correct formalism is that of second quantisation,
which we introduce next.

```{note}
In many applications, people do still work with the framework of first quantisation, and
consider $N$-particle states constructed by symmetrising or antisymmetrising the tensor
product of $N$ single-particle states, in a so-called independent particle model or
approximation. Such states are quite cumbersome to work with. As can already be seen, the
antisymmetric case is slightly easier and is known as a Slater determinant. Indeed, the
antisymmetrisation formula is reminiscent of the Leibniz formula of a determinant, and for
example the inner product between two Slater determinants constructed from
$\{\ket{\psi_n},n=1,\ldots,N\}$ and $\{\ket{\varphi_n},n=1,\ldots,N\}$ is given by the
determinant of the matrix containing all overlaps $\braket{\varphi_m \vert \psi_n}$. Slater
determinants form the basis of Hartree-Fock theory for approximating the state of electrons
in an atom or molecule.

The bosonic version occurs in the context of Bose-Einstein condensation and cold atom
systems more generally. In that case, the inner product between two such states gives rise
to a determinant-like formula, but without the minus signs. This construction is known as
the permenant, but unlike the determinant it is very hard to compute in general and really
requires to explicitly sum up all $N!$ terms.
```

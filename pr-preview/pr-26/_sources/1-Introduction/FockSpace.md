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

(fock_space)=
# Fock Space and Second Quantisation

When working with basis vectors using the occuppation number representation, we might
consider dropping the overall constraint $\sum_{j=1}^L n_j = N$. This amounts to working in
a larger Hilbert space, which is known as the Fock space, and consists of the direct sum of
all physical (symmetrised or antisymmetrised) Hilbert spaces $\mathbb{H}^{(N)}$ for
different particle numbers $N$, going all the way from $N=0$:

```{math}
\mathbb{H} = \bigoplus_{N=0}^{+\infty} \mathbb{H}^{(N)}
```

In the case of fermions and with a finite-dimensional single-particle Hilbert space
$\mathbb{H}^{(1)} \cong \mathbb{C}^L$, the upper limit in the direct sum is $N=L$, i.e.
there are no states with $N > L$ and so the associated Hilbert spaces are zero-dimensional.
This direct sum furthermore also contains the case $N=0$, which we have not discussed
before. In the previous subsection we started the construction of $\mathbb{H}^{(N)}$ from a
given single particle Hilbert space $\mathbb{H}^{(1)}$. When there are no particles in the
system, there is only a single state in which it can be, corresponding to having all
occupation numbers $n_j = 0$ for all $j$. Hence, for $N=0$ particles, the Hilbert space
$\mathbb{H}^{(0)}$ is spanned by a single state, which we typically denote as
$\ket{\Omega}=\ket{0,0,\dots,0}$ and refer to as the *vacuum state*. Note that this vacuum
state is normalised, and is thus very different from an actual zero vector of the vector
space, which has norm zero.

The Fock space becomes a Hilbert space simply by incorporating the inner product from each
of its summands. States within the different summands of this direct sum are defined to be
orthogonal, i.e. $\braket{\varphi^{(M)} \vert \psi^{(N)}}=0$ for all $M$-particle states
$\ket{\varphi^{(M)}}$ and $N$-particle states $\ket{\psi^{(N)}}$ with $M \neq N$.

The main benefits of using the formalism of second quantisation are not about losing the
overall particle number constraint, but for working with operators, in particular to
describe (interacting) Hamiltonians. In first quantisation, we need to specify a Hamiltonian
given a particular number of particles, i.e.\ the number of particles is an external
parameter of the system. Using the Fock space, we can now define operators in such a way
that their action is immediately defined for states with an arbitrary number of particles,
including even states which are superpositions over different particle numbers.

Hereto, we first introduce operators that enable us to connect the different particle number
sectors, by creating (adding) or annihilating (removing) particles in the system. In
particular, we denote with $\hat{a}_j^+$ the operator that adds a new particle in the mode
$j$ in the system and with $\hat{a}_j^-$ the operator that removes a particle that is in
mode $j$ from the system. As it turns out that both operators are related via the adjoint,
i.e. $\braket{\Phi| \hat{a}_j^+ \Psi} = \braket{\hat{a}_j^- \Phi | \Psi}$, we use the
simpler notation $\hat{a}_j$ for the *annihilation operator* and $\hat{a}_j^\dagger$ for the
*creation operator*. To construct these operators in a mathematically precise and
constructive way is actually somewhat tedious (but see
[Wikipedia](https://en.wikipedia.org/wiki/Second_quantization#Creation_and_annihilation_operators)).
We just summarize their main properties. In particular, we want to have the property that
the (anti)symmetrized states satisfy

```{math}
\ket{j_1,j_2, \ldots, j_N} = \hat{a}_{j_1}^\dagger \hat{a}_{j_2}^\dagger \cdots
\hat{a}_{j_N}^\dagger\ket{\Omega}.
```

It is immediately clear that, because of the (anti)symmetry, this requires that

```{math}
[\hat{a}_i^\dagger, \hat{a}_j^\dagger] = 0\ (\text{bosons})\quad\text{or}\quad
\{\hat{a}_i^\dagger,\hat{a}_j^\dagger\} = 0\ (\text{fermions}).
```

From the normalisation of these states, it also follows that

```{math}
[\hat{a}_i, \hat{a}_j^\dagger] = \delta_{i,j}\ (\text{bosons})\quad\text{or}\quad
\{\hat{a}_i,\hat{a}_j^\dagger\} = \delta_{i,j}\ (\text{fermions}).
```

With respect to the normalized basis vectors, using the occupation representation, we have

```{math}
\ket{n_1, n_2, \ldots, n_L} = \frac{1}{\sqrt{n_1! n_2! \cdots n_L!}}
(\hat{a}_1^\dagger)^{n_1} (\hat{a}_2^\dagger)^{n_2} \cdots (\hat{a}_L^\dagger)^{n_L}
\ket{\Omega}
```

which can be summarized using

```{math}
\hat{a}_j \ket{n_1, n_2, \ldots, n_j, \ldots, n_L} &= (\pm 1)^{n_1 + n_2 + \ldots + n_{j-1}} \sqrt{n_j} \ket{n_1, n_2, \ldots, n_j - 1, \ldots, n_L},\\
\hat{a}_j^\dagger \ket{n_1, n_2, \ldots, n_j, \ldots, n_L} &= (\pm 1)^{n_1 + n_2 + \ldots + n_{j-1}} \sqrt{n_j+1} \ket{n_1, n_2, \ldots, n_j + 1, \ldots, n_L}.
```

It then follows easily that the operator $\hat{n}_j = \hat{a}_j^\dagger \hat{a}_j$ satisfies

```{math}
\hat{n}_j \ket{n_1, n_2, \ldots, n_j, \ldots, n_L} = n_j \ket{n_1, n_2, \ldots, n_j, \ldots, n_L}
```

and thus measures the number of particles in mode $j$. The operators $\hat{n}_j$ are
referred to as *number operators*. The total number of particles can then be measured using

```{math}
\hat{N} = \sum_{j=1}^{L} \hat{n}_j
```

but the Fock space does of course contain states which are superpositions over different
particle numbers (and which are thus not eigenstates of $\hat{N}$).

Furthermore, by studying how single particle states $\ket{j} \equiv \hat{a}_j^\dagger
\ket{\Omega}$ change under a change of single particle basis, or thus, a transformation to a
new set of modes, we can deduce how the associated creation and annihilation operators
transform. Suppose we have a different single-particle basis, which for clarity we label
with greek letters $\kappa = 1,\ldots,L$. We then find

```{math}
\ket{\kappa} = \hat{a}_\kappa^\dagger \ket{\Omega} = \sum_{j} \ket{j} \braket{j\vert\kappa}
= \sum_{j} \braket{j\vert \kappa} \hat{a}_j^\dagger \ket{\Omega}
```

from which we obtain

```{math}
\hat{a}_\kappa^\dagger = \sum_{j} \braket{j\vert \kappa} \hat{a}_j^\dagger,
\qquad\hat{a}_\kappa = \sum_{j} \braket{\kappa\vert j} \hat{a}_j.
```

Note that the transformation matrix $\braket{j \vert \kappa}$ between two orthonormal bases
correspond to a unitary matrix. These transformation rules will be employed often, for
example, to switch between a position and momentum space representation.

````{note}
The bosonic creation and annihilation operators are of course reminiscent from the operators
introduced for diagonalising the single particle harmonic oscillator model. Indeed, out of
the bosonic creation and annihilation operator associated to every mode $j$ we can build two
Hermitian operators

```{math}
\hat{X}_j = \frac{1}{\sqrt{2}}(\hat{a}_j + \hat{a}_j^\dagger),\quad \hat{P}_j =
\frac{-\mathrm{i}}{\sqrt{2}}(\hat{a}_j - \hat{a}_j^\dagger)
```

which than satisfy the well known commutation relations $\left[\hat{X}_j, \hat{P}_k\right] =
\mathrm{i} \delta_{j,k}$. In second quantisation, the Fock space of bosons built from a
single particle system with $L$ modes can equivalently be thought of as a regular tensor
product space of $L$ distinguishable quantum particles moving on the real line, or
technically, as $\left(L^2(\mathbb{R})\right)^{\otimes L}$.
````
````{note}
For fermions, we can also construct Hermitian operators out of the creation and annihilation
operators, which we denote as

```{math}
\hat{\eta}^{(1)}_j = \frac{1}{\sqrt{2}}(\hat{a}_j + \hat{a}_j^\dagger),\quad
\hat{\eta}^{(2)}_j = \frac{-\mathrm{i}}{\sqrt{2}}(\hat{a}_j - \hat{a}_j^\dagger).
```

In this case, we find that they satisfy the anticommutation relation

```{math}
\{ \hat{\eta}^{(\alpha)}_j, \hat{\eta}^{(\beta)}_k \} = \delta_{\alpha,\beta} \delta_{j,k}
```

so that the $\hat{\eta}^{(1)}$ type operators and $\hat{\eta}^{(2)}$ type operators behave
similarly. In that case, one often uses a different notation by setting

```{math}
\hat{\chi}_{2j-1} = \hat{\eta}^{(1)}_j = \frac{1}{\sqrt{2}}(\hat{a}_j +
\hat{a}_j^\dagger),\quad \hat{\chi}_{2j} = \hat{\eta}^{(2)}_j =
\frac{-\mathrm{i}}{\sqrt{2}}(\hat{a}_j - \hat{a}_j^\dagger)
```

and thus $\{\hat{\chi}_k, \hat{\chi}_l\} = \delta_{k,l}$ for all $k, l = 1,\ldots, 2L$.
These Hermitian fermionic operators are referred to as **Majorana operators**.

Note furthermore that the Fock space of fermions built from a single particle system with
$L$ modes looks remarkably like a system of $L$ qubits, i.e. the tensor product
$(\mathbb{C}^2)^{\otimes L}$. While this is true for how the occupation number basis vectors
are labelled, one important fact is that the operators $\hat{a}_j$ and $\hat{a}_j^\dagger$
should not be thought of as local operators that act nontrivially on the single site $j$,
and as the identity elsewhere, since they do not mutually commute, but rather anticommute.
It is possible to map these fermionic creation and annihilation operators to 'nonlocal'
qubit operators using the
[Jordan-Wigner transformation](https://en.wikipedia.org/wiki/Jordan–Wigner_transformation).
````

With these creation and annihilation operators, we can now represent general operators in a
way that does not depend on the precise number of particles in the system. The simplest case
are 'single-particle' operators, i.e. operators that were defined with respect to the
single-particle Hilbert space $\mathbb{H}^{(1)}$. The easiest case are operators which are
are diagonal with respect to the chosen single-particle basis. In that case, every particle
one of the eigenmodes of the single-particle operator will give a contribution that equals
the associated eigenvalue. Hence, the many-body representation of such an operator is given
by

```{math}
\hat{O}^{(1)} = \sum_{j} \lambda_j \ket{j}\bra{j} \quad \rightarrow \quad \hat{O} = \sum_{j}
\lambda_j \hat{a}_j^\dagger \hat{a}_j.
```

However, we can easily transform away from the basis of eigenmodes to a general set of
modes, and then find

```{math}
\hat{O}^{(1)} = \sum_{j,k} O_{j,k} \ket{j}\bra{k} \quad \rightarrow \quad \hat{O} =
\sum_{j,k} O_{j,k} \hat{a}_j^\dagger \hat{a}_k.
```

Vice versa, if you are given an operator that only contains terms where every term contains
exactly one creation and one annilation operator, then it is especially easy to diagonalise
this operator, since one only needs to diagonalise the corresponding single-particle version
of the operator. When the Hamiltonian of the many-body system is of this form, the system is
said to be *free* or noninteracting.

```{note}
There is a larger class of operators that can easily be diagonalised, namely operators where
every term is quadratic in the creation and annihilation operators. This means that every
term contains either a creation and an annilation operator, or two creation operators, or
two annihilation operators. Such Hamiltonians are said to be quadratic or Gaussian, and can
be diagonalised using a
[Bogoliubov transformation](https://en.wikipedia.org/wiki/Bogoliubov_transformation).
```

Similarly, there exist two-particle operators, in particular, typical interaction terms in
the Hamiltonian such as the Coulomb interaction between electrons. Such operators take the
form

```{math}
\hat{O}^{(2)} = \sum_{j,k,l,m} O_{j,k; l,m} \ket{j,k} \bra{l,m}
```

and can be translated to act on the full Fock space as

```{math}
\hat{O} = \sum_{j\leq k; l\leq m} O_{(j,k); (l,m)} \hat{a}_{j}^\dagger \hat{a}_k^\dagger
\hat{a}_m \hat{a}_l = \frac{1}{4} \sum_{j, k; l, m} O_{(j,k); (l,m)} \hat{a}_{j}^\dagger
\hat{a}_k^\dagger \hat{a}_m \hat{a}_l.
```

As soon as such type of operators are present in the Hamiltonian (which thus contain more
than two creation of annihilation operators), it becomes impossible to diagonalise the
Hamiltonian based on a simple calculation in the single-particle Hilbert space, and the
exponentially large many-body Hilbert space need to be considered.

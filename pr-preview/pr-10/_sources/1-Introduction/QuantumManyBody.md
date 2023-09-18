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

 # Quantum Many-Body Physics

## Quantum mechanics and its postulates

While the energy levels of the hydrogen atom played an important role in the historical
development of quantum mechanics, it became almost immediately clear that the true challenge
is in applying the laws of quantum mechanics to systems with many interacting particles
or fields. Note that the formalism of quantum mechanics, and in particular its postulates,
are generically valid and not restricted to the description of a single particle. These
postulates characterise the mathematical model by which quantum mechanics describes physical
systems, and more specifically how it represents states, observables, measurements and
dynamics. We briefly reiterate this postulates and base our discussion on the wonderfull
lecture notes ``Quantum Information and Computation'' by John Preskill. 

###	**Postulate 1: States**

The state of an isolated quantum system is associated to a ray of vectors in a complex
Hilbert space $\mathbb{H}$.

A Hilbert space is a metric complete inner product space. Let us unpack this definition:
1.  $\mathbb{H}$ is a vector space in this case over the complex numbers. We will denote
    elements of this vector space with Dirac's ket notation $\ket{\psi}$.
    In particular, we can build linear combinations

    ```{math}
    \ket{\psi} = a \ket{\psi_1} + b \ket{\psi_2}
    ```

    for all $a, b \in \mathbb{C}$ and all $\ket{\psi_1}, \ket{\psi_2} \in \mathbb{H}$.

2.  $\mathbb{H}$ has an inner product, which maps to vecors $\ket{\psi}$ and
    $\ket{\varphi}$ onto a scalar $\braket{\varphi|\psi} \in \mathbb{C}$
    with the properties that
    * Linearity: $\bra{\varphi} ( a \ket{\psi_1} + b \ket{\psi_2}) = a \braket{ \varphi | \psi_1} + b \braket{ \varphi | \psi_2}$
    * Skew-symmetry: $\braket{ \varphi | \psi} = \braket{ \psi | \varphi}^\ast$
	* Positivity: $\braket{ \psi | \psi} \geq 0$ with equality only if $\ket{\psi} = 0$.
	
	This last property enables us to define a norm $\lVert \psi \rVert = \lVert \ket{\psi} \rVert = \sqrt{\braket{\psi|\psi}}$,
	which satisfies known properties such as $\lVert \psi \rVert = 0 \Leftrightarrow \ket{\psi} = 0$
	$\lVert a \ket{\psi} \rVert = a \lVert \psi \rVert$ and the triangle inequality
	$\lVert \ket{\varphi} + \ket{\psi} \rVert \leq \lVert \varphi \rVert + \lVert \psi \rVert$.

3.	The final property of metric completeness is a technical requirement that is only relevant
	in infinite-dimensional Hilbert spaces. Firstly, a metric is a notation of distance between
	the elements in $\mathbb{H}$, which is provided by the norm of the difference, i.e.
	$d(\varphi, \psi) = \lVert \varphi - \psi \rVert.
	
	Completeness of the metric is a specific property that guarantees that certain sequences
	of vectors are guaranteed to have a limit value that also exists in $\mathbb{H}$. This is
	necessary to make sense of e.g. Fourier series.

The state of a quantum system is associated to a ray of vectors, which is the one-dimensional
space $\{ a \ket{\psi} , \forall a \in \mathbb{C}\}$ spanned by a single (nonzero)
vector $\ket{\psi} \in \mathbb{H}$. We will describe the state of the system using a
single representative $\ket{\psi}$ of this ray, which we typically choose such that
$\braket{ \psi | \psi} = 1$. However, this does not fix the representative completely,
as we can still add arbitrary phases  $\exp(\mathrm{i}\alpha)$, i.e. $\ket{\psi}$ and
$\mathrm{e}^{\mathrm{i}\alpha} \ket{\psi}$ describe the same state.

The best known Hilbert space from your courses on single-particle quantum mechanics is probably
$L^2(\mathbb{R}^n)$, the Hilbert space for a single quantum particle moving in the $n$-dimensional
coordinate space $\mathbb{R}^n$ (typically $n=1,2,3$). This Hilbert space corresponds to the
space of all square-integrable functions $\psi:\mathbb{R}^d \to \mathbb{C}: x \mapsto \psi(x)$
and the inner product is given by

$$\braket{\varphi | \psi} = \int_{\mathbb{R}^n} \varphi(x)^\ast \psi(x)\,\mathrm{d} x$$

However, this is already a complicated Hilbert space from a technical perspective. Hilbert
spaces can also be finite-dimensional, i.e. $\mathbb{C}^d$, the space of column vectors of
length $d$, with the standard Euclidean inner product

$$\braket{\varphi | \psi} = \sum_{i=1}^d \varphi_i^\ast \psi_i$$

These Hilbert spaces will be very important in our discussion. The simplest nontrivial case
corresponds to $d=2$ and the associated quantum system is known under various names. It is
often referred to as a qubit in the context of quantum information theory. There are various
ways in which qubits can be physically realised. Another common example of a two-dimensional
Hilbert space is for describing the spin degree of freedom of an electron, or another particle
with spin quantum number 1/2. Such a reduced description (forgetting about the position) is
possible if the electron is localised in space, for example when it is strongly bound to an atom.

If we do want to describe a particle that moves in space, we might also consider it to exist
only at discrete positions in space, i.e. on a lattice. For example, on a one-dimensional 
lattice (a.k.a. a chain) with $L$ sites, the Hilbert space would also correspond to
$\mathbb{H} = \mathbb{C}^L$ and the standard basis vectors $\vert j \rangle, j=1,\dots,L$ 
correspond to the state of the system if the particle is exactly localised on site $j$. We
can also consider infinitely large lattices, e.g. the one-dimensional chain where there is
a site associated with every $j \in \mathbb{Z}$ (or the n-dimensional hypercubic lattice $\mathbb{Z}^n$).
The resulting Hilbert space is then spanned by the states $\vert j \rangle$ for all
$j \in \mathbb{Z}$, and is thus infinite-dimensional but with a straightforward
countably infinite basis.

Of course, our goal is to find the Hilbert space of a many body system. We return to this
question below and devote a complete section to it.

###	**Postulate 2: Observables**

Physical observables of the system correspond to self-adjoint (a.k.a. Hermitian) linear
operators on the Hilbert space $\mathbb{H}$.

An operator $\hat{A}$ on $\mathbb{H}$ is a linear map $\hat{A}:\mathbb{H} \to \mathbb{H}$,
i.e. a map from vectors to vectors that satisfies

$$ \hat{A}( a \ket{\varphi} + b \lvert \psi \rangle) = a \hat{A}(\ket{\varphi}) + b \hat{A}(\ket{\psi}) $$

The adjoint of an operator $\hat{A}$ is a new operator $\hat{A}^\dagger$ that is constructed
such that

$$ \bra{\varphi} \hat{A} \psi \rangle = \langle \hat{A}^\dagger \varphi \ket{\psi} $$

for all $\ket{\varphi}, \ket{\psi} \in \mathbb{H}$ and where
$\vert \hat{A}\psi \rangle = \hat{A} \ket{\psi}$. This definition requires that
$(a_1 \hat{A}_1 + a_2 \hat{A}_2)^\dagger = a_1^\ast \hat{A}_1^\dagger + a_2^\ast \hat{A}_2^\dagger$
and $(\hat{A}_1 \hat{A}_2)^\dagger = \hat{A}_2^\dagger \hat{A}_1^\dagger$.

A self-adjoint operator is an operator such that $\hat{A}^\dagger = \hat{A}$ or thus

$$ \braket{\varphi | \hat{A} \psi } = \braket{ \hat{A}^\dagger \varphi | \psi} $$

for all $\ket{\varphi}, \ket{\psi} \in \mathbb{H}$. Linear combinations of
self-adjoint operators with real coefficients are self-adjoint. The composition of two self-adjoint
linear operators $\hat{A}_1 \hat{A}_2$ is self-adjoint if and only if

$$ \left[ \hat{A}_1 , \hat{A}_2 \right] = \hat{A}_1 \hat{A}_2 - \hat{A}_2 \hat{A}_1 = 0, $$

i.e. if the operators also commute. Self-adjoint operators have real eigenvalues, and eigenvectors 
associated to distinct  eigenvalues are orthogonal. In a finite-dimensional Hilbert space,
all self-adjoint operators are represented as Hermitian matrices and admit a spectral decomposition

$$\hat{A} = \sum_{n} \lambda_n \hat{P}_n$$

where $\hat{P}_n$ is the spectral projector onto the eigenspace associated with $\lambda_n$.
The spectral projectors satisfy $\hat{P}_n \hat{P}_m = \delta_{n,m}$, $\hat{P}_n^\dagger = \hat{P}_n$
and $\sum_{n} \hat{P}_n = \mathbb{1}$, the identity operator. If $\lambda_n$ has one-dimensional
eigenspace spanned by the eigenvector $\vert\phi_n\rangle$, then

$$\hat{P}_n = \frac{\vert \phi_n \rangle \langle \phi_n \vert}{\langle \phi_n \vert \phi_n \rangle}$$

where the denominator can be omitted if the eigenvector is normalised.

### **Postulate 3: Measurements, expectation values and collapse**

TBC

### **Postulate 4: Dynamics**

TBC

## The Hilbert space of many body physics

All of the previous axioms remain valid for a composite system consisting of several quantum
degrees of freedom. However, we need to know how to describe the state of the system, and 
thus more specifically, how to define the Hilbert space associated to such a system. It
turns out that quantum mechanics forces us to distinguish two cases.

### Distinguisable particles and tensor products

Consider a quantum system composed out of two subsystems, which we call $A$ and $B$, sometimes
referred to as Alice and Bob in quantum information context. These can themselves already
be many-body systems. Suppose we know the Hilbert space $\mathbb{H}^A$ in which to describe
states of subsystem $A$ when considered as an isolated system on itself, and analoguously
for $\mathbb{H}^B$. Now consider both systems together, but where they do not interact,
so that we can still treat them independently. In particular, we can prepare subsystem $A$
in a state $\ket{\psi^A}$ and subsystem $B$ in a state $\ket{\psi^B}$. We should also be
able to describe these two independent subsystems jointly, so that there must exist a map
from the two arguments $(\ket{\psi^A}, \ket{\varphi^B}) \in \mathbb{H}^A \times \mathbb{H}^B$
to a single state which we denote as $\ket{ \psi^A} \otimes \ket{\varphi^B}$ and that lives
in a joint Hilbert space $\mathbb{H}^{AB}$ that we have yet to determine.

Now, it makes sense that, if we build superpositions in one of the two subsystems, while keeping
the other fixed, this also correspond to a superposition in the joint description of both
systems together. This leads to

$$ \left (a_1 \ket{\psi^A_1} + a_2 \ket{\psi^A_2}\right) \otimes \ket{\varphi^B} = a_1 \ket{\psi^A_1 }\otimes \ket{\varphi^B} + a_2 \ket{\psi^A_1 }\otimes \ket{\varphi^B}$$

and similarly

$$ \ket{\psi^A} \otimes \left(b_1 \ket{\varphi^B_1} + b_2 \ket{\varphi^A_2}\right) = b_1 \ket{\psi^A }\otimes \ket{\varphi^B_1} + b_2 \ket{\psi^A }\otimes \ket{\varphi^B_2}.$$

Hence, the Hilbert space $\mathbb{H}^{AB}$ that we are trying to construct must contain
all states $\ket{\psi^A} \otimes \ket{\varphi^B}$ for all $\ket{\psi^A} \in \mathbb{H}^A$ and
all $\ket{\varphi^B} \in \mathbb{H}^B$, all possible linear combinations thereof (in order to
be a vector space), but in such a way that the above equalities hold. This construction, which
can be made mathematically precise, is known as the tensor product of vector spaces
$\mathbb{H}^{AB} = \mathbb{H}^A \otimes \mathbb{H}^B$. 

We have also denoted the output of the map from two states $(\ket{\psi^A}, \ket{\varphi^B}) \in \mathbb{H}^A \times \mathbb{H}^B$
to $\mathbb{H}^A \otimes \mathbb{H}^B$ using the same tensor product symbol, and refer to
such a state as a (tensor) product state $\ket{ \psi^A} \otimes \ket{\varphi^B}$. Importantly,
however, the tensor product space $\mathbb{H}^A \otimes \mathbb{H}^B$ certainly contains vectors
which are not product states, such as

$$ a_1 \ket{ \psi_1^A} \otimes \ket{\varphi_1^B} + a_2 \ket{ \psi_2^A} \otimes \ket{\varphi_2^B}. $$

This forms the basis for quantum correlations and the concept of (quantum) entanglement, which
will be a fundamental property of quantum many-body systems. That the Hilbert space of a
composite system is given by the tensor product of the individual Hilbert spaces is often
introduced as a separate axiom. The deductive (but informal) argument just given can however
be turned into a proof that depends only on the axioms given above (in fact only on the first two).

As expected (and required), it can be shown that the tensor product of two Hilbert spaces is 
again a Hilbert space, if we define its inner inner product in the following way.
We first define the inner product for product states as

$$ \braket{\psi_1^A \otimes \varphi_1^B | \psi_2^A \otimes \varphi_2^B } = \braket{\psi_1^A | \psi_2^A} \braket{\varphi_1^B | \varphi_2^B} $$

and then extend this definition by linearity (in the second argument and antilinearity in the first argument).

In practice, given two finite-dimensional Hilbert spaces $\mathbb{H}^A \cong \mathbb{C}^{d_A}$
and $\mathbb{H}^B \cong \mathbb{C}^{d_B}$ with a basis $\{ \ket{j}, j=1,\dots, d_A\}$
and $\{\ket{k}, k=1,\dots, d_B\}$, the tensor product space is spanned by a basis composed
of all products

$$\{ \ket{j,k} = \ket{j} \otimes \ket{k}, j=1,\dots, d_A, k=1,\dots, d_B\}$$

and thus has dimension $d_A \dot d_B$. A general state $\ket{\Psi} \in \mathbb{H}^{A}\otimes \mathbb{H}^B$
can then be expanded as

$$ \ket{\Psi} = \sum_{j=1}^{d_A }\sum_{k=1}^{d_B} \Psi_{jk} \ket{j,k} $$

The expansion coefficients $\Psi_{jk}$ thus have two indices, and it is often useful to think
of them as a matrix. Note that we will almost always use this product basis, also sometimes
referred to as the computational basis, for working with tensor product spaces. However, one
can certainly also use more complicated basis choices, where the basis vectors are not simple
product states. One well known choice that you might remember from your quantum mechanics course
is in the case of two spin-1/2 systems. If we denote the basis for a single spin-1/2 system as
$\{\ket{\uparrow},\ket{\downarrow}\}$, then the product basis for a system consisting of two
spin-1/2 systems is given by $\{\ket{\uparrow,\uparrow}, \ket{\downarrow,\uparrow}, \ket{\uparrow,\downarrow}, \ket{\downarrow,\downarrow}, \}$.
However, in the context of spin coupling (see Section on Symmetries), one also uses the coupled
basis

```{math}
\ket{0,0} = \frac{1}{\sqrt{2}} \left(\ket{\uparrow,\downarrow} - \ket{\downarrow,\uparrow}\right)\\
\ket{1,1} = \ket{\uparrow,\uparrow}\\
\ket{1,0} = \frac{1}{\sqrt{2}} \left(\ket{\uparrow,\downarrow} + \ket{\downarrow,\uparrow}\right)\\
\ket{1,-1} = \ket{\downarrow,\downarrow}
```

The tensor product construction extends to systems with multiple subsystems. Consider for example
a system consisting of qubits, where every individual qubit has an associated Hilbert
space $\mathbb{C}^2$ with basis denoted as $\{\ket{0},\ket{1}\}$. The Hilbert space $\mathbb{H}^N$ of $N$
qubits is then spanned by a computational basis which we can denote as

$$\{\ket{s_1, s_2, \ldots, s_N} = \ket{s_1} \otimes \ket{s_2} \otimes \cdots \otimes \ket{s_N}; s_1 =0,1; s_2 =0,1; \ldots; s_n =0,1\}.$$

Hence, the Hilbert space thus has dimension $2^N$, and a general state $\ket{\Psi}$ has expansion coefficients

$$\Psi_{s_1,s_2, \ldots, s_N}$$

which can be interpreted as a single vector of length $2^n$, or as a $n$-dimensional tensor,
where every tensor index ranges over the two values 0 and 1.

### Identical particles and Pauli's exclusion principle

The tensor product construction needs to be revised when discussing the Hilbert space of a
system composed of identical particles. Consider for example a system made out of $N$ 
identical particles. To every individual particle we can associate a particular Hilbert
space, which we denote as $\mathbb{H}^{(1)}$, for example $\mathbb{H}^{(1)} = L^2(\mathbb{R})$
for a particle moving on the real line, or $\mathbb{H}^{(1)} = \mathbb{C}^L$ for a particle living
on the sites of a chain of length $L$.

If we temporarily assign each of the $N$ particles a label $n=1,\dots, N$, then the Hilbert
space of the composite system would be given by the $N$-fold tensor product
$\tilde{\mathbb{H}}^{(N)} = \left(\mathbb{H}^{(1)}\right)^{\otimes N}$. However, for identical particles,
our labeling is completely arbitrary. For the case of $N=2$ particles on a chain of $L$ sites,
we cannot distinguish between the state $\ket{j_1, j_2}$ where particle $1$ is on site $j_1$
and particle $2$ is on site $j_2$ versus the state $\ket{j_2, j_1}$ where site $j_1$ is occupied
by the particle that we gave label $2$ and site $j_2$ is occupied by the particle with label $1$.
A general redefinition of the particle labels amounts to a permutation, and we have to require
that no physical measurement can distinguish between such permutations. Hence, this permutation
invariance does not behave like a regular symmetry (like e.g. rotation symmetry, one can still construct 
observables along preferred directions such that they can detect rotations).

We are forced to restrict our tensor product Hilbert space $\left(\mathbb{H}^{(1)}\right)^{\otimes N}$
to the subspace $\mathbb{H}^{(N)}$ of physical states which are not affected by acting with
such permutations. Note that, due to the fact that quantum states actually correspond to rays
of vectors, it is still allowed that the vectors in $\mathbb{H}^{(N)}$ pick up a phase factor
when applying certain permutations. It is a result in the representation theory of the
permutation group that there are only two possibilities. Either the phase factor is always 
absent (or thus 1), or the phase factor is (-1) for odd permutations and (+1) for even
permutations, i.e. the phase factor equals the sign(ature) of the permutation. Identical particles
for which the phase factor is always one are known as *bosons*, whereas those with the nontrival
phase factor choice correspond to *fermions*. Indeed, the nontrivial phase factor automatically
gives rise to Pauli's exclusion principle: two fermions cannot be in the same quantum state,
since $P_{12} \ket{j_1,j_2} = \ket{j_2,j_1} = -\ket{j_1,j_2}$ and for $j_1=j_2$ we would thus
find $\ket{j,j} = -\ket{j,j}$.

TBC



## Fock space and second quantisation

TBC
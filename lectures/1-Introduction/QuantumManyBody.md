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
is in applying the laws of quantum mechanics to systems with many interacting particles or
fields. Note that the formalism of quantum mechanics, and in particular its postulates, are
generically valid and not restricted to the description of a single particle. Quantum field
theory also follows these postulates and is thus not a generalisation of quantum mechanics,
but rather a specific case of it. These postulates characterise the mathematical model by
which quantum mechanics describes physical systems, and more specifically how it represents
states, observables, measurements and dynamics. We briefly reiterate this postulates and
base our discussion on the wonderfull lecture notes "Quantum Information and Computation" by
John Preskill.

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
	$\lVert a \ket{\psi} \rVert = \vert a\vert  \lVert \psi \rVert$ and the triangle inequality
	$\lVert \ket{\varphi} + \ket{\psi} \rVert \leq \lVert \varphi \rVert + \lVert \psi \rVert$.

3.	The final property of metric completeness is a technical requirement that is only relevant
	in infinite-dimensional Hilbert spaces. Firstly, a metric is a notation of distance between
	the elements in $\mathbb{H}$, which is provided by the norm of the difference, i.e.
	$d(\varphi, \psi) = \lVert \varphi - \psi \rVert$.
	
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
self-adjoint operators admit a spectral decomposition

$$\hat{A} = \sum_{n} \lambda_n \hat{P}_n$$

where $\hat{P}_n$ is the spectral projector onto the eigenspace associated with $\lambda_n$.
The spectral projectors satisfy $\hat{P}_n \hat{P}_m = \delta_{n,m} \hat{P}_n$, $\hat{P}_n^\dagger = \hat{P}_n$
and $\sum_{n} \hat{P}_n = \mathbb{1}$, the identity operator. If $\lambda_n$ has one-dimensional
eigenspace spanned by the eigenvector $\vert\phi_n\rangle$, then

$$\hat{P}_n = \frac{\vert \phi_n \rangle \langle \phi_n \vert}{\langle \phi_n \vert \phi_n \rangle}$$

where the denominator can be omitted if the eigenvector is normalised.

In the language of matrices, these properties can be rephrased as follows: With respect to an
orthonormal basis choice, self-adjoint operators are represented as hermitian matrices. Such
matrices can be diagonalised by a unitary transformation, or thus, we can construct a complete
basis consisting of eigenvectors. With respect to this basis, the self-adjoint operator is
represented by a diagonal matrix with real values on the diagonal.

### **Postulate 3: Measurements, expectation values and collapse**

Given an observable to which we associate the operator $\hat{A}$, we now need to prescribe
the result of measuring this observable with respect to a system that is in a state $\ket{\psi}$.
The most compact way of describing the result is by stating that, the expectation value $\braket{\hat{A}}$
(= the  mean value of the measurement when averaging over an ensemble of identical copies of the system)
is given by

$$ \braket{\hat{A}} = \frac{\braket{\psi \vert \hat{A} \vert \psi}}{\braket{\psi \vert \psi}} $$

By exploiting the fact that this also prescribes the expectation value of all higher moments
$\braket{\hat{A}^k}$, this determines the full probability distribution of the measurement
outcome, and yields the more familiar result: The only possible measurement outcomes are 
given by the eigenvalues $\lambda_n$ of $\hat{A}$, and for a system in state $\ket{\psi}$ 
(now assumed normalized), the probability of obtaining $\lambda_n$ is given by
$p_n = \braket{\psi \vert \hat{P}_n \vert \psi}$ with $\hat{P}_n$ the spectral projector from
above. In the case that $\lambda_n$ has a single (linearly independent) eigenvector
$\ket{\phi_n}$ (also assumed normalised), this amounts to $p_n = \vert \braket{\phi_n|\psi}\vert^2$.

There is a second part to the measurement postulate, which states that, if the measurement
is immediately repeated (without intermediate dynamcis, as described by the next postulate),
then the same measurement outcome is obtained. Because the measurement outcome with respect
to the initial state $\ket{\psi}$ is probabilistic and can yield different results, this
requires that after the first measurement, the state changes is changed. This is the
well-known **collapse** of the wave function. More specifically, if a measurement of
observable $\hat{A}$ is performed in a system with state $\ket{\psi}$ and the measurement
value $\lambda_n$ is obtained, then the state of the system changes to

$$ \ket{\psi} \longrightarrow \frac{\hat{P}_n \ket{\psi}}{\lVert \hat{P}_n \ket{\psi}\rVert}. $$

Note that the denominator cannot vanish, as otherwise the probability of having obtained
measurement outcome $\lambda_n$ would have been zero in the first place.

### **Postulate 4: Dynamics**

During time intervals without measurements, the state of an isolated quantum system evolves 
unitarily according to the (first order linear) differential equation

$$\frac{\mathrm{d}\ }{\mathrm{d} t} \ket{\psi(t)} = - \mathrm{i} \hat{H}(t) \ket{\psi(t)}$$

known as the Schr\"{o}dinger equation, where $\hat{H}(t)$ is the Hamiltonian of the system,
which may itself be time-dependent. In the case of a time-independent Hamiltonian, we can
define the evolution operator

$$ U(t, t') = \exp\left(-\mathrm{i}(t-t')\hat{H}) $$

which relates states at different times via $\ket{\psi(t)} = \hat{U}(t, t') \ket{\psi(t')}$
and is clearly a unitary operator. Clearly, we need to know the Hamiltonian of a system in
order to even start thinking about modelling its quantum properties. We will always assume
that the Hamiltonian is given. In practice, however, the situation can be much more
complicated. Typically, we want to build only an effective quantum description of the system
(e.g. only the electrons, only certain electrons, $\ldots$) and not start all the way down
at the level of fundamental particles and the standard model (which is also only an
effective model valid up to some energy scale).

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

and thus has dimension $d_A \cdot d_B$. A general state $\ket{\Psi} \in \mathbb{H}^{A}\otimes \mathbb{H}^B$
can then be expanded as

$$ \ket{\Psi} = \sum_{j=1}^{d_A }\sum_{k=1}^{d_B} \Psi_{jk} \ket{j,k} $$

The expansion coefficients $\Psi_{jk}$ thus have two indices, and it is often useful to
think of them as a matrix. Note that we will almost always use this product basis, also
sometimes referred to as the computational basis, for working with tensor product spaces.
However, one can certainly also use more complicated basis choices, where the basis vectors
are not simple product states. One well known choice that you might remember from your
quantum mechanics course is in the case of two spin-1/2 systems. If we denote the basis for
a single spin-1/2 system as $\{\ket{\uparrow},\ket{\downarrow}\}$, then the product basis
for a system consisting of two spin-1/2 systems is given by 
$\{\ket{\uparrow,\uparrow}, \ket{\downarrow,\uparrow}, \ket{\uparrow,\downarrow}, \ket{\downarrow,\downarrow}, \}$.
However, in the context of spin coupling (see Section on Symmetries), one also uses the
coupled basis

```{math}
\ket{0,0} = \frac{1}{\sqrt{2}} \left(\ket{\uparrow,\downarrow} - \ket{\downarrow,\uparrow}\right)\\
\ket{1,1} = \ket{\uparrow,\uparrow}\\
\ket{1,0} = \frac{1}{\sqrt{2}} \left(\ket{\uparrow,\downarrow} + \ket{\downarrow,\uparrow}\right)\\
\ket{1,-1} = \ket{\downarrow,\downarrow}
```

Note that we also use teh same tensor product notation as an operation to map operators from
the subsystems into operators acting on the full tensor product Hilbert space. In
particular, the process of measuring operator $\hat{A}$ in subsystem $A$ and simultaneously
operator $\hat{B}$ in subsystem $B$ is associated with an operator $\hat{A}\otimes \hat{B}$
acting on $\mathbb{H}^A \otimes \mathbb{H}^B$, the action of which is first defined on the
product states as

$$\left(\hat{A} \otimes \hat{B}\right) \left(\ket{\psi^A}\otimes \ket{\varphi^B}\right) = \left(\hat{A}\ket{\psi^A}\right) \otimes \left(\hat{B}\ket{\varphi^B}\right) $$

and then extended by linearity. With respect to a product basis, the matrix representation
of $\left(\hat{A} \otimes \hat{B}\right)$ is given by the
[Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product).

When we are only interested in an operator $\hat{O}$ acting on subsystem $A$ without doing
anything on subsystem $B$, we should create the operator $\hat{O} \otimes \hat{1}_B$, with
$\hat{1}_B$ the identity operator of the Hilbert space $\mathbb{H}^B$. Often, we will omit
this explicit tensor product with the identity operator, and simply use some notation which
indicates that an operator acts on a certain subsystem, such as
$\hat{O}^{(A)} = \hat{O} \otimes \hat{1}_B$.

The tensor product construction extends to systems with multiple subsystems. Consider for
example a system consisting of qubits, where every individual qubit has an associated
Hilbert space $\mathbb{C}^2$ with basis denoted as $\{\ket{0},\ket{1}\}$. The Hilbert space
$\mathbb{H}^N$ of $N$ qubits is then spanned by a computational basis which we can denote as

$$\{\ket{s_1, s_2, \ldots, s_N} = \ket{s_1} \otimes \ket{s_2} \otimes \cdots \otimes \ket{s_N}; s_1 =0,1; s_2 =0,1; \ldots; s_n =0,1\}.$$

Hence, the Hilbert space thus has dimension $2^N$, and a general state $\ket{\Psi}$ has
expansion coefficients

$$\Psi_{s_1,s_2, \ldots, s_N}$$

which can be interpreted as a single vector of length $2^N$, or as a $N$-dimensional tensor,
where every tensor index ranges over the two values 0 and 1. This exponential increase of
the Hilbert space dimension with the number of particles is exactly why the quantum
many-body problem is so difficult, but also essential for providing a quantum computer with
its speed-up.

The Hamiltonian of a many-body system typically takes the form of a sum of terms, where
every individual term acts nontrivially on only a few subsystems. One important example that
will reappear throughout these tutorials is the "Quantum Ising Model with transverse
magnetic field", which acts on a system composed of qubits or spin-1/2 particles, and is
defined as

$$ \hat{H} = - J \sum_{\langle i, j \rangle} \sigma^z_i \otimes \sigma^z_j - h \sum_i \sigma^x_i $$

Here, the summation variables $i$ and $j$ correspond to the sites of a lattice. The notation
$\sum_{\langle i,j \rangle}$ denotes a sum over pairs of neighbouring lattice sites $i$ and
$j$. The second sum contains terms $\sigma^x_i$ which act nontrivially only on the site $i$,
and as the identity operator elsewhere. If, for example, we enumerate the sites from $1$ to
$N$, it would act as

$$ \sigma^x_i = \underbrace{1 \otimes 1 \otimes \ldots \otimes 1}_{\text{$i-1$ factors}} \otimes \sigma^x \otimes \underbrace{1 \otimes \ldots \otimes 1}_{\text{$N-i-1$ factors}}$$

with $\sigma^x = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}$ the Pauli x matrix, and $1$
the $2 x 2$ unit matrix. The first set of terms in $\hat{H}$ acts nontrivially on two sites,
and is defined analoguously, using the Pauli z matrices
$\sigma^z = \begin{bmatrix} 1 & 0 \\ 0 & -1 \end{bmatrix}$.


### Identical particles and Pauli's exclusion principle

The tensor product construction needs to be revised when discussing the Hilbert space of a
system composed of identical particles. Consider for example a system made out of $N$ 
identical particles. To every individual particle we can associate a particular Hilbert
space, which we denote as $\mathbb{H}^{(1)}$, for example $\mathbb{H}^{(1)} = L^2(\mathbb{R})$
for a particle moving on the real line, or $\mathbb{H}^{(1)} = \mathbb{C}^L$ for a particle living
on the sites of a chain of length $L$.

If we temporarily assign each of the $N$ particles a label $n=1,\dots, N$, then the Hilbert
space of the composite system would be given by the $N$-fold tensor product
$\widetilde{\mathbb{H}}^{(N)} = \left(\mathbb{H}^{(1)}\right)^{\otimes N}$. However, for identical particles,
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
gives rise to *Pauli's exclusion principle*: two fermions cannot be in the same quantum state,
since $P_{12} \ket{j_1,j_2} = \ket{j_2,j_1} = -\ket{j_1,j_2}$ and for $j_1=j_2$ we would thus
find $\ket{j,j} = -\ket{j,j}$.

Bosons are thus described by states which are symmetric under permutations, whereas fermions
are described by states which are called antisymmetric. We can define an operator on
$\tilde{\mathbb{H}}^{(N)} = $\left(\mathbb{H}^{(1)}\right)^{\otimes N}$ that maps any given
state onto such a (anti)symmeric state, namely by first defining its action on product states
as

$$ \hat{S}^{\pm} \ket{\psi_1} \otimes \ket{\psi_2} \otimes \cdots \otimes \ket{\psi_N} = \frac{1}{\sqrt{N!}} \sum_{\sigma \in S_N} \epsilon_\sigma \ket{\psi_{\sigma(1)}} \otimes \ket{\psi_{\sigma(2)}} \otimes \cdots \otimes \ket{\psi_{\sigma(N)}} $$

and then extending it by linearity. Here, $S_N$ is the symmetric group containing all permutations
$\sigma$ of $N$ elements, where the permutation $\sigma$ is a bijective map from integers $j \in \{1,\dots,N\}$
to a new number $\sigma(j) \in \{1,\dots,N\}$. The sign(ature) $\epsilon_\sigma$ of the permutation
takes the value $+1$ or $-1$, depending on whether the permutation $\sigma$ can be obtained
by composing an even or odd number of elementary transpositions. An elementary transposition 
$\tau_{i,j}$ is a permutation which only interchanges the two numbers $i$ and $j \neq i$:

$$\tau_{i,j}(i) =j, \tau_{i,j}(j) = i, \tau_{i,j}(k) =k, \forall k\neq i \land k \neq j $$

Note that $\hat{S}^{\pm}$ does not necessarily yield a normalised state, and can indeed even
map a state to zero, in order to give rise to Pauli's exclusion principle: $\hat{S}^- \ket{j,j} = 0$.
The image of $\hat{S}^{\pm}$ contains all states with the proper behaviour under relabeling
permutations, and thus correspond to the physical Hilbert space for bosons or fermions:

$$ \mathbb{H}^{(N)} = \hat{S}^{\pm} \widetilde{\mathbb{H}}^{(N)} = \hat{S}^{\pm} \left(\mathbb{H}^{(1)}\right)^{\otimes N} $$

Note that in this case, the physical Hilbert space is not a tensor product. We had to define
a larger auxiliary Hilbert space, $ \widetilde{\mathbb{H}}^{(N)}$, which is a tensor product.
The physical Hilbert space can then be thought of as a subspace thereof. This situation also
occurs in other places, most notably, in the case of gauge theories.

Now consider a single particle Hilbert space $\mathbb{H}^{(1)}$ with an orthonormal basis
$\{\ket{j}, j=1,\ldots,L\}$, for example where $\ket{j}$ corresponds to the particle being
positioned on site $j$ of a lattice with $L$ sites. We also refer to these single particle
states as modes. To construct a basis for $\mathbb{H}^{(N)}$, we can start from the tensor
product basis of $\widetilde{\mathbb{H}}^{(N)}$ and apply $\hat{S}^{\pm}$ to each of its
$L^N$ elements. Let us henceforth denote these states as

$$ \ket{j_1,j_2,\ldots ,j_N} = \hat{S}^{\pm} \left(\ket{j_1} \otimes \ket{j_2} \otimes \cdots \otimes \ket{j_N}\right)$$

The application of $\hat{S}^{\pm}$ will create certain linear dependences. In particular, states
$ \ket{j_1,j_2, \ldots, j_N}$ that contain the same set of modes $j_k$, i.e. for which the $j_k$'s
are related by a permutation, are equal (up to a sign in the case of $\hat{S}^-$). We can thus select
a single state by ordering the $j_k$ arguments. Furthermore, in the case of $\hat{S}^{-}$, the
state is mapped to zero as soon as two $j_k$ values coincide, so we can eliminate such states. If we
thus restrict the set to states $\ket{j_1,j_2,\ldots ,j_N}$ which are such that the modes are
ordered as $j_1 < j_2 < \ldots < j_N$ (for fermions) or $j_1 \leq j_2 \leq \ldots \leq j_N$ (for bosons),
then we have a linearly indepenent set of states. For fermions, this implies in particular that
we need to have $N \leq L$, there cannot be more fermions in the system then there are linearly
independent modes (single particle states). 

Finally, one can wonder about the normalisation of these states. For fermions, the superposition 
created by $\hat{S}^-$ contains $N!$ terms, which are mutually orthogonal, so that the 
resulting state is normalised (because of the $1/\sqrt{N!}$ prefactor in the definition of 
$\hat{S}^{-}$. More generally, one then finds

$$\braket{i_1 < i_2 < \ldots < i_N | j_1 < j_2 < \ldots < j_N} = \delta_{i_1,j_1} \delta_{i_2,j_2} \cdots \delta_{i_N,j_N}$$

For bosons, the situation is more complicated in the case that some $j_k$ values coincide.
Some of the $N!$ terms created by $\hat{S}^+$ are then equal and contribute differently to 
the norm. If we denote with $n_1, n_2, \ldots, n_L$ the number of $j$ values that equal 
the value $1, 2, \ldots, L$, i.e. the number of particles in mode $1, 2, \ldots, L$,
 then we find

$$\braket{i_1 \leq i_2 \leq \ldots \leq i_N | j_1 \leq j_2 \leq \ldots \leq j_N} = (n_1! n_2! \cdots n_L!) \delta_{i_1,j_1} \delta_{i_2,j_2} \cdots \delta_{i_N,j_N}$$

This more general exprression is also valid for fermions, where every $n_j$ is restricted to be
zero or one. In fact, the values $n_j$ for $j=1,\ldots,L$ completely characterise the state, and
can thus be used to relabel the basis. Instead of specifying the mode $j_k$
that each particle $k=1,\ldots,N$ occuppies (where the labeling of the particles is arbitrary
because they are identical), we can move to a mode-based description and thus specify the 
number of particles in each mode, also known as the mode occupation number. We can then
refer to the basis vectors as

$$\ket{n_1, n_2, \ldots, n_L}$$

where $n_j = 0, 1$ (fermions) or $n_j = 0,1,2, \ldots $ (bosons) and furthermore
$\sum_{j=1}^{L} n_j = N$. Furthermore, we define these states to be normalised to 1, i.e. we
absorb a suitable normalisation factor when defining $\ket{n_1, n_2, \ldots, n_L}$ in terms
of the construction above. 

This way of labelling the basis states now is again reminiscent of a tensor product structure,
i.e.\ we could think of $\ket{n_1, n_2, \ldots, n_L}$ as the tensor product of states $\ket{n_j}$
associated to every mode, and where the Hilbert space associated with such a mode is two-dimensional 
in the case of fermions, or infinite-dimensional in the case of bosons.  However, there is still
a global constraint $\sum_{j=1}^{L} n_j = N$ so that we cannot let the different $n_j$ values
vary completely independently from each other. Furthermore, some caution is now needed as to
what it means to have operators acting on these different "mode Hilbert spaces". The correct
formalism is that of second quantisation, which we introduce next.

```{note}
In many applications, people do still work with the framework of first quantisation, and consider
$N$-particle states constructed by symmetrising or antisymmetrising the tensor product of
$N$ single-particle states, in a so-called independent particle model or approximation.
Such states are quite cumbersome to work with. As can already be seen,
the antisymmetric case is slightly easier and is known as a Slater determinant. Indeed, the 
antisymmetrisation formula is reminiscent of the Leibniz formula of a determinant, and for 
example the inner product between two Slater determinants constructed from $\{\ket{\psi_n},n=1,\ldots,N\}$ and
$\{\ket{\varphi_n},n=1,\ldots,N\}$ is given by the determinant of the matrix containing all
overlaps $\braket{\varphi_m \vert \psi_n}$. Slater determinants form the basis of Hartree-Fock
theory for approximating the state of electrons in an atom or molecule.

The bosonic version occurs in the context of Bose-Einstein condensation and cold atom systems
more generally. In that case, the inner product between two such states gives rise to a
determinant-like formula, but without the minus signs. This construction is known as the
permenant, but unlike the determinant it is very hard to compute in general and really requires
to explicitly sum up all $N!$ terms.
```

## Fock space and second quantisation

However, there is an easier formalism that is furthermore required when dealing with systems
in which the precise number of particles $N$ might fluctuate. While non-relativistic processes
do typically not create new particles, this can still be useful for providing an effective
or approximate description. For example, the BCS theory of superconductivity is constructed
by transforming into a description where the number of Cooper pairs can fluctuate.

We start by defining the Fock space, which is the direct sum of all physical (symmetrised or
antisymmetrised) Hilbert spaces $\mathbb{H}^{(N)}$ for different particle numbers $N$, going
all the way from $N=0$:

$$\mathbb{H} = \bigoplus_{N=0}^{+\infty} \mathbb{H}^{(N)}$$

Note that we have not discussed the case $N=0$ before, as in the previou subsection we started
the construction of $\mathbb{H}^{(N)}$ from a given single particle Hilbert space $\mathbb{H}^{(1)}$.
When there are no particles in the system, there is only a single state in which it can be.
Hence, for $N=0$ particles, the Hilbert space $\mathbb{H}^{(0)}$ is spanned by a single state,
which we typically denote as $\ket{0}$ or $\ket{\Omega}$ and refer to as the *vacuum state*.
We choose this state to be normalised, and thus note that $\ket{0}$ is very different from
a the zero vector of the vector space, which has norm zero.

The Fock space becomes a Hilbert space simply by incorporating the inner product from each
of its summands. States within the different summands of this direct sum are defined
to be orthogonal, i.e. $\braket{\varphi^{(M)} \vert \psi^{(N)}}=0$ for all $M$-particle states
$\ket{\varphi^{(M)}}$ and $N$-particle states $\ket{\psi^{(N)}}$ with $M \neq N$.

Given the above mode occupation numbers, we can now use the states

$$\{\ket{n_1,n_2,\ldots,n_L}, n_j = 0,1\ \text{or}\ n_j = 0,1,2,\ldots\}

as basis, where the global constraint can now be omitted. Aside from losing the global
constraint, we have so far not yet gained anything and the different particle number sectors 
cannot yet easily be related to each other. Hereto, we now introduce a set of operators on
$\mathbb{H}$ which relate these different sectors by creating or annihilating a particle.

TBC

## Quantum-to-classical mapping



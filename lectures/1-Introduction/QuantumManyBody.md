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

## Quantum Mechanics and its Postulates

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
    elements of this vector space with Dirac's ket notation $\ket{\psi}$. In particular, we
    can build linear combinations

    ```{math}
    \ket{\psi} = a \ket{\psi_1} + b \ket{\psi_2}
    ```

    for all $a, b \in \mathbb{C}$ and all $\ket{\psi_1}, \ket{\psi_2} \in \mathbb{H}$.

2.  $\mathbb{H}$ has an inner product, which maps to vecors $\ket{\psi}$ and $\ket{\varphi}$
    onto a scalar $\braket{\varphi|\psi} \in \mathbb{C}$ with the properties of
    * Linearity: $\bra{\varphi} ( a \ket{\psi_1} + b \ket{\psi_2}) = a \braket{ \varphi | \psi_1} + b \braket{ \varphi | \psi_2}$
    * Skew-symmetry: $\braket{ \varphi | \psi} = \braket{ \psi | \varphi}^\ast$
	  * Positivity: $\braket{ \psi | \psi} \geq 0$ with equality only if $\ket{\psi} = 0$.
	
 This last property enables us to define a norm $\lVert \psi \rVert = \lVert \ket{\psi}
 \rVert = \sqrt{\braket{\psi|\psi}}$, which satisfies known properties such as $\lVert \psi
 \rVert = 0 \Leftrightarrow \ket{\psi} = 0$ $\lVert a \ket{\psi} \rVert = \vert a\vert
 \lVert \psi \rVert$ and the triangle inequality $\lVert \ket{\varphi} + \ket{\psi} \rVert
 \leq \lVert \varphi \rVert + \lVert \psi \rVert$.

3.  The final property of metric completeness is a technical requirement that is only
    relevant in infinite-dimensional Hilbert spaces. Firstly, a metric is a notation of
    distance between the elements in $\mathbb{H}$, which is provided by the norm of the
    difference, i.e. $d(\varphi, \psi) = \lVert \varphi - \psi \rVert$.
    
    Completeness of the metric is a specific property that guarantees that certain sequences
    of vectors are guaranteed to have a limit value that also exists in $\mathbb{H}$. This
    is necessary to make sense of e.g. Fourier series.

The state of a quantum system is associated to a ray of vectors, which is the
one-dimensional space $\{ a \ket{\psi} , \forall a \in \mathbb{C}\}$ spanned by a single
(nonzero) vector $\ket{\psi} \in \mathbb{H}$. We will describe the state of the system using
a single representative $\ket{\psi}$ of this ray, which we typically choose such that
$\braket{ \psi | \psi} = 1$. However, this does not fix the representative completely, as we
can still add arbitrary phases $\exp(\mathrm{i}\alpha)$, i.e. $\ket{\psi}$ and
$\mathrm{e}^{\mathrm{i}\alpha} \ket{\psi}$ describe the same state.

The best known Hilbert space from your courses on single-particle quantum mechanics is
probably $L^2(\mathbb{R}^n)$, the Hilbert space for a single quantum particle moving in the
$n$-dimensional coordinate space $\mathbb{R}^n$ (typically $n=1,2,3$). This Hilbert space
corresponds to the space of all square-integrable functions $\psi:\mathbb{R}^d \to
\mathbb{C}: x \mapsto \psi(x)$ and the inner product is given by

$$\braket{\varphi | \psi} = \int_{\mathbb{R}^n} \varphi(x)^\ast \psi(x)\,\mathrm{d} x$$

However, this is already a complicated Hilbert space from a technical perspective. Hilbert
spaces can also be finite-dimensional, i.e. $\mathbb{C}^d$, the space of column vectors of
length $d$, with the standard Euclidean inner product

$$\braket{\varphi | \psi} = \sum_{i=1}^d \varphi_i^\ast \psi_i$$

These Hilbert spaces will be very important in our discussion. The simplest nontrivial case
corresponds to $d=2$ and the associated quantum system is known under various names. It is
often referred to as a qubit in the context of quantum information theory. There are various
ways in which qubits can be physically realised. Another common example of a two-dimensional
Hilbert space is for describing the spin degree of freedom of an electron, or another
particle with spin quantum number 1/2. Such a reduced description (forgetting about the
position) is possible if the electron is localised in space, for example when it is strongly
bound to an atom.

If we do want to describe a particle that moves in space, we might also consider it to exist
only at discrete positions in space, i.e. on a lattice. For example, on a one-dimensional
lattice (a.k.a. a chain) with $L$ sites, the Hilbert space would also correspond to
$\mathbb{H} = \mathbb{C}^L$ and the standard basis vectors $\vert j \rangle, j=1,\dots,L$
correspond to the state of the system if the particle is exactly localised on site $j$. We
can also consider infinitely large lattices, e.g. the one-dimensional chain where there is a
site associated with every $j \in \mathbb{Z}$ (or the n-dimensional hypercubic lattice
$\mathbb{Z}^n$). The resulting Hilbert space is then spanned by the states $\vert j \rangle$
for all $j \in \mathbb{Z}$, and is thus infinite-dimensional but with a straightforward
countably infinite basis.

Of course, our goal is to find the Hilbert space of a many body system. We return to this
question below and devote a complete section to it.

###	**Postulate 2: Observables**

Physical observables of the system correspond to self-adjoint (a.k.a. Hermitian) linear
operators on the Hilbert space $\mathbb{H}$.

An operator $\hat{A}$ on $\mathbb{H}$ is a linear map $\hat{A}:\mathbb{H} \to \mathbb{H}$,
i.e. a map from vectors to vectors that satisfies

```{math}
\hat{A}( a \ket{\varphi} + b \lvert \psi \rangle) = a \hat{A}(\ket{\varphi}) + b
\hat{A}(\ket{\psi})
```

The adjoint of an operator $\hat{A}$ is a new operator $\hat{A}^\dagger$ that is constructed
such that

```{math}
\bra{\varphi} \hat{A} \psi \rangle = \langle \hat{A}^\dagger \varphi \ket{\psi}
```

for all $\ket{\varphi}, \ket{\psi} \in \mathbb{H}$ and where $\vert \hat{A}\psi \rangle =
\hat{A} \ket{\psi}$. This definition requires that $(a_1 \hat{A}_1 + a_2 \hat{A}_2)^\dagger
= a_1^\ast \hat{A}_1^\dagger + a_2^\ast \hat{A}_2^\dagger$ and $(\hat{A}_1
\hat{A}_2)^\dagger = \hat{A}_2^\dagger \hat{A}_1^\dagger$.

A self-adjoint operator is an operator such that $\hat{A}^\dagger = \hat{A}$ or thus

```{math}
\braket{\varphi | \hat{A} \psi } = \braket{ \hat{A}^\dagger \varphi | \psi}
```

for all $\ket{\varphi}, \ket{\psi} \in \mathbb{H}$. Linear combinations of self-adjoint
operators with real coefficients are self-adjoint. The composition of two self-adjoint
linear operators $\hat{A}_1 \hat{A}_2$ is self-adjoint if and only if

```{math}
\left[ \hat{A}_1 , \hat{A}_2 \right] = \hat{A}_1 \hat{A}_2 - \hat{A}_2 \hat{A}_1 = 0,
```

i.e. if the operators also commute. Self-adjoint operators have real eigenvalues, and
eigenvectors associated to distinct eigenvalues are orthogonal. In a finite-dimensional
Hilbert space, self-adjoint operators admit a spectral decomposition

$$\hat{A} = \sum_{n} \lambda_n \hat{P}_n$$

where $\hat{P}_n$ is the spectral projector onto the eigenspace associated with $\lambda_n$.
The spectral projectors satisfy $\hat{P}_n \hat{P}_m = \delta_{n,m} \hat{P}_n$,
$\hat{P}_n^\dagger = \hat{P}_n$ and $\sum_{n} \hat{P}_n = \mathbb{1}$, the identity
operator. If $\lambda_n$ has one-dimensional eigenspace spanned by the eigenvector
$\vert\phi_n\rangle$, then

```{math}
\hat{P}_n = \frac{\vert \phi_n \rangle \langle \phi_n \vert}{\langle \phi_n \vert \phi_n \rangle}
```

where the denominator can be omitted if the eigenvector is normalised.

In the language of matrices, these properties can be rephrased as follows: With respect to
an orthonormal basis choice, self-adjoint operators are represented as hermitian matrices.
Such matrices can be diagonalised by a unitary transformation, or thus, we can construct a
complete basis consisting of eigenvectors. With respect to this basis, the self-adjoint
operator is represented by a diagonal matrix with real values on the diagonal.

### **Postulate 3: Measurements, Expectation Values and Collapse**

Given an observable to which we associate the operator $\hat{A}$, we now need to prescribe
the result of measuring this observable with respect to a system that is in a state
$\ket{\psi}$. The most compact way of describing the result is by stating that, the
expectation value $\braket{\hat{A}}$ (= the mean value of the measurement when averaging
over an ensemble of identical copies of the system) is given by

```{math}
\braket{\hat{A}} = \frac{\braket{\psi \vert \hat{A} \vert \psi}}{\braket{\psi \vert \psi}}
```

By exploiting the fact that this also prescribes the expectation value of all higher moments
$\braket{\hat{A}^k}$, this determines the full probability distribution of the measurement
outcome, and yields the more familiar result: The only possible measurement outcomes are
given by the eigenvalues $\lambda_n$ of $\hat{A}$, and for a system in state $\ket{\psi}$
(now assumed normalized), the probability of obtaining $\lambda_n$ is given by $p_n =
\braket{\psi \vert \hat{P}_n \vert \psi}$ with $\hat{P}_n$ the spectral projector from
above. In the case that $\lambda_n$ has a single (linearly independent) eigenvector
$\ket{\phi_n}$ (also assumed normalised), this amounts to $p_n = \vert
\braket{\phi_n|\psi}\vert^2$.

There is a second part to the measurement postulate, which states that, if the measurement
is immediately repeated (without intermediate dynamics, as described by the next postulate),
then the same measurement outcome is obtained. Because the measurement outcome with respect
to the initial state $\ket{\psi}$ is probabilistic and can yield different results, this
requires that after the first measurement, the state changes is changed. This is the
well-known **collapse** of the wave function. More specifically, if a measurement of
observable $\hat{A}$ is performed in a system with state $\ket{\psi}$ and the measurement
value $\lambda_n$ is obtained, then the state of the system changes to

```{math}
\ket{\psi} \longrightarrow \frac{\hat{P}_n \ket{\psi}}{\lVert \hat{P}_n \ket{\psi}\rVert}.
```

Note that the denominator cannot vanish, as otherwise the probability of having obtained
measurement outcome $\lambda_n$ would have been zero in the first place.

### **Postulate 4: Dynamics**

During time intervals without measurements, the state of an isolated quantum system evolves 
unitarily according to the (first order linear) differential equation

```{math}
\frac{\mathrm{d}\ }{\mathrm{d} t} \ket{\psi(t)} = - \mathrm{i} \hat{H}(t) \ket{\psi(t)}
```

known as the Schr\"{o}dinger equation, where $\hat{H}(t)$ is the Hamiltonian of the system,
which may itself be time-dependent. In the case of a time-independent Hamiltonian, we can
define the evolution operator

```{math}
U(t, t') = \exp\left(-\mathrm{i}(t-t')\hat{H}\right)
```

which relates states at different times via $\ket{\psi(t)} = \hat{U}(t, t') \ket{\psi(t')}$
and is clearly a unitary operator. Clearly, we need to know the Hamiltonian of a system in
order to even start thinking about modelling its quantum properties. We will always assume
that the Hamiltonian is given. In practice, however, the situation can be much more
complicated. Typically, we want to build only an effective quantum description of the system
(e.g. only the electrons, only certain electrons, $\ldots$) and not start all the way down
at the level of fundamental particles and the standard model (which is also only an
effective model valid up to some energy scale).

## The Hilbert Space of Many-Body Physics

All of the previous axioms remain valid for a composite system consisting of several quantum
degrees of freedom. However, we need to know how to describe the state of the system, and 
thus more specifically, how to define the Hilbert space associated to such a system. It
turns out that quantum mechanics forces us to distinguish two cases.

### Distinguisable Particles and Tensor Products

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

### Identical Particles and Pauli's Exclusion Principle

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
$\ket{j_1,j_2,\ldots ,j_N}$ which are such that the modes are ordered as
$j_1<j_2<\ldots<j_N$ (for fermions) or $j_1 \leq j_2 \leq \ldots \leq j_N$ (for bosons),
then we have a linearly independent set of states. For fermions, this implies in particular
that we need to have $N \leq L$, there cannot be more fermions in the system then there are
linearly independent modes (single particle states).

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
each particle $k=1,\ldots,N$ occuppies (where the labeling of the particles is arbitrary
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

## Fock Space and Second Quantisation

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
i.e. $\braket{\Phi| hat{a}_j^+ \Psi} = \braket{hat{a}_j^- \Phi | \Psi}$, we use the simpler
notation $\hat{a}_j$ for the *annihilation operator* and $\hat{a}_j^\dagger$ for the
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

With respect to the normalized basis vectors using number occupation representation, we have

```{math}
\ket{n_1, n_2, \ldots, n_L} = \frac{1}{\sqrt{n_1! n_2! \cdots n_L!}}
(\hat{a}_1^\dagger)^{n_1} (\hat{a}_2^\dagger)^{n_2} \cdots (\hat{a}_L^\dagger)^{n_L}
\ket{\Omega}
```

which can be summarized using

```{math}
\begin{align}
\hat{a}_j \ket{n_1, n_2, \ldots, n_j, \ldots, n_L} &= (\pm 1)^{n_1 + n_2 + \ldots + n_{j-1}} \sqrt{n_j} \ket{n_1, n_2, \ldots, n_j - 1, \ldots, n_L},\\
\hat{a}_j^\dagger \ket{n_1, n_2, \ldots, n_j, \ldots, n_L} &= (\pm 1)^{n_1 + n_2 + \ldots + n_{j-1}} \sqrt{n_j+1} \ket{n_1, n_2, \ldots, n_j + 1, \ldots, n_L}.
\end{align}
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

## Interesting States and Observables in Quantum Many-Body Physics

Having introduced the Hilbert space and Hamiltonian of quantum many-body systems, we still
need to define which states we are actually interested in, and which type of observables we
want to compute for such states. So far, we have only mentioned that isolated systems have a
quantum state which corresponds to a vector (or rather a ray of vectors) in its Hilbert
space $\mathbb{H}$. Before answering this question, we first need to generalize our concept
of a quantum state.

### Quantum States Revisited

More abstractly and generally, the quantum state of a system can be introduced as a map from
observables (operators) to numbers (expectation values). This is typically denoted as $\rho:
\mathrm{End}(\mathbb{H}) \mapsto \mathbb{C}:\hat{A} \to \hat{A}$. Here,
$\mathrm{End}(\mathbb{H})$ is the set of linear operators (a.k.a endomorphisms) on
$\mathbb{H}$. This set is itself a vector space, as we can consider linear combinations of
linear operators. Furthermore, as we can compose two linear operators and obtain a new
linear operator, we have a product operation, which makes $\mathrm{End}(\mathbb{H})$ into an
algebra. Finally, we have defined the concept of the adjoint of an operator, which in
mathematics terminology gives $\mathrm{End}(\mathbb{H})$ the structure of a
$C^\ast$-algebra.

The map $\rho$ that represents a quantum state should have a number of properties, that
generalise those of the case we have encountered so far, where $\rho(\hat{A}) =
\frac{\braket{\Psi\vert \hat{A} \vert \Psi}}{\braket{\Psi | \Psi}}$. In particular, this map
is linear with respect to linear combinations of operators. This implies that it can be
written as $\rho(\hat{A}) = \mathrm{Tr}\left[\hat{\rho}\hat{A}\right]$, where $\hat{\rho}$
is now itself an element of $\mathrm{End}(\mathbb{H})$ (technically, $\rho$ is an element
from the dual space of $\mathrm{End}(\mathbb{H})$). Furthermore, we must have that our state
gives rise to nonnegative and normalised probabilities, which implies that

*   $\rho(\hat{1}) = \mathrm{Tr}\left[\hat{\rho}\right] = 1$
*   $\rho(\hat{P}) \geq 0$ for any projector, and more generally, for any positive definite
    operator $\hat{P}$. This implies that the associated operator $\hat{\rho}$, known as the
    **density operator** or density matrix (when expressed with respect to a chosen basis),
    is itself a positive (and thus self-adjoint) operator, which is furthermore normalised
    to have trace one.

The particular case where the state of the system was given by a vector
$\ket{\Psi}\in\mathbb{H}$ corresponds to
$\hat{\rho}=\frac{\ket{\Psi}\bra{\Psi}}{\braket{\Psi\vert \Psi}}$ and thus satisfies
$\hat{\rho}^2=\hat{\rho}$, i.e. $\hat{\rho}$ is itself a projector. Such states are called
*pure states*. All density operators which do not have this property are called *mixed
states*.

Being positive definite operators, any density operator admits a spectral decomposition of
the form

```{math}
\hat{\rho} = \sum_{n} p_n \ket{\Phi_n}\bra{\Phi_n}
```

where the states $\{\ket{\Phi_n}\}$ form an orthonormal set and the eigenvalues $p_n$
satisfy $\sum_{n} p_n =1 $ and $p_n \geq 0$ (which together also yields $p_n < 1$).

Mixed states arise in the quantum world in two scenarios:

1.	If the system is not isolated, but is rather a subsystem of a larger system and
    interacting with its complement therein. This is discussed in the next section.

2.  Even for an isolated system, it can happen that the state is not exaclty known and one
    must deal with classical uncertaintity and probability. Indeed, a mixed state can be
    interpreted as a statistical ensemble. If the system can be prepared into different (not
    necessarily orthogonal) states $\{\ket{\Psi_1}, \ket{\Psi_2}, \ldots\}$ with
    probabilities $p_1, p_2, \ldots$ that sum up to one, then the state of the system is
    given by
	
	 ```{math}
   \hat{\rho} = p_1 \ket{\Psi_1}\bra{\Psi_1} + p_2 \ket{\Psi_2}\bra{\Psi_2} + \ldots
   ```

Note that this does not necessarily correspond to the spectral decomposition of
$\hat{\rho}$, as the states $\ket{\Psi_i}$ are not necessarily orthogonal. It is nonetheless
a valid density operator. More generally, given two density operator $\hat{\rho}_1$ and
$\hat{\rho}_2$, aany convex combination $\hat{\rho} = p \hat{\rho}_1 + (1-p) \hat{\rho}_2$
with thus $0 \leq p \leq 1$ is a valid density operator.

To a mixed state, we can associate the Von Neumann entropy

```{math}
S(\hat{\rho}) = - \mathrm{Tr}\left[\hat{\rho}\log \hat{\rho}\right] = - \sum_{n} p_n \log(p_n)
```

with $p_n$ the eigenvalues of $\hat{\rho}$. For a pure state, the Von Neumann entropy
evaluates to zero (using $\lim_{x\to 0} x \log x = 0$). The maximal value of the Von Neumann
entropy is obtained when all values $p_n$ are equal so that $\hat{\rho} \sim \hat{1}$.
Because of normalisation, we then have $p_n = 1/d$ with $d$ the Hilbert space dimension and
thus obtain

```{math}
0 \leq S(\hat{\rho}) \leq \log d.
```

In a many-body system, the Hilbert space dimension scales exponentially with the number of
sites or number of degrees of freedom in the system. If we call this quantity the "volume"
of the system, than we can conclude that the maximal value of the Von Neumann entropy is
thus proportional to the volume of the system. 

Depending on the context, the interpretation and meaning of the Von Neumann entropy can
differ, as we discuss below.

### From Tensor Products to Mixed States and Entanglement

Consider a bipartite system composed of two subsystems $A$ and $B$, with thus $\mathbb{H} =
\mathbb{H}^{(A)} \otimes \mathbb{H}^{(B)}$. Now suppose that we are only interested in
measuring observables that act non-trivially on subsystem $A$. This might be the case if
subsystem $A$ is the actual quantum system that we want to model, but it is not isolated and
instead interacting with an environment, corresponding to subsystem $B$. With the axioms so
far, we are forced to include the environment into our discussion. Only the combined system
+ environment can be assigned a pure state $\ket{\Psi}$. However, this seems complete
overkill, as the environment might extend the whole universe and so it will be impossible to
know the complete state $\ket{\Psi}$. Since we are only interested in observables that act
nontrivially on the system $A$, i.e. all observables that we want to measure take the form
$\hat{O}^A = \hat{O} \otimes \hat{1}_B$, and thus we expect that a reduced and simplified
description must exist.

Let us now assume that the Hilbert space of subsystem $A$ is spanned by a basis
$\{\ket{\psi_k}, k=1,\dots, d^A\}$ and the Hilbert space of subsystem $B$ is spanned by a
basis $\{\ket{\varphi_l}, l=1,\ldots, d^B\}$. A reduced description for the system $A$ can
be obtained by observing that we can write

```{math}
\braket{\Psi \vert \hat{O}^A \vert \Psi} &= \mathrm{Tr}\left[\hat{O} \otimes \hat{1}_B \ket{\Psi} \bra{\Psi}\right]\\
&= \sum_{k = 1}^{d^A}\sum_{l = 1}^{d^B} \left(\bra{\psi_k} \otimes \bra{\varphi_l}\right) \left(\hat{O} \otimes \hat{1}_B\right) \ket{\Psi}\bra{\Psi}\left(\ket{\psi_k} \otimes \ket{\varphi_l}\right)\\
&= \sum_{k = 1}^{d^A}] \bra{\psi_k} \hat{O} \left[ \sum_{l=1}^{d^B} \bra{\varphi_l} \ket{\Psi} \bra{\Psi} \ket{\varphi_l}\right] \ket{\psi_k}\\
&= \mathrm{Tr}_A \left[\hat{O} \mathrm{Tr}_B\left(\ket{\Psi}\bra{\Psi} \right)\right]\\
&= \mathrm{Tr}\left[\hat{O} \hat{\rho}^{(A)}\right]
```

Hence, subsystem $A$ an be described in terms of a mixed state that is obtained as

```{math}
\hat{\rho}^A = \mathrm{Tr}_B \ket{\Psi}\bra{\Psi} = \sum_{l=1}^{d^B} \bra{\varphi_l}
\ket{\Psi} \bra{\Psi} \ket{\varphi_l}
```

This construction is known as a *partial trace* and the resulting mixed state of subsystem
$A$ as the *reduced density operator*. It is based on the fact that by using a tensor
product basis for the joint Hilbert space $\mathbb{H} = \mathbb{H}^{(A)} \otimes
\mathbb{H}^{(B)}$, a trace operation leads to a double sum, namely one over all basis
vectors for $\mathbb{H}^A$ and one over all basis vectors for $\mathbb{H}^B$. Hence, the
complete trace can be interpreted as the composition of two partial traces, one over
subsystem $A$ and one over subsystem $B$. If all relevant operators act trivially on $B$,
the partial trace over $B$ can be performed directly on the state $\hat{\rho}^{(AB)} =
\ket{\Psi}\bra{\Psi}$ and gives rise to the reduced density matrix $\hat{\rho}^{(A)}$. Some
notes are in order.
* While we made reference to a specific tensor product basis to define this construction,
  the concepts of reduced density operator and partial trace do not depend on the specific
  choice of basis for $\mathbb{H}^{(A)}$ and $\mathbb{H}^{(B)}$. The construction only
  requires a tensor product basis to expose the tensor product structure of $\mathbb{H}$.
* While we have assumed that the total system is described by a pure state
  $\hat{\rho}^{(AB)} = \ket{\Psi}\bra{\Psi}$. However, for the construction of the reduced
  density opeator as $\hat{\rho}^{(A)} = \mathrm{Tr}_B \hat{\rho}^{(AB)}$ this is not
  necessary.
* With can expand the whole construction with respect to an explicitly chosen basis. If
  $\ket{\Psi} = \sum_{k=1}^{d^A} \sum_{l=1}^{d^B} \Psi_{k,l} \ket{k}\otimes \ket{l}$, we
  find
  
  ```{math}
  \hat{\rho}^{(AB)} = \sum_{k,k'=1}^{d^A} \sum_{l,l'=1}^{d^B} \Psi_{k,l} \Psi_{k',l'}^\ast
  \left(\ket{k}\otimes \ket{l}\right) \left(\bra{k'}\otimes \bra{l'}\right)
  ```
   
  and
  
  ```{math}
  \hat{\rho}^{(A)} = \sum_{k,k'=1}^{d^A} \sum_{l=1}^{d^B} \Psi_{k,l} \Psi_{k',l}^\ast
  \ket{k}\bra{k'}
  ```

If the reduced density operator $\hat{\rho}^A$ is pure, this indicates that the state
$\ket{\Psi}$ was itself a tensor product. In all other cases, the subsystems $A$ and $B$ are
said to be entangled. This entanglement can be quantified by computing the Von Neumann
entropy $S(\hat{\rho}^A)$, which is then called the **entanglement entropy** of the combined
system $A$ and $B$. Indeed, that this entropy is a property of how both subsystems are
entangled follows from the fact that $S(\hat{\rho}^A) = S(\hat{\rho}^B)$, i.e.\ it doesn't
matter whether the Von Neumann entropy of the reduced density operator of subsystem $A$ or
of subsystem $B$ is computed. This is only true if the total system is in a pure state
$\ket{\Psi}$. When also the total system is in a mixed state, because of classical
randomness, then it is harder to differentiate between true quantum entanglement and
classical probability.

To conclude, we analyze the case where the combined system is in a pure state a bit more.
If we again expand $\ket{\Psi}$ with respect to the tensor product basis as 

```{math}
\ket{\Psi} = \sum_{k=1}^{d^A} \sum_{l=1}^{d^B} \Psi_{k,l} \ket{k}\otimes \ket{l}
```

and interpret its expansion coefficients $\Psi_{k,l}$ as the entries of a $d^A \times d^B$
matrix $C$. The reduced density operators can now be written as

```{math}
\hat{\rho}^A = \sum_{k,k'=1}^{d^A} [C C^\dagger]_{k,k'} \ket{k}\bra{k'}\quad
\text{and}\quad\hat{\rho}^A = \sum_{k,k'=1}^{d^A} [C^\dagger C]_{l,l'} \ket{l}\bra{l'}
```

Hence, the reduced density matrices for subsystems $A$ and $B$ are related by the fact that
they correspond to the two different ways in which we can multiply the matrix $C$ with its
Hermitian conjugate $C^\dagger$. It is a well-known result from linear algebra that for two
matrices $A \in \mathbb{C}^{d_1 \times d_2}$ and $B \in \mathbb{C}^{d_2 \times d_1}$, the
square matrices $A B \in \mathbb{C}^{d_1 \times d_1}$ and $BA \in \mathbb{C}^{d_2 \times
d_2}$ have the same set of nonzero eigenvalues, counted with degeneracy. If $d_1 \neq d_2$,
the larger of the two matrices will have additional eigenvalues zero. This result already
proofs the equality $S(\hat{\rho}^{(A)}) = S(\hat{\rho}^{(B)})$.

However, we can even make this more explicit. We can decompose the matrix $C \in
\mathbb{C}^{d^A \times d^B}$ as $C = U S V^\dagger$ with $U$ and $V$ unitary matrices (of
size $d^A \times d^A$ and $d^B \times d^B$ respectively), and $S$ a $d^A \times d^B$ matrix
which only has nonzero entries on the diagonal. Furthermore, the nonzero entries of $S$ can
be chosen positive in descending order. This decomposition is known as the **singular value
decomposition**. For further reference, we denote the diagonal elements of $S$ as $s_i =
S_{i,i}$ for $i=1,\ldots, \mathrm{min}(d^A,d^B)$.

The unitary matrices $U$ and $V$ can be interpreted as basis transforms in the subsystems
$A$ and $B$ respectively, i.e. they define a new basis (which is thus specific to the chosen
state $\ket{\Psi}$), which we denote as $\{\ket{\psi^{(A)}_k}, k= 1,\ldots, d^A\}$ and
$\{\ket{\psi^{(B)}_l}, l= 1,\ldots, d^B\}$. We can then write

```{math}
\ket{\Psi} &= \sum_{i=1}^{\min(d^A,d^B)} s_i \ket{\psi^A_i} \otimes \ket{\psi^B_i}\\
\hat{\rho}^{(A)} &= \sum_{i=1}^{d^A} (s_i)^2 \ket{\psi^A_i} \bra{\psi^A_i}\\
\hat{\rho}^{(B)} &= \sum_{i=1}^{d^b} (s_i)^2 \ket{\psi^B_i} \bra{\psi^B_i}
```

Hence, the reduced matrices appear immidiately in diagonalised form. The singular values
$s_i$, or rather their squares $p_i = (s_i)^2$ are referred to as Schmidt coefficients, and
together make up the **entanglement spectrum**. The entanglement entropy is then given by

```{math}
S = - \sum_{i} p_i \log(p_i)
```

In all of this, it is clear that subsystems $A$ and $B$ were treated on equal footing, and
it does in fact not matter which of the two is chosen to probe the entanglement structure of
the state.

If the entanglement entropy evaluates to zero, the two subsystems $A$ and $B$ are said to be
unentangled, and the state $\ket{\Psi}$ actually factorises as a tensor product
$\ket{\psi^A} \otimes \ket{\psi^B}$. As soon as the entanglement entropy is nonzero, the two
subsystems are entangled. As stated above, the entropy is upper bounded, in this case by the
logarithm of the smallest of the two Hilbert space dimensions $d^A$ or $d^B$. Thus, if $d^A
\leq d^B$, we find that the entanglement entropy satisfies

```{math}
0 \leq S \leq \log(d^A).
```

It turns out that states that are randomly selected from the Hilbert space typically have
an entanglement entropy that is close to maximal.

### Quantum Many-Body Physics at Finite Temperature

With the concept of mixed states and quantum entanglement at hand, we can now discuss
physically interesting states. Firstly, when considering a system that is in contact with a
large environment that acts as a heat bath at temperature $T$, it will equilibrate. The
state of the system at equilibrium is then given by the so-called Gibbs state

```{math}
\hat{\rho} = \frac{1}{Z(\beta)} \mathrm{e}^{-\beta \hat{H}}
```

where the normalization factor

```{math}
Z(\beta) = \mathrm{Tr}\left[\mathrm{e}^{-\beta \hat{H}}\right]
```

is typically referred to as the **partition function**. Here, $\beta = \frac{1}{k_B T}$ with
$T$ de temperature and $k_B$ Boltzmann's constant. Henceforth, we simply refer to $\beta$ as
*inverse temperature*.

The Von Neumann entropy of the Gibbs state corresponds to the thermodynamical notion of
entropy. In a many-body system, the thermal entropy at finite temperature will be extensive
and thus scale with the volume of the system, just like the energy expectation value, so
that together the free energy $E - T S$ is minimised.

At infinite temperature ($\beta=0$), we obtain $\hat{\rho} \sim \hat{1}$ and the Von Neumann
entropy reaches its upper bound. At zero temperature ($\beta \to +\infty$), we obtain
$\hat{\rho} \sim \hat{P}_0$, with $\hat{P}_0$ the projector onto the eigenspace of lowest
energy. Hence, in that case the Von Neumann entropy is given by $\log(d_0)$, with $d_0$ the
lowest energy eigenvalue, i.e. the dimension of the ground state subspace. Most quantum
lattice systems have a single or at least a small number of linearly indepenent ground
states, so that $d_0$ is a small number independent of the system size. However, there are
also cases where the number of ground states scales exponentially with the system size, and
the thermal entropy remains extensive at zero temperature. This then constitutes a violation
of the infamous third law of thermodynamics.

While the heat bath or environment with which the system interacts is in practice typically
much "larger" (in terms of number of degrees of freedom and thus Hilbert space dimension),
we can use a property that any mixed state can be obtained as the reduced density operator
from a pure state in a Hilbert space that is the tensor product of two copies of the
system's Hilbert space, or thus, where the environment is just exactly as large as the
system. Writing a mixed state $\hat{\rho}$ as the reduced density operator of a pure state
$\ket{\Psi}$ in a Hilbert space $\mathbb{H} = \mathbb{H}^S \otimes \mathbb{H}^E$ is known as
a *purification*. With respect to the purification, all expectation values can be obtained
as

```{math}
\mathrm{Tr}\left[\hat{O} \hat{\rho}\right] = \braket{\Psi \vert \hat{O} \otimes \hat{1}_E \vert \Psi}
```

which can be an advantage if one has an efficient (mathematical, computational, …) formalism
for working with pure states. We can thus always construct such a purification by just
taking the environment to be a copy of the system. Note that the environment in this
construction is merely an auxiliary tool, and has no physical meaning or relation to the
actual environment. If the Hilbert space of the system is described by a basis $\{\ket{j},
j=1,\dots, d\}$, then we can first build a purification if the infinite temperature state as

```{math}
\ket{\Psi_0} = \sum_{j} \ket{j}_S \otimes \ket{j}_E
```

In particular, if the system is a many-body system and $\ket{j}$ is itself already a tensor
product basis state $\ket{j_1} \otimes \ket{j_2} \otimes \ldots$, we can organise the
environment so that matching tensor product factors between system and environment are taken
together. This has the advantage that the infinite temperature state can be written as

```{math}
\ket{\Psi_0} = \left(\sum_{j_1} \ket{j_1}_S \otimes \ket{j_1}_E\right) \otimes \left(\sum_{j_2} \ket{j_2}_S \otimes \ket{j_2}_E\right) \otimes \ldots
```

and still has an overall tensor product structure. A purification of the finite temperature
state can then be obtained as

```{math}
\ket{\Psi_\beta} &= \frac{1}{\sqrt{Z(\beta)}} \sum_{j=1}^{d} \left(\mathrm{e}^{-\frac{\beta}{2} \hat{H}} \ket{j}_S \right)\otimes \ket{j}_E \\
&= \frac{1}{\sqrt{Z(\beta)}} \exp\left(-\frac{\beta}{2} \hat{H} \otimes \hat{1}_E \right) \ket{\Psi_0}
```

This expression now looks remarkibly similar to how to compute a time-evolved state, by
replacing $-\mathrm{i} t \mapsto -\beta/2$. Hence, methods that solve Schrödinger's equation
and are sufficiently general to also work with imaginary values of the time coordinate can
be used to prepare thermal states.

```{note}
Purifications of thermal states are also referred to as *thermofield double states*,
especially in the high energy physics literature, i.e. in the context of quantum field
theory, holography and quantum gravity.
```

### Quantum Many-Body Physics at Zero Temperature

As quantum effects are most pronounced at zero temperature, we typically assume to be
operating in this regime. It then follows that we are mostly interested in the lowest energy
states of the Hamiltonian, and in particular in the ground state(s).

A trivial but nonetheless important property of ground states is that they can easily be
characterised as states that minimise the expectation value $\braket{\Psi \vert \hat{H} |
\Psi}$. Indeed, this forms the basis for the variational principle. If we have a set of
trial states, in which there are a number of free parameters, than we can construct an
approximation to the ground state by 'simply' minimising the energy expectation value of the
Hamiltonian with respect to these free parameters. How good this approximation is in
practice depends on the properties of both the Hamiltonian and the trial states. However,
one way to quantify the quality of the ground state approximation is by computing the energy
variance

```{math}
\braket{\Psi \vert (\hat{H} - \braket{\Psi \vert \hat{H} \vert \Psi})^2 \vert \Psi} =
\braket{\Psi \vert \hat{H}^2 \vert \Psi} - \braket{\Psi \vert \hat{H} \vert \Psi}^2.
```

The Hamiltonians for quantum lattice systems that we are interested in, will typically
contain a sum of terms where every individual term acts nontrivially only in a small patch
of the lattice. For such Hamiltonians, the low-energy states have a special property.
Hereto, we consider arbitrary bipartitions of the system, where one subsystem corresponds to
a (connected) region of sites, whereas the complement, i.e.\ the remaining sites in the
lattice, make up the second subsystem. The ground state will in general not factorize into a
tensor product, as it contains correlations and entanglement between these two (arbitrarily
chosen) subsystems. As pointed out above, the maximal value of entanglement entropy is given
by the logarithm of the smallest of the two Hilbert space dimensions. Assuming that our
chosen region of sites is smaller than its complement, its Hilbert space will itself scale
exponentially with the number of sites in that region, i.e.\ with its volume. Hence, the
entanglement entropy computed for such a bipartition has an upper bound that is proportional
to the volume of the subsystem. As was also stated above, random states typically satisfy
this upper bound. However, the special property of low-energy states of locally interacting
Hamiltonians is exactly that they have much less entanglement. They typically have an
entanglement entropy that only scales with the common area between the subsystem and its
complement. This scaling behavior is referred to as the **area law of entanglement entropy**
and provides the key motivation for approximating such low energy states using a tensor
network decomposition. It indicates that the most important quantum correlations are short
range (just like the interactions that generate them), and are thus situated across the
boundary connecting the subsystem and its complement. However, this does not exclude that
there is also a small amount of nontrivial long-range correlations in the system.

Aside from the ground state, one might also be interested in the first excited states, as
these will be important for understanding how the system at zero temperature reacts to
external perturbations. In a macroscopically large many-body system, one should not expect
that the energy spectrum consists of a number of discrete levels with gaps in between. The
lowest-energy excited states in a quantum lattice system can typically be given a
particle-like interpretation, as is well known from quantum field theory. They correspond to
small bumps of energy, i.e. they can be thought of as perturbations of the ground state in a
small region and thus have an anergy cost that does not scale with the system size. However,
because of kinetic energy-like terms, actual eigenstates will not correspond to having this
energy bump in a localized region, but will rather be in a superposition where the
"particle" is delocalized. In particular, in the case of a translation invariant system, the
eigenstates will also be momentum eigenstates, and thus describe a particle that is in a
momentum superposition across the lattice. The energy (surplus) of such a particle like
excitation will thus be a number $\epsilon(k)$ that is of order $1$ independent of the
system size, and that depends on the specific momentum. If $\epsilon(k)$ is everywhere lower
bounded by some value $\Delta$ (again indepenent of system size), then the system is said to
be gapped. However, in some systems, $\epsilon(k)$ can become zero for particular values of
$k$. Such systems are called gapless, and they often correspond to phase transition points,
where the nature of the ground state radically changes if the parameters in the Hamiltonian
are varied.

The energy spectrum of a quantum lattice system will then consist of one or a few ground
states, the energy of which is an extensive number that is most easily expressed as some
energy density per site. To study the excited states, it is then convenient to shift the
energy scale such that the ground state energy is zero. The lowest excited states will then
correspond to particles which can be created at a certain momentum. In an energy-momentum
diagram, their dispersion relation $\epsilon(k)$ will apear as an isolated band. Note that a
system can have different types of such particle-like excitations, each with their own
dispersion relation. Higher up in the energy spectrum we start to obtain regions
corresponding to states with two- or more particles, that are travelling independent from
each other. For such states, the energy can be obtained simply as the sum of the indivual
particles in the state, and since the relative momentum of the particles can change while
keeping the total momentum fixed, the energy in such states can also vary continuously (at
least in the thermodynamic limit).

### Quantum Dynamics and Quenches

Aside from low-energy eigenstates of the Hamiltonian, we are often also interested in states
that have a non-trivial time depence. One particular use case is where one starts from the
ground state $\ket{\Psi_0}$ of a certain Hamiltonian $\hat{H}_0$, and then some parameters
in the Hamiltonian are suddenly changed, so that the Hamiltonian now corresponds to a new
operator $\hat{H}_1$. This sudden change is reminiscent of quenching a system, and such a
setup is called a global quench. We then want to compute

```{math}
\ket{\Psi(t)} = \exp(-\mathrm{i} t \hat{H}_1) \ket{\Psi_0}.
```

With respect to $\hat{H}_1$, the state $\ket{\Psi_0}$ will no longer be a ground state and
most likely not even be an eigenstate. However, it will have a certain (extensive)
energy expectation value that is preserved throughout the evolution.

In terms of entanglement and correlations, even when $\ket{\Psi_0}$ is a state with an area
law entanglement scaling (because it is a low-energy state of another local Hamiltonian,
namely $\hat{H}_0$), the entanglement in the state will grow rapidly with time. Indeed, the
bipartite entanglement entropy for will tend to grow linearly with time. From the
perspective of a given subsystem, its entropy will grow from an initial value proportional
to the area of the subsystem, until it saturates at a value that is propertional to the
volume of the subsystem. This process is known as thermalisation, as it turns out that at
that point, the subsystem is locally indistinguishable from a Gibbs state with a temperature
set by the energy density of the initial state. The subsystem has thermalized with respect
to its complement behaving as an environment or heat bath. The larger the subsystem, the
longer it will take before thermalisation is complete, and the overall state of the global
system remains a pure state, albeit a highly entangled one.

### Observables and Static and Dynamic Correlation Functions

TO BE WRITTEN!

## Quantum-to-Classical Mapping

In this final section, we introduce a general technique that essentially enables us to map
any quantum lattice system in $d$ dimensions to a classical partition function in $d+1$
dimensions, up to some caveats that we will return to at the end of this section.

TO BE WRITTEN!

<!-- 
### Suzuki-Trotter decomposition

Remember that thermal expectation values are given by

$$\braket{\operator{O}} = \tr\left[\operator{O} \mathrm{e}^{-\beta \hat{H}}\right]/Z(\beta)= \tr\left[\mathrm{e}^{-\beta \hat{H}/2} \operator{O} \mathrm{e}^{-\beta \hat{H}/2}\right]/Z(\beta)$$

with the thermal partition function $Z(\beta)$ given by

$$Z(\beta) = \tr \mathrm{e}^{-\beta \hat{H}} = \sum_{s= \pm 1} \braket{s|\mathrm{e}^{-\beta \hat{H}}|s}.$$

The ground state physics is encoded in the limit $\beta \to \infty$. Note that, if the
system has a unique ground state, we can obtain the ground state $\ket{\Psi}$ of a quantum
system by starting from essentially a random state $\ket{\Phi}$ and evolving it in imaginary
time $\tau = -\ic t$ for sufficiently long

$$\ket{\psi_0} \sim \lim_{\tau \to \infty} \mathrm{e}^{-\tau \hat{H}} \ket{\phi}$$

Expanding the initial state $\ket{\phi}$ in the energy eigenbasis of $\hat{H}$, we see that the only condition is that it is not orthogonal to the ground state (subspace). In addition, the ground state will be well approximated if $\tau \Delta E \gg 1$, with $\Delta E=E_1 - E_0$ the energy gap. This imaginary time evolution also forms the basic ingredient of several numerical algorithms for approximating ground states of quantum many body systems, often in combination with the Suzuki-Trotter decomposition which is introduced below.

Using this approach, the following expression for the ground state expectation value of on
operator $\operator{O}$ is obtained

$$\braket{\hat{O}} = \lim_{\tau\to\infty} \braket{\phi\vert \mathrm{e}^{-\tau \hat{H}} \operator{O} \mathrm{e}^{-\tau \hat{H}}\vert \phi}/\braket{\phi\vert\mathrm{e}^{-2\tau \hat{H}} \vert\phi}$$

This expression can be compared to the thermal expectation value with $\beta = 2\tau$; the
only difference is in the boundary conditions.

 For a quantum many body system, taking the exponential is as hard as determining the full
diagonalisation of the hamiltonian, which is impossible due to the exponentially large
Hilbert space. If the hamiltonian is a sum of local terms, each of these terms can be
exponentiated easily, but for arbitrary $\tau$ there is no relation ship between $\exp(-\tau
\sum_{i} \operator{h}_i)$ and the individual $\exp(-\tau \operator{h}_i)$, unless the
different $\operator{h}_i$ commute. However, for an infinitesimal time step $\epsilon$, we
can use to Baker-Campbell-Hausdorff formula (or better yet, the Zassenhaus formula) to
obtain $\exp(-\epsilon \sum_{i} \operator{h}_i) = \prod_i \exp(-\epsilon \operator{h}_i) +
\mathcal{O}(\epsilon^2)$. This then leads to the \emph{Suzuki-Trotter decomposition}
\begin{equation} \exp\left(-\tau \sum_{i} \operator{h}_i\right) = \lim_{M\to\infty} \left(
\mathrm{e}^{-\frac{\tau}{M} \sum_i \operator{h}_i} \right)^M =\lim_{M\to\infty}
\left(\prod_i \mathrm{e}^{-\frac{\tau}{M} \operator{h}_i} + \order(M^{-1})\right)^M
\end{equation} Note that splitting the time interval $[0,\tau]$ into small segments
$\epsilon = \tau/M$ is also the starting point for deriving a path integral representation
of the quantum partition function. The next step is to insert a resolution of the identity
in between the $N$ different factors, where the labels of the basis will behave as classical
degrees of freedom. For obtaining a path integral, the basis should be labeled by a number
of continuous degrees of freedom, which can then become continuous functions of time in the
limit $\epsilon\to 0$. Here, instead, we will keep $\epsilon$ small but finite, and use a
discrete basis. \subsection{From quantum to statistical mechanics} Let's start with a
quantum system in $d=0$, i.e.\ a small number of spins, or in particular, a single spin,
described by a hamiltonian \begin{equation} \hat{H} = - h_x \sigma^x - h_z \sigma^z
\end{equation} While we could in principle exponentiate $\hat{H}$ directly, we will treat it
using the Suzuki-Trotter decomposition. Throughout the remainder of this section, we will
use the $\sz$ basis, which we denote as $\ket{1} = \ket{\uparrow}$ and $\ket{-1} =
\ket{\downarrow}$. Inserting resolutions of the identity, we write \begin{equation} Z(\beta)
= \tr \mathrm{e}^{-\beta \hat{H}} = \sum_{\{s_k\}=\pm 1} \braket{s_1|\mathrm{e}^{-\epsilon
\hat{H}}|s_2}\cdots \braket{s_M-1| \mathrm{e}^{-\epsilon \hat{H}}|s_{M}} \cdots \braket{s_M|
\mathrm{e}^{-\epsilon \hat{H}}|s_{M+1}} \end{equation} with $s_{M+1} = s_1$, $M\epsilon
=\beta$, and where \begin{align*} \braket{s_{i}|\mathrm{e}^{-\epsilon H}|s_{i+1}} &=
\braket{s_{i}|\mathrm{e}^{-\epsilon H}|s_{i+1}} \approx \braket{s_{i}|\mathrm{e}^{\epsilon
h_z \sz}\mathrm{e}^{\epsilon h_x \sx}|s_{i+1}}\\ &= \mathrm{e}^{\epsilon h_z s_i}
\braket{s_{i}|\cosh(\epsilon h_x) \one + \sinh(\epsilon h_x) \sx |s_{i+1}}\\ &=
\mathrm{e}^{K s_{i}s_{i+1} + h s_i + f_0}. \end{align*} \begin{exercise} Prove that the
equality on the last line is obtained with , $K = -\frac{1}{2}\log \tanh(\epsilon h_x)$, $h
= \epsilon h_z$ and $f_0 = \frac{1}{2}\log[\cosh(\epsilon h_x)\sinh(\epsilon h_x)]$.
\end{exercise}\\ We thus obtain $Z(\beta) = \sum_{s_k} \mathrm{e}^{\sum_{i=1}^{M}K s_i
s_{i+1} + h s_i}$, the partition function of the one-dimensional classical Ising model with
periodic boundary conditions. Indeed, $\braket{s_{i}|\mathrm{e}^{-\epsilon H}|s_{i+1}}$ does
exactly correspond to the transfer matrix, and diagonalising the transfer matrix is the most
straightforward approach to solving the one-dimensional classical Ising model.

\subsection{Higher dimensional generalisation}

We now apply the same approach to the transverse field Ising model in $d$ dimensions, on a hypercubic lattice. We separate the hamiltonian in two parts according to
\begin{equation}
	\hat{H}_{\text{TFIM}} = \left(-J\sum_{\braket{i,j}} \sz_i \sz_j  - h_z \sum_{i} \sz_i\right) + \left(- h_x \sum_{i} \sx_i\right)= \operator{H}_1 + \operator{H}_2
\end{equation}
Note that $\operator{H}_1$ and $\operator{H}_2$ in itself contain commuting terms, but of course don't mutually commute. We follow the same strategy, and will in every (imaginary) time step introduce a resolution of the identity using the tensor product $\sz$ basis. We now denote the basis at time step $k$ as $\ket{\{s_{i,k}\}}$, where $i$ labels a site in the $d$ dimensional lattice hosting the quantum degrees of freedom, and $k$ labels points along the imaginary time axis, which emerges as a new dimension in the problem. We find
\begin{align*}
\exp(-\epsilon \operator{H}_1)\ket{\{s_{i,k}\}} = \exp(\epsilon J \sum_{\braket{i,j}} s_{i,k} s_{j,k}+\epsilon h \sum_i s_{i,k}) \ket{\{s_{i,k}\}}
\end{align*}
as $\operator{H}_1$ is diagonal in this basis, and 
\begin{align*}
\braket{\{s_{i,k}\}|\exp(-\epsilon \operator{H}_2)|\{s_{i,k+1}\}} = \prod_i \braket{s_{i,k}|\mathrm{e}^{-\epsilon h_x \sx_i} | s_{i,k+1}}
 \sim \exp(K_\perp\sum_{i} s_{i,k} s_{i,k+1})
\end{align*}
with $K_\perp = \log \tanh(\epsilon h_x)$ as before. Here, we have now ignored an overall proportionality factor, which is irrelevant when using the partition function to compute expectation values. With this, we find
\begin{equation}
	Z(\beta) = \sum_{\{s_i,k\}} \exp\left(\sum_{k=1}^{M}\sum_{i} K_\perp s_{i,k} s_{i,k+1}+\sum_{k=1}^M\sum_{\braket{i,j}} K_\parallel s_{i,k} s_{j,k} + \sum_{k=1}^{M}\sum_{i} h s_{i,k} \right)
\end{equation}
with $K_{\parallel} = \epsilon J$ and $h = \epsilon h_z$. We thus find the partition function of the classical Ising model in $d+1$ dimensions with anisotropic interaction strengths, periodic boundary condition in the imaginary time direction and a number of sites in the time direction given by $M = \beta/\epsilon$. Hence, the ground state regime $\beta\to \infty$ corresponds to the thermodynamic limit in this additional time direction of the corresponding classical system, and we will see that there are many similarities (or actually, equivalences) between between quantum phenomena in $d$ dimensions and classical phenomena in $D=d+1$ dimensions. On the other hand, when the quantum system is at finite temperature $\beta$, the additional dimension is finite and never in the thermodynamic limit. In that case, this extra dimension cannot cause new non-analyticities in the partition function and finite temperature quantum systems in $d$ dimensions are very similar to classical systems in $d$ dimensions.

In particular, the ground state of the quantum Ising chain has an extended symmetry broken
regime $g\in[0,1)$, matching the ferromagnetic phase at low temperatures in the
two-dimensional classical Ising model. One can show that also the order parameter between
the quantum chain and the classical model is equivalent. At finite temperature however, the
quantum Ising chain does not exhibit symmetry breaking, exactly like the one-dimensional
classical Ising model. It is not surprising that the $d=1$ quantum phase transition itself
is also completely equivalent to the $D=2$ finite temperature phase transition, which also
fits within the broader context of universality, which we briefly discuss in the next
section. The duality mapping in the quantum Ising chain is furthermore equivalent to the
well-known Kramers-Wannier duality in the classical Ising model.

This quantum to classical mapping can also be inverted.\footnote{This does not imply that
the mapping is unique. Choosing a different time step $\epsilon$ or a different basis, e.g.\
the $\sx$ instead of $\sz$ basis will result in a different classical model.} Taking a
codimension $1$ slice out of a $D$-dimensional classical partition function, one obtains a
transfer matrix which can be interpreted as the exponential of a quantum hamiltonian acting
on the Hilbert space of a $d=D-1$ dimensional quantum system. In fact, the modern
perspective on the exact solution of the two-dimensional classical Ising model by Onsager is
exactly by this transfer matrix approach, which is subsequently diagonalised using the
Jordan Wigner and free fermion approach from Section~\ref{s:isingexact}.

It is clear that the quantum to classical mapping is not specific to the quantum Ising model
and can be applied to any hamiltonian. The path integral representation of the partition
function fits within the same scheme, and only differs in the fact that the limit
$\epsilon\to 0$ is taken such that the additional dimension becomes continuous. This is
particularly natural if also the spatial dimensions of the quantum dimension are continuous,
i.e.\ if we have a quantum field theory. In this case, a $D=d+1$ dimensional classical
statistical field theory is obtained. For relativistic quantum field theories, where it is
common practice to explicitly count the time dimension together with the space dimensions,
imaginary time evolution leads to an action, which is equivalently a hamiltonian of a
classical field theory in $D=d+1$ spatial dimensions, with full Euclidean invariance.

One might thus wonder if there is anything new to be learned from studying quantum ground
states. First of all, there is one important catch which we have overlooked so far. There is
no guarantee that the above process yields a classical partition function with Boltzmann
weights which are positive, or even real. While this may seem like a technical detail, it is
of major importance. The quantum to classical mapping is the basis behind the Quantum Monte
Carlo method, one the most successful numerical methods for studying quantum many body
systems. One maps the quantum problem to a classical partition function and then uses one of
the many flavours of Monte Carlo sampling. However, with non-positive Boltzmann weights, the
interpretation of a probability distribution is lost and no efficient sampling procedure can
be designed, as samples might annihilate each other. This is known as the \emph{sign}
problem.

Secondly, for many non-relativistic quantum systems, the anisotropy between the imaginary
time direction and the spatial dimensions in the corresponding classical system cannot be
ignored, even at the critical point. In those cases, critical correlations behave
differently in the spatial and the time direction, which is characterised by a dynamical
critical exponent $z$. The case $z=1$ corresponds to the case where the critical point has
(emergent) rotation/Lorentz invariance between time and space.

A final reason to study quantum systems directly is that certain concepts are more natural
in that setting. In particular, the last 15 years, ideas from quantum information theory,
and in particular the concept of entanglement, have made their way into the standard toolbox
to study and characterize quantum many body systems.
 -->

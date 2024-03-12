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

(observables)=
# Interesting States and Observables in Quantum Many-Body Physics

Having introduced the Hilbert space and Hamiltonian of quantum many-body systems, we still
need to define which states we are actually interested in, and which type of observables we
want to compute for such states. So far, we have only mentioned that isolated systems have a
quantum state which corresponds to a vector (or rather a ray of vectors) in its Hilbert
space $\mathbb{H}$. Before answering this question, we first need to generalize our concept
of a quantum state.

## Quantum States Revisited

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

2.  Even for an isolated system, it can happen that the state is not exactly known and one
    must deal with classical uncertaintity and probability. Indeed, a mixed state can be
    interpreted as a statistical ensemble. If the system can be prepared into different (not
    necessarily orthogonal) states $\{\ket{\Psi_1}, \ket{\Psi_2}, \ldots\}$ with
    probabilities $p_1, p_2, \ldots$ that sum up to one, then the state of the system is
    given by
	
    ```{math}
    \hat{\rho} = p_1 \ket{\Psi_1}\bra{\Psi_1} + p_2 \ket{\Psi_2}\bra{\Psi_2} +
    \ldots
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

(entanglement)=
## From Tensor Products to Mixed States and Entanglement

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

This particular way of writing the bipartite state $\ket{\Psi}$ is known as the Schmidt decomposition.
As a result, the reduced matrices appear immidiately in diagonalised form. The singular values
$s_i$, or rather their squares $p_i = (s_i)^2$ are referred to as **Schmidt coefficients**, and
together make up the **entanglement spectrum**. The entanglement entropy is then given by

```{math}
:label: entanglement_entropy
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

## Quantum Many-Body Physics at Finite Temperature

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

(zero_temp)=
## Quantum Many-Body Physics at Zero Temperature

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

## Quantum Dynamics and Quenches

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

## Observables and Static and Dynamic Correlation Functions

There are a number of typical observables that we want to measure for a given quantum state
(ground state or thermal state) of a quantum lattice system. The first are operators which
have the same structure as the Hamiltonian, in being given by a sum of terms where every
individual term acts nontrivially only on a single site, or a small patch of neighbouring
sites. Furthermore, these terms all act identically, except that they are translated to the
different patches that make up the lattice. Typical examples include the energy itself or
specific contributions to it (kinetic energy, interaction energy, … ) as well as the total
number of particles, a total magnetisation, and other similar quantities.

Let us use the transverse field Ising model as an example. An interesting quantity in the
transverse field Ising model is the longitudinal magnetisation, given by 
$\hat{S}^z=\frac{1}{2} \sum_{n} \sigma^z_n$. The expectation value of such operators is
extensive, and so we are typically interested in the associated density, which is obtained
by dividing out the volume factor. With respect to a translation-invariant state, this is
equivalent to simply measuring the expectation value of a single term, i.e.
$\frac{1}{2} \sigma^z$ in the case of the longitudinal magnetisation, which is thus a local
operator. The precise position where it acts does then not really matter.

Such quantities come up for example to measure potential symmetry breaking, where an operator that should
have zero expectation value given the symmetries of the Hamiltonian, actually acquires a
nonzero expectation value. Indeed, given that the Ising Hamiltonian has
the property that $\left[\hat{H},\hat{U}\right] = 0$, where $\hat{U} = \bigotimes_{n} \sigma^x_n$,
we also expect the ground state $\ket{\Psi_0}$ to satisfy 
$\hat{U} \ket{\Psi_0} \sim \ket{\Psi_0}$, where the proportionality factor can only be plus
or minus one, due to $\hat{U}^2 = \hat{1}$. The magnetisation in the $z$-direction, on the
other hand, satisfies $\hat{U}^\dagger \hat{S}^z \hat{U} = - \hat{S}^z$, so that we expect

$$\braket{\Psi_0 \vert \hat{S}^z \vert \Psi_0} = 0$$

Indeed, if the ground state is unique, this is precisely what happens. However, it can
happen that there are multiple linearly independent ground states. In that case, the
restriction of $\hat{S}^z$ into the ground subspace can be nontrivial, and there exist
specific ground state choices for which the expectation value is nonzero. For reasons that
go beyond what can be explained here, it are typically these states which are easiest to
create or approximate (they have lower entanglement). Such operators that can characterise
the presence of symmetry breaking are referred to as *order parameters*.

Nonetheless, the use of local operators to probe the ground state properties is somewhat
limited. A different class of observables that are of typical interest are (static)
correlation functions, which take the form

$$ C^{A,B}_{i,j} = \braket{\Psi\vert \hat{A}_i^\dagger \hat{B}_j\vert \Psi} $$

where $\hat{A}_i$ is a local operator acting on or in the neighbourhood of site (or unit
cell) $i$ and $\hat{B}_j$ is a local operator acting on or in the neighbourhood of site $j$.
Typically, the operators $\hat{A}$ and $\hat{B}$ are chosen such that their local expectation
value is zero. If this is not the case, we can subtract these local expectation values by
redefining

$$ C^{A,B}_{i,j} = \braket{\Psi\vert (\hat{A}_i - \braket{\Psi\vert\hat{A}_i\vert \Psi})^\dagger (\hat{B}_j - \braket{\Psi\vert\hat{B}_j\vert \Psi})\vert \Psi} = \braket{\Psi\vert \hat{A}_i \hat{B}_j\vert \Psi} -  \braket{\Psi\vert\hat{A}_i\vert \Psi} \braket{\Psi\vert\hat{B}_j\vert \Psi}.$$

We are then interested in the dependence of this correlation function on the positions $i$
and $j$. In particular, in a translation invariant system, it is only the relative lattice
vector from site $i$ to site $j$ on which this quantity depends. It can be proven that if
$\ket{\Psi}$ is the unique ground state of a gapped local Hamiltonian, then the asymptotic
behaviour of such correlation functions is that they decay exponentially in the distance
between the two sites. This exponential thus defines a length scale $\xi$ via $\exp(-d/\xi)$
with $d$ the relevant distance. The length scale $\xi$ is known as the correlation length of
the system.

When the system is gapless, static correlation functions still go to zero in the limit of
infinite separation distance, but rather decay as an algebraic function of the distance,
i.e. they give rise to power laws. The exponents that appear in these power laws do
typically have universal values that are set by general properties such as the number of
spatial dimensions, the global symmetries in the system, etc. In particular for the case of
one-dimensional system, there is a rich literature and well developed framework for
analysing such gapless systems using methods from conformal field theory.

Finally, in systems with potential symmetry breaking, the static correlation function of the
order parameter with itself is an extremely useful diagnostic. In particular, when the system
has symmetry breaking, the large distance limit of the correlation function does not vanish
and the system is said to contain *long range order*. Unlike the expectation value of the local
order parameter, which can have a nonzero expectation value for particularly chosen
symmetry breaking ground states but is still zero for other choices of ground states, the
value of the correlation function and its large distance limit is insensitive to the specifically
chosen ground state out of the ground subspace that it is computed with. For the transverse-field
Ising model, symmetry breaking will thus be present whenever

$$ C_{i,j} = \braket{\Psi \vert \sigma^z_i \sigma^z_j \vert \Psi} $$

does not decay to zero limit for large distance between sites $i$ and $j$. The limiting
value of this correlation function can then be considered as $m^2$, i.e. the square of the
local magnetisation that would be measured in some states of the ground subspace.

Strictly speaking, the ground state static correlation function does not provide information
about excited states or other dynamical information of the Hamiltonian. In most physical
system, it however does provide some qualitative information. Since the static correlation
function, considered as a matrix with rows $i$ and columns $j$, has a particular structure
resulting from translation invariance, it can be diagonalised by a (multidimensional)
discrete Fourier transform. The resulting eigenvalues depend on the lattice momentum
$\kappa$ and are known as the static structure factor $S(\kappa)$. In particular, in the
case of a gapped system with unique ground state, these values are well defined for all
$\kappa$. Nonetheless, it can be argued (using different techniques) that maxima for
$S(\kappa)$ will correspond to momenta where the single particle excitations have minima in
their dispersion relations. For critical systems, $S(\kappa)$ can also have algebraic
divergences, whereas in the case of long range order, $S(\kappa)$ will contain a Dirac-delta
type of divergence, typically at zero momentum, unless there is some spatially repeating
pattern in the way symmetry is broken (and thus also translation invariance is broken).

More detailed quantitative information about the spectrum of excited states is contained in
the time-dependent correlation function 

```{math}
G^{A,B}_{i,j}(t) &= \braket{\Psi_0\vert  \hat{A}_i(t)^\dagger \hat{B}_j(0)\vert \Psi_0}\\
&=\braket{\Psi_0\vert  \mathrm{e}^{+\mathrm{i} \hat{H} t}  \hat{A}_i^\dagger \mathrm{e}^{-\mathrm{i} \hat{H} t} \hat{B}_j\vert \Psi_0} \\
&=\braket{\Psi_0\vert  \hat{A}_i \mathrm{e}^{-\mathrm{i} (\hat{H}-E_0) t} \hat{B}_j\vert \Psi_0} 
```

On the first line, we have used operator 
$\hat{A}(t) = \mathrm{e}^{+\mathrm{i} t \hat{H}}\hat{A} \mathrm{e}^{-\mathrm{i} t \hat{H}}$
in the Heisenberg picture. In going from the second to the third line, we have used that
$\ket{\Psi_0}$ is the ground state of $\hat{H}$ with ground state energy $E_0$. Again, this
quantity will depend on the relative lattice vector connecting sites $i$ and $j$ in a
translation-invariant system. Once again, we can diagonalise the spatial dependence using a
multidimensional discrete Fourier transform. If we now furthermore also perform a Fourier
transform of the time-dependence into frequency space, we obtain the *dynamical structure
factor* given by

$$ S^{A,B}(\kappa, \omega) = \sum_{n} \delta(\omega- (E_n - E_0))\braket{\Psi_0 \vert \hat{A}(\kappa)^\dagger | \Psi_n}\braket{\Psi_n \ \hat{B}(k\kappa) \vert \Psi_0} $$

Here, $\hat{A}(k\kappa)$ and $\hat{B}(\kappa)$ correspond to the discrete Fourier transforms
of $\hat{A}_i$, which amounts to the momentum superposition. As a consequence,
$\hat{B}(\kappa) \ket{\Psi_0}$ is a state with definite momentum $k$ (provided the ground
state is translation invariant), and thus only has overlap with excited states
$\ket{\Psi_n}$ with momentum $\kappa$. Because of the factor $\delta(\omega - (E_n -E_0))$,
the dynamical structure factor $S^{A,B}(\kappa, \omega)$ can be nonzero only if there exist
eigenstates with momentum $\kappa$ and excitation energy $\omega$ in the spectrum of the
Hamiltonian $\hat{H}$. By studying $S^{A,B}(\kappa, \omega)$ for different choices of
operators $\hat{A}$ and $\hat{B}$, we can detect all eigenstates and map out the full
(low-energy) spectrum of $\hat{H}$.

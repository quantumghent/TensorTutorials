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

# Matrix Product Operators

If Matrix Product States are a tensor network way of representing quantum states in one
dimensions, we can similarly use tensor networks to represent the operators that act on
these states. Matrix Product Operators (MPOs) form a structured and convenient description
of such operators, that can capture most (if not all) relevant operators. Additionally, they
also form a natural way of representing the transfer matrix of a 2D statistical mechanical
system, and can even be used to study higher dimensional systems by mapping them to quasi-1D
systems. 

In this lecture, we will discuss the construction of MPOs, as well as showcase their use
through [MPSKit.jl](https://github.com/maartenvd/MPSKit.jl) and
[MPSKitModels.jl](https://github.com/maartenvd/MPSKitModels.jl).

```{code-cell}
using TensorKit
using MPSKit
using MPSKitModels
```

In general, an MPO is a chain of tensors, where each tensor has two physical indices and two
virtual indices:

<!-- Insert image MPO -->

## Statistical Mechanics in 2D

Before discussing one-dimensional transfer matrices, let us first consider how partition
functions of two-dimensional classical many-body systems can be naturally represented as a
tensor network. To this end, consider the partition function of the
[classical Ising model](https://en.wikipedia.org/wiki/Ising_model),

```{math}
\mathcal Z = \sum_{\{s_i\}} \text{e}^{-\beta H(\{s_i\})},
```

where $s_i$ denotes a configuration of spins, and $H(\{s_i\})$ is the corresponding
energy, as determined by the Hamiltonian:

```{math}
H(\{s_i\}) = -J \sum_{\langle i,j \rangle} s_i s_j
```

where the first sum is over nearest neighbors. 

### Partition Functions as Tensor Networks

As the expression for the partition function is an exponential of a sum, we can also write
it as a product of exponentials, which can be reduced to the following network:

```{image} /_static/figures/mpo/partition_function_1.svg
:scale: 12%
:name: partfunc
:align: center
```

Here, the black dots at the vertices represent Kronecker $\delta$-tensors,

```{image} /_static/figures/mpo/kronecker.svg
:scale: 12%
:name: kronecker
:align: center
```

and the matrices $t$ encode the Boltzmann weights associated to each nearest-neighbor interaction,

```{image} /_static/figures/mpo/boltzmann.svg
:scale: 12%
:name: boltzmann
:align: center
```

It is then simple, albeit somewhat involved to check that contracting this network gives
rise to the partition function, where the sum over all configurations is converted into the
summations in the contractions of the network. Finally, it is more common to absorb the edge
tensors into the vertex tensors by explicitly contracting them, such that the remaining
network consists of tensors at the vertices only:

```{image} /_static/figures/mpo/partition_function.svg
:scale: 12%
:name: partfunc
:align: center
```

````{note}
Because there are two edges per vertex, an intuitive way of absorbing the edge tensors is to
absorb for example the left and bottom edge tensors into the vertex tensor. However, this
leads to a slightly asymmetric form, and more commonly the square root $q$ of the Boltzmann
matrices is taken, such that each vertex tensor absorbs such a factor from each of the
edges, resulting in a rotation-invariant form.

```{image} /_static/figures/mpo/boltzmann_mpo.svg
:scale: 12%
:name: boltzmann_mpo
:align: center
```
````

In particular, the construction of the operator that makes up the MPO can be achieved in a
few lines of code, through the use of TensorKit:

```{code-cell} julia
β = 1.0

# construct edge tensors
t = TensorMap(ComplexF64[exp(β) exp(-β); exp(-β) exp(β)], ℂ^2, ℂ^2)
q = sqrt(t)

# construct vertex tensors
δ = TensorMap(zeros, ComplexF64, ℂ^2 ⊗ ℂ^2, ℂ^2 ⊗ ℂ^2)
δ[1, 1, 1, 1] = 1.0
δ[2, 2, 2, 2] = 1.0

# absorb edge tensors
@tensor O[-1 -2; -3 -4] := δ[1 2; 3 4] * q[-1; 1] * q[-2; 2] * q[3; -3] * q[4; -4]
```

### Transfer Matrices

In order to then evaluate the partition function, we can use the
[Transfer-matrix method](https://en.wikipedia.org/wiki/Transfer-matrix_method), which is a
technique that splits the two-dimensional network into rows (or columns) of so-called
transfer matrices, which are already represented as MPOs. In fact, this method has even led
to the famous exact solution of the two-dimensional Ising model by Onsager.
{cite}`onsager1944crystal`.

<!-- Insert image of transfer matrix -->

In the context of tensor networks, this technique is even useful beyond exactly solvable
cases, as efficient algorithms exist to determine the product of an MPO with an MPS in an
approximate manner. This allows us to efficiently split the computation of the partition
function in a sequence of one-dimensional contractions, thus reducing the complexity of the
problem by solving it layer by layer. 

### Thermodynamic Limit

Importantly, this technique is not limited to finite systems, and in fact allows for the
computation of the partition function of systems directly in the thermodynamic limit,
alleviating the need to consider finite-size effects and extrapolation techniques. The key
insight that allows for this is that the partition function may be written as

```{math}
\mathcal Z = \lim_{N \to \infty} \mathrm{Tr} \left( T^N \right)
```

where $T$ is the row-to-row transfer matrix, and $N$ is the number of rows (or columns) in
the network. If we then consider the spectral decomposition of the transfer matrix, we can
easily show that as the number of rows goes to infinity, the largest eigenvalue of the
transfer matrix dominates, and the partition function is given by

```{math}
\mathcal Z = \lim_{N \to \infty} \lambda_{\mathrm{max}}^N \braket{\psi}{\psi}
```

where $\lambda_{\mathrm{max}}$ is the largest eigenvalue of the transfer matrix, and
$\ket{\psi}$ is the corresponding (MPS) eigenvector. In other words, the partition function
can be computed if it is possible to find the largest eigenvalue of the transfer matrix, for
which efficient algorithms exist.

For example, one can resort to many types of _boundary MPS techniques_
{cite}`zauner-stauber2018variational`, which are a generic class of algorithms to
numerically solve these kinds of problems. In particular, they all rely on an efficient way
of finding an (approximate) solution to the following problem:

```{image} /_static/figures/alg/boundary_mps.svg
:scale: 12%
:name: boundary_mps
:align: center
```

### Expectation Values

In order to compute relevant quantities for such systems, we can verify that the expectation
value of an operator $O$ is given by the weighing the value of that operator for a given
microstate, with the probability of that microstate:

```{math}
\langle O \rangle = \frac{1}{\mathcal Z} \sum_{\{s_i\}} O(\{s_i\})\text{e}^{-\beta
H(\{s_i\})}
```

For a local operator $O_i$, this can again be written as a tensor network, where a single
Kronecker tensor at a vertex is replaced with a tensor measuring the operator, and then
absorbing the remaining edge tensors:

<!-- Insert image of expectation value -->

For example, in the case of the magnetisation $O = \sigma_z$, the tensor $M$ can be
explicitly constructed as follows:

```{code-cell} julia
Z = TensorMap(ComplexF64[1.0 0.0; 0.0 -1.0], ℂ^2, ℂ^2)
@tensor M[-1 -2; -3 -4] := δ[1 2; 3 4] * Z[4; 5] * q[-1; 1] * q[-2; 2] * q[3; -3] * q[5; -4]
```

Using this network, the expectation value can be computed by first contracting the top and
bottom part, replacing them by their fixed-point MPS representations, and then contracting
the remaining MPS-MPO-MPS sandwich. This is achieved by similarly contracting the left and
right part, replacing them by their fixed-point tensors, which are commonly called the
_environments_ $G_L$ and $G_R$, respectively. The final resulting network is then just a
local network, which can be contracted efficiently.

<!-- Insert image of expectation value with environments -->

```{note}
This process of sequentally reducing the dimensionality of the network can even be further
extended, where 3D systems can be studied by first determining a 2D boundary PEPS, for which
a 1D boundary MPS can be determined, which admits 0D boundary tensors. This kind of
algorithms are commonly referred to as _boundary methods_.
```

## Quantum Mechanics in 1+1D

For quantum systems in one spatial dimension, the construction of MPOs boils down to the
ability to write a sum of local operators in MPO-form. The resulting operator has a very
specific structure, and is often referred to as a _Jordan block MPO_.

### Jordan Block MPOs

For example, if we consider the
[Transverse-field Ising model](https://en.wikipedia.org/wiki/Transverse-field_Ising_model),

```{math}
H = -J \sum X_j X_{j+1} - h \sum Z_j
```

it can be represented as an MPO through the (operator-valued) matrix, 

```{math}
W = \begin{pmatrix}
1 & X & -hZ \\ 
0 & 0 & -JX \\
0 & 0 & 1
\end{pmatrix}
```

along with the boundary vectors,

```{math}
v_L = \begin{pmatrix}
1 & 0 & 0
\end{pmatrix}
, \qquad 
v_R = \begin{pmatrix}
0 \\ 0 \\ 1
\end{pmatrix}
```

The Hamiltonian on $N$ sites is then given by the contraction

```{math}
H = V_L W^{\otimes N} V_R
```

```{note}
While the above example can be constructed from building blocks that are strictly local
operators, this is not always the case, especially when symmetries are involved. In those
cases, the elements of the matrix $W$ have additional virtual legs that are contracted
between different sites.
```

### Finite-State Machines

An intuitive approach to construct such MPOs is to consider the sum of local
terms by virtue of a
[finite-state machine](https://en.wikipedia.org/wiki/Finite-state_machine). This is a
mathematical model of computation that consists of a finite set of states, and a set of
transitions between those states. In the context of MPOs, this is realised by associating
each _virtual level_ with a state, and each transition then corresponds to applying a local
operator. In that regard, the MPO is then a representation of the state of the finite-state
machine, and the matrix $W$ is the transition matrix of the machine.

In general, the matrix $W$ can then be thought of as a block matrix with entries

```{math}
\begin{pmatrix}
1 & C & D \\
0 & A & B \\
0 & 0 & 1
\end{pmatrix}
```

which corresponds to the finite-state diagram:

<!-- Insert image of finite-state diagram -->

It can then be shown that this MPO generates all single-site local operators $D$, two-site
operators $CB$, three-site operators $CAB$, and so on. In other words, the MPO is a
representation of the sum of all local operators, and by carefully extending the structure
of the blocks $A$, $B$, $C$, and $D$, it is possible to construct MPOs that represent sums
of generic local terms, and even approximate long-range interactions by a sum of
exponentials.

To gain a bit more understanding of this, we can use the following code to reconstruct the
total sum of local terms, starting from the Jordan MPO construction:

```{code-cell} julia
using Symbolics

L = 4
# generate W matrices
@variables A[1:L] B[1:L] C[1:L] D[1:L]
Ws = map(1:L) do l
    return [1 C[l] D[l]
            0 A[l] B[l]
            0 0    1]
end

# generate boundary vectors
Vₗ = [1, 0, 0]'
Vᵣ = [0, 0, 1]

# expand the MPO
expand(Vₗ * prod(Ws) * Vᵣ)
```

### Expectation Values

In order to compute expectation values of such MPOs, we can use the same technique as
before, and sandwich the MPO between two MPSs.

<!-- Insert image of expectation value -->

However, care must be taken when the goal is to determine a local expectation value density,
as this is not necessarily well-defined. In fact, the MPO represents the sum of all local
terms, and sandwiching it will always lead to the total energy. In order to consistently
define local contributions, a choice must be made how to *distribute* this among the sites.
For example, even in the case of two-site local operators, it is unclear if this local
expectation value should be accredited to the left, or right site, or even split between
both sites. In the implementation of MPSKit, the chosen convention is to distribute the
expectation value evenly among its starting and ending point, in order to not overcount
contributions of long-range interactions.

Typically this is achieved by renormalizing the environment tensors in a particular way,
such that then local expectation values can be obtained by either contracting the first row
of $W$ with the right regularized environment, or the last column of $W$ with the left
regularized environment. This respectively yields the expectation value of all terms
starting at that site, or all terms ending at that site.

Again, it can prove instructive to write this out explicitly for some small examples to gain
some intuition. Doing this programatically, we get all terms starting at some site as
follows:

```{code-cell} julia
Ws_reg_right = Ws .- Ref([1 0 0; 0 0 0; 0 0 0])
expand(Vₗ * Ws_reg_right[end-2] * Ws_reg_right[end-1] * Ws_reg_right[end] * Vᵣ)
```

and similarly all terms ending at some site as follows:

```{code-cell} julia
Ws_reg_left = Ws .- Ref([0 0 0; 0 0 0; 0 0 1])
expand(Vₗ * Ws_reg_left[1] * Ws_reg_left[2] * Ws_reg_left[3] * Vᵣ)
```


### Jordan MPOs in the Thermodynamic Limit

In the thermodynamic limit, the same MPO construction can be used to represent the infinite
sum of local terms. However, special care must be taken when considering expectation values,
as now only local expectation values are well-defined, and the total energy diverges with
the system size.

This is achieved by considering the same regularization of the environment tensors, such
that the divergent parts are automatically removed. This construction can be found in more
detail in {cite}`hubig17generic`.

### Quasi-1D Systems

Finally, it is worth noting that the MPO construction can also be used to study
two-dimensional systems, by mapping them to quasi-one-dimensional systems. This is typically
achieved by imposing periodic boundary conditions in one of the spatial directions, and then
_snaking_ an MPS through the resulting lattice. In effect, this leads to a one-dimensional
model with longer-range interactions, which can then be studied using the standard MPS
techniques. However, the
[no free lunch theorem](https://en.wikipedia.org/wiki/No_free_lunch_theorem) applies here as
well, and the resulting model will typically require a bond dimension that grows
exponentially with the periodic system size, in order to achieve the area law of
entanglement in two-dimensional systems.


### MPSKitModels and the `@mpoham` Macro

While the above construction of MPOs is quite general, it is also quite cumbersome to
manually construct, especially when dealing with complicated lattices or non-trivial unit
cells. To this end, the package
[MPSKitModels.jl](https://github.com/maartenvd/MPSKitModels.jl) offers a convenient way of
constructing these MPOs automatically, by virtue of the `@mpoham` macro. This macro allows
for the construction of MPOs by specifying the local operators that are present in the
Hamiltonian, and the lattice on which they act. For example, we can construct the MPO for
the Heisenberg models with nearest- or next-nearest-neighbor interactions as follows:

```{code-cell} julia
:tags: ["hide-output"]
J₁ = 1.2
SS = S_exchange() # predefined operator in MPSKitModels

lattice = InfiniteChain(1)
H₁ = @mpoham begin
    sum(J₁ * SS{i, j} for (i, j) in nearest_neighbours(lattice))
end
```

```{code-cell} julia
:tags: ["hide-output"]
lattice = InfiniteCylinder(4)
H₂ = @mpoham begin
    sum(J₁ * SS{i, j} for (i, j) in nearest_neighbours(lattice))
end
```

```{code-cell} julia
:tags: ["hide-output"]
J₂ = 0.8
lattice = InfiniteCylinder(4)
H₃ = @mpoham begin
    sum(J₁ * SS{i, j} for (i, j) in nearest_neighbours(lattice)) +
    sum(J₂ * SS{i, j} for (i, j) in next_nearest_neighbours(lattice))
end
```

## Conclusion

In conclusion, Matrix Product Operators are a powerful tool to represent quantum operators
as well as transfer matrices. They allow for efficient and versatile expressions of
expectation values, and form the building block for many tensor network algorithms, both in
(1+1) or (2+0) dimensions, as well as in higher dimensions.

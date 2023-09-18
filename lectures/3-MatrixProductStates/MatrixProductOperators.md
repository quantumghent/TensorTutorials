# Matrix Product Operators

If Matrix Product States are a tensor network way of representing quantum states in one
dimensions, we can similarly use tensor networks to represent the operators that act on
those states. Matrix Product Operators (MPOs) form a structured and convenient description
of such operators, that can capture most (if not all) relevant operators. Additionally, they
also form a natural way of representing the transfer matrix of a 2D statistical mechanical
system, and can even be used to study higher dimensional systems by mapping them to quasi-1D
systems.

In general, an MPO is a chain of tensors, where each tensor has two physical indices and two
virtual indices:

<!-- Insert image MPO -->

## Statistical Mechanics in 2D

Before discussing one-dimensional transfer matrices, let us first consider how partition
functions of two-dimensional classical many-body systems can be naturally represented as
a tensor network. To this end, consider the partition function of the [classical Ising model](https://en.wikipedia.org/wiki/Ising_model),

```{math}
Z(\beta) = \sum_\sigma e^{-\beta H(\sigma)}
```

where $\sigma$ is a configuration of spins, and $H(\sigma)$ is the corresponding energy, as determined by the Hamiltonian:

```{math}
H(σ) = -J \sum_{\langle i,j \rangle} \sigma_i \sigma_j - h \sum_i \sigma_i
```

where the first sum is over nearest neighbors, and the second sum is over all sites. 

### Partition Functions as Tensor Networks

As the expression for the partition function is an exponential of a sum, we can also write
it as a product of exponentials, and it is an easy exercise to show that this is equivalent
to contracting the following network, where the tensors at the vertices are Kronecker
deltas, which denote the spins at each site, and the tensors on the edges are exponentials
of the local Hamiltonian that represent the interactions between spins:

<!-- Insert image of tensor network -->

It is then easy to check that contracting this network gives rise to the partition function,
where the sum over all configurations is converted into the summations in the contractions
of the network. Finally, it is more common to absorb the edge tensors into the vertex
tensors by explicitly contracting them, such that the remaining network consists of tensors
at the vertices only:

<!-- Insert image of tensor network -->

```{note}
Because there are two edges per vertex, an intuitive way of absorbing the edge tensors is to
absorb for example the left and bottom edge tensors into the vertex tensor. However, this
leads to a slightly asymmetric form, and more commonly the square root of the edge tensors
is taken, such that each vertex tensor absorbs such a factor from each of the edges,
resulting in a rotation-symmetric form.
```

### Transfer Matrices

In order to then evaluate the partition function, we can use the
[Transfer-matrix method](https://en.wikipedia.org/wiki/Transfer-matrix_method), which is a
technique that splits the two-dimensional network into rows (or columns) of so-called
transfer matrices, which are easily seen to be represented as MPOs. In fact, this method has
led to the famous exact solution of the two-dimensional Ising model by Onsager.
{cite}`Onsager1944`.

In the context of tensor networks, this technique is particularly useful as efficient
algorithms exist to determine the product of an MPO with an MPS, either exactly or
approximately. This allows us to efficiently split the computation of the partition function
in a sequence of one-dimensional contractions, thus reducing the complexity of the problem
by solving it layer by layer.

### Thermodynamic Limit

Importantly, this technique is not limited to finite systems, and in fact allows for the
computation of the partition function of systems directly in the thermodynamic limit,
alleviating the need to consider finite-size effects and extrapolation techniques. The key
insight that allows for this is that the partition function may be written as

```{math}
Z = \lim_{N \to \infty} \mathrm{Tr} \left( T^N \right)
```

where $T$ is the transfer matrix, and $N$ is the number of rows (or columns) in the network.
If we then consider the spectral decomposition of the transfer matrix, we can easily show
that as the number of rows goes to infinity, the largest eigenvalue of the transfer matrix
dominates, and the partition function is given by

```{math}
Z = \lim_{N \to \infty} \lambda_{\mathrm{max}}^N
```

where $\lambda_{\mathrm{max}}$ is the largest eigenvalue of the transfer matrix. In other
words, the partition function can be computed if it is possible to find the largest
eigenvalue of the transfer matrix, for which efficient algorithms exist.

### Expectation Values

In order to compute relevant quantities for such systems, we can verify that the expectation
value of an operator $O$ is given by the weighing the value of that operator for a given
microstate, with the probability of that microstate:

```{math}
\langle O \rangle = \frac{1}{Z} \sum_\sigma O(\sigma) e^{-\beta H(\sigma)}
```

Again, this can be written as a tensor network, where the network is the same as before, but
with a single vertex tensor exchanged for one that represents the operator $O$:

<!-- Insert image of expectation value -->

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
specific structure, and is often referred to as a "Jordan block" MPO.

### Jordan Block MPOs

For example, if we consider the
[Transverse-field Ising model](https://en.wikipedia.org/wiki/Transverse-field_Ising_model),

```{math}
H = -J \sum X_j X_{j+1} - h \sum Z_j
```

it can be represented as an MPO through the (operator-valued) matrix, 

```{math}
M = \begin{pmatrix}
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
H = V_L M^{\otimes N} V_R
```

```{note}
While the above example can be constructed from building blocks that are just local
operators, this is not always the case, especially when symmetries are involved. In those
cases, the elements of the matrix $M$ have additional virtual legs that are contracted
between different sites.
```

### Finite-State Machines

In fact, an intuitive approach of constructing such MPOs is to consider the sum of local
terms by virtue of a
[finite-state machine](https://en.wikipedia.org/wiki/Finite-state_machine). This is a
mathematical model of computation that consists of a finite set of states, and a set of
transitions between those states. In the context of MPOs, this is realised by associating
each _virtual level_ with a state, and each transition then corresponds to applying a local
operator. In that regard, the MPO is then a representation of the state of the finite-state
machine, and the matrix $M$ is the transition matrix of the machine.

In general, the matrix $M$ can then be thought of as a block matrix with entries

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

### Expectation Values

In order to compute expectation values of such MPOs, we can use the same technique as
before, and sandwich the MPO between two MPSs.

<!-- Insert image of expectation value -->

However, care must be taken when the goal is to determine a local expectation value density,
as this is not necessarily well-defined. In fact, the MPO-MPS sandwich will always lead to
the total energy, and in order to consistently define local contributions, a choice must be
made how to *distribute* this among the sites. For example, even in the case of two-site
local operators, it is unclear if this local expectation value should be accredited to the
left, right, or both sites. In the implementation of MPSKit, the chosen convention is to
distribute the expectation value evenly among its starting and ending point, in order to not
overcount contributions of long-range interactions.

Thus, the computation of local expectation values is done by first contracting the left and
right ends of the network, and denoting the resulting tensors as $G_L$ and $G_R$,
respectively. These tensors are typically referred to as the _environments_ of the network,
as they depict the part of the network that is not under consideration, and that is thus
considered as constant.

<!-- Insert image of environments -->

Then, the local expectation value is given by the mean of the contraction with the first row
of $M$, which is equivalent to all terms starting at that site, and the last column of $M$,
which is equivalent to all terms ending at that site.

<!-- Insert image of local expectation value -->

Again, it is instructive to write this out explicitly for some small examples to gain some
intuition.

### Jordan MPOs in the Thermodynamic Limit

In the thermodynamic limit, the same MPO construction can be used to represent the infinite
sum of local terms. However, special care must be taken when considering expectation values,
as now only local expectation values are well-defined, and the total energy diverges with
the system size.

This is achieved by considering a regularization of the environment tensors, such that the
divergent parts are already removed. This construction can be found in more detail in
{cite}`mcculloch(?)`.

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


## Conclusion

In conclusion, Matrix Product Operators are a powerful tool to represent quantum operators
as well as transfer matrices. They allow for efficient and versatile expressions of
expectation values, and form the building block for many tensor network algorithms, both in
(1+1) or (2+0) dimensions as well as even beyond.

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

# A Symmetric Tensor Deep Dive: Constructing Your First Tensor Map

In this tutorial, we will demonstrate how to construct specific `TensorMap`s which are relevant to some common physical systems, with an increasing degree of complexity. We will assume the reader has gone through the tutorial sections on [tensor network theory](tensor_networks) and [symmetries in tensor networks](symmetric_tensors). In going through these examples we aim to provide a relatively gently introduction to the meaning of [symmetry sectors](https://jutho.github.io/TensorKit.jl/latest/man/sectors/#ss_sectors) and [vector spaces](https://jutho.github.io/TensorKit.jl/latest/man/sectors/#ss_rep) within the context of [TensorKit.jl](https://jutho.github.io/TensorKit.jl/latest/), [how to initialize a `TensorMap` over a given vector space](https://jutho.github.io/TensorKit.jl/latest/man/tensors/#ss_tensor_construction) and finally how to manually set the data of a symmetric `TensorMap`. We will keep our discussion as intuitive and simple as possible, only adding as many technical details as strictly necessary to understand each example. When considering a different physical system of interest, you should then be able to adatpt these recipes and the intuition behind them to your specific problem at hand.

```{note}
Many of these examples are already implemented in the [MPSKitModels.jl package](https://github.com/maartenvd/MPSKitModels.jl), in which case we basically provide a narrated walk-through of the corresponding code.
```

```{code-cell} julia
using LinearAlgebra
using TensorKit
using MPSKitModels
using WignerSymbols
using SUNRepresentations
using Test # for showcase testing
```


## Level 0: The 1-Dimensional Ising Model

As the most basic example, we will consider the [1-dimensional transverse-field Ising model](https://en.wikipedia.org/wiki/Transverse-field_Ising_model), whose Hamiltonian is given by

```{math}
:label: eq:isingham
H = -J \left (\sum_{\langle i, j \rangle} Z_i Z_j + g \sum_{i} X_i\right).
```

Here, $X_i$ and $Z_i$ are the [Pauli operators](https://en.wikipedia.org/wiki/Pauli_matrices) acting on site $i$, and the first sum runs over pairs of nearest neighbors $\langle i, j \rangle$. This model has a global $\mathbb{Z}_2$ symmetry, as it is invariant under the transformation $U H U^\dagger = H$ where the symmetry transformation $U$ is given by a global spin flip,

```{math}
:label: eq:z2gen
U = \prod_i X_i.
```

We will circle back to the implications of this symmetry later.

As a warmup we will implement the Hamiltonian {eq}`eq:isingham` in the standard way by encoding the matrix elements of the single-site operators $X$ and $Z$ into aan array of complex numbers, and then combine them in a suitable way to get the Hamiltonian terms. Instead of using plain Julia arrays, we will use a representation in terms of `TensorMap`s over complex vector spaces. These will essentially just be wrappers around base arrays at this point, but their construction requires some consideration of the notion of *spaces*, which generalize the notion of `size` for arrays. Each of the operators $X$ and $Z$ acts on a local 2-dimensional complex vector space. In the context of TensorKit.jl such a space can be represented as `ComplexSpace(2)`, or using the convenient shorthand `ℂ^2`. A single-site Pauli operator maps from a domain physical space to a codomain physical space, and can therefore be represented as instances of a `TensorMap(..., ℂ^2 ← ℂ^2)`. The corresponding data can then be filled in by hand according to the familiar Pauli matrices in the following way:

```{code-cell} julia
:tags: [hide-output]
# initialize numerical data for Pauli matrices
x_mat = ComplexF64[0 1; 1 0]
z_mat = ComplexF64[1 0; 0 -1]

# construct physical Hilbert space
V = ℂ^2

# construct the physical operators as TensorMaps
X = TensorMap(x_mat, V ← V)
Z = TensorMap(z_mat, V ← V)

# combine single-site operators into two-site operator
ZZ = Z ⊗ Z
```

We can easily verify that our operators have the desired form by checking their data in the computational basis:

```{code-cell} julia
ZZ
```

```{code-cell} julia
X
```

```{note}
In order to combine these local operators into a concrete Hamiltonian that can be used in [MPSKit.jl](https://github.com/maartenvd/MPSKit.jl) we can make use of the convenient `@mpoham` macro exported by [MPSKitModels.jl](https://github.com/maartenvd/MPSKitModels.jl). For an infinite translation invariant Ising chain, we can use the following piece of code which produces the Hamiltonian in an interesting-looking form (see MPSKit.jl for details on this format).
```
```{code-cell} julia
:tags: [hide-output]
lattice = InfiniteChain(1)
H = @mpoham begin
    sum(nearest_neighbours(lattice)) do (i, j)
        return ZZ{i,j}
    end + sum(vertices(lattice)) do i
        return X{i}
    end
end
```


## Level 1: The $\mathbb{Z}_2$-symmetric Ising Model

### The irrep basis and block sparsity

Let us now return to the global $\mathbb{Z}_2$ invariance of the Hamiltonian {eq}`eq:isingham`, and consider what this implies for its local terms $ZZ$ and $X$. Representing these operators as `TensorMap`s, the invariance of $H$ under a global $\mathbb{Z}_2$ transformation implies the following identities for the local tensors:

```{figure} ../_static/SymmetricTensors/ZZX_symm.svg
:scale: 12%
:name: ZZX_symm
:align: center
```

Recalling the [discussion on symmetries in tensor networks](symmetric_tensors), we recognize that these identitities precisely mean that these local tensors transform trivially under a tensor product representation of $\mathbb{Z}_2$. This implies that, in an appropriate basis for the local physical vector space, our local tensors would become block-diagonal where each so-called *matrix block* is labeled by a $\mathbb{Z}_2$ irrep. From the same discussion, we recall that the appropriate local basis transformation is precisely the one that brings the local representation $X$ into block-diagonal form. Clearly, this transformation is nothing more than the Hadamard transformation which maps the computational basis of $Z$ eigenstates $\{\ket{\uparrow}, \ket{\downarrow}\}$ to that of the $X$ eigenstates $\{\ket{+}, \ket{-}\}$ defined as $\ket{+} = \frac{\ket{\uparrow} + \ket{\downarrow}}{\sqrt{2}}$ and $\ket{-} = \frac{\ket{\uparrow} - \ket{\downarrow}}{\sqrt{2}}$. In the current context, this basis is referred to as the *irrep basis* of $\mathbb{Z}_2$, where the local basis state $\ket{+}$ corresponds to the trivial representation of $\mathbb{Z}_2$ while $\ket{-}$ corresponds to the sign representation.

Next, let's make the statement that 'the matrix blocks of the local tensors are labeled by $\mathbb{Z}_2$ irreps' more concrete. To this end, consider the action of $ZZ$ in the irrep basis, which is given by the four nonzero matrix elements

```{math}
\begin{align}
ZZ : \mathbb C^2 \otimes \mathbb C^2 &\to \mathbb C^2 \otimes \mathbb C^2 : \\
\ket{+} \otimes \ket{+} &\mapsto \ket{-} \otimes \ket{-},  \\
\ket{+} \otimes \ket{-} &\mapsto \ket{-} \otimes \ket{+},  \\
\ket{-} \otimes \ket{+} &\mapsto \ket{+} \otimes \ket{-},  \\
\ket{-} \otimes \ket{-} &\mapsto \ket{+} \otimes \ket{+}.
\end{align}
```

If we denote the trivial $\mathbb{Z}_2$ irrep by $'0'$, corresponding to a local $\ket{+}$ state, and the sign irrep by $'1'$, corresponding to a local $\ket{-}$ state, and recall that in this notation the fusion rules of $\mathbb{Z}_2$ are given by addition modulo 2, we can associate each of the above matrix elements to a so-called *fusion tree* of $\mathbb{Z}_2$ irreps with a corresponding coefficient of 1,

```{figure} ../_static/SymmetricTensors/Z2_fusiontrees.svg
:scale: 12%
:name: Z2_fusiontrees
:align: center
```

From this we can observe our previous statement very clearly: the $ZZ$ operator indeed consists of two distinct two-dimensional matrix blocks, each of which are labeled by the value of the *coupled irrep* on the middle line of each fusion tree. The first block corresponds to the even coupled irrep '0', and acts within the two-dimensional subspace spanned by $\{\ket{+,+}, \ket{-,-}\}$, while the second block corresponds to the odd coupled irrep '1', and acts within the two-dimensional subspace spanned by $\{\ket{+,-}, \ket{-,+}\}$. In TensorKit.jl, this block-diagonal structure of a symmetric tensor is explicitly encoded into its representation as a `TensorMap`, where only the matrix blocks corresponding to each coupled irrep are stored.

For our current purposes however, we never really need to explicitly consider these matrix blocks. Indeed, when constructing a `TensorMap` it is sufficient to set its data by manually assigning a matrix element to each [fusion tree of the form above](Z2_fusiontrees) labeled by a given tensor product of irreps. This matrix element is then automatically inserted into the appropriate matrix block. So, for the purpose of this tutorial **we will interpret a symmetric `TensorMap` simply as a list of fusion trees, to each of which corresponds a certain reduced matrix element**.

```{note}
In general, such a reduced matrix element is not necessarily a scalar, but rather an array whose size is determined by the degeneracy of the irreps in the codomain and domain of the fusion tree. For this reason, a reduced matrix element associated to a given fusion tree is also referred to as an *array block*. In the following we will use terms 'reduced matrix element', 'array block' or just 'block' interchangeably. However, it should be remembered that these are distinct from the matrix blocks in the block-diagonal decomposition of the tensor.
```


(fusion_trees)=
### Fusion trees and how to use them

This view of the underlying symmetry structure in terms of fusion trees and corresponding array blocks is a very convenient way of working with the `TensorMap` type. Consider a generic fusion tree of the form

```{figure} ../_static/SymmetricTensors/fusiontree.svg
:scale: 12%
:name: fusiontree
:align: center
```

which can be used to label a block of a `TensorMap` corresponding to a two-site operator. This object should actually be seen as a *pair of fusion trees*. The first member of the pair, related to the codomain of the `TensorMap`, is referred to as the *splitting tree* and encodes how the *coupled charge* $c$ splits into the *uncoupled charges* $s_1$ and $s_2$. The second member of the pair, related to the domain of the `TensorMap`, is referred to as the *fusion tree* and encodes how the uncoupled charges $f_1$ and $f_2$ fuse to the coupled charge $c$. Both the splitting and fusion tree can be represented as a [`TensorKit.FusionTree`](https://jutho.github.io/TensorKit.jl/latest/lib/sectors/#TensorKit.FusionTree) instance. You will find such a `FusionTree` has the following properties encoded into its fields:

- `uncoupled::NTuple{N,I}`: a list of `N` uncoupled charges of type `I<:Sector`
- `coupled::I`: a single coupled charge of type `I<:Sector`
- `isdual::NTuple{N,Bool}`: a list of booleans indicating whether the corresponding uncoupled charge is dual
- `innerlines::NTuple{M,I}`: a list of inner lines of type `I<:Sector` of length `M = N - 2`
- `vertices::NTuple{L,T}`: list of fusion vertex labels of type `T` and length `L = N - 1`

For our current application only `uncoupled` and `coupled` are relevant, since $\mathbb{Z}_2$ irreps are self-dual and have Abelian fusion rules. We will come back to these other properties when discussion more involved applications. Given some `TensorMap`, the method `TensorKit.fusiontrees(t::TensorMap)` returns an iterator over all pairs of splitting and fusion trees that label the blocks of `t`.

We can now put this into practice by directly constructing the $ZZ$ operator in the irrep basis as a $\mathbb{Z}_2$-symmetric `TensorMap`. We will do this in three steps:

- First we construct the physical space at each site as a $\mathbb{Z}_2$-graded vector space.
- Then we initialize an empty `TensorMap` with the correct domain and codomain vector spaces built from the previously constructed physical space.
- And finally we iterate over all splitting and fusion tree pairs and manually fill in the corresponding nonzero blocks of the operator.

After the basis transform to the irrep basis, we can view the two-dimensional complex physical vector space we started with as being spanned by the trivial and sign irrep of $\mathbb{Z}_2$. In the language of TensorKit.jl, this can be implemented as a `Z2Space`, an alias for a $\mathbb{Z}_2$-[graded vector space](https://jutho.github.io/TensorKit.jl/latest/lib/spaces/#TensorKit.GradedSpace) `Vect[Z2Irrep]`, which contains the trivial irrep `Z2Irrep(0)` with degeneracy 1 and the sign irrep `Z2Irrep(1)` with degeneracy 1. We can define this space in the following way and check its dimension:

```{code-cell} julia
V = Z2Space(0 => 1, 1 => 1)
dim(V)
```

Given this physical space, we can initialize the $ZZ$ operator as an empty `TensorMap` with the appropriate structure.

```{code-cell} julia
ZZ = TensorMap(zeros, ComplexF64, V ⊗ V ← V ⊗ V)
```

The output of this command again demonstrates the underlying structure of a symmetric tensor. We see that all eight valid fusion trees with two incoming irreps and two outgoing irreps [of the type above](fusiontree) are listed with their corresponding block data. Each of these blocks is an array of shape $(1, 1, 1, 1)$ since each irrep occuring in the space $V$ has degeneracy 1. Using the `fusiontrees` method and the fact that we can index a `TensorMap` using a splitting/fusion tree pair, we can now fill in the nonzero blocks of the operator by observing that the $ZZ$ operator flips the irreps of the uncoupled charges in the domain with respect to the codomain, as shown in the diagrams above. Flipping a given `Z2Irrep` in the codomain can be implemented by fusing them with the sign irrep `Z2Irrep(1)`, giving:

```{code-cell} julia
flip_charge(charge::Z2Irrep) = only(charge ⊗ Z2Irrep(1))
for (s, f) in fusiontrees(ZZ)
    if s.uncoupled == map(flip_charge, f.uncoupled)
        ZZ[s, f] .= 1
    end
end
ZZ
```

Indeed, the resulting `TensorMap` exactly encodes the matrix elements of the $ZZ$ operator shown in {ref}`the diagrams above <Z2_fusiontrees>`. The $X$ operator can be constructed in a similar way. Since it is by definition diagonal in the irrep basis with blocks directly corresponding to the trivial and sign irrep, its construction is particularly simple:

```{code-cell} julia
X = TensorMap(zeros, ComplexF64, V ← V)
for (s, f) in fusiontrees(X)
    if only(f.uncoupled) == Z2Irrep(0)
        X[s, f] .= 1
    else
        X[s, f] .= -1
    end
end
X
```

Given these local operators, we can use them to construct the full manifestly $\mathbb{Z}_2$-symmetric Hamiltonian.

```{note}
An important observation is that when explicitly imposing the $\mathbb{Z}_2$ symmetry we directly constructed the full $ZZ$ operator as a single symmetric tensor. This in contrast to the case without symmetries, where we constructed a single-site $Z$ operator and then combined them into a two-site operator. Clearly this can no longer be done when imposing $\mathbb{Z}_2$, since a single $Z$ is not invariant under conjugation with the symmetry operator $X$. One might wonder whether it is still possible to construct a two-site Hamiltonian term by combining local objects. This is possible if one introduces an auxiliary index on the local tensors that carries a non-trivial charge. The intuition behind this will become more clear in the next example.
```


## Level 2: The $\mathrm{U}(1)$ Bose-Hubbard Model

For our next example we will consider the [Bose-Hubbard model](https://en.wikipedia.org/wiki/Bose%E2%80%93Hubbard_model), which describes interacting bosons on a lattice. The Hamiltonian of this model is given by
```{math}
:label: eq:bhh
H = -t \sum_{\langle i,j \rangle} \left( a_{i}^+ a_{j}^- + a_{i}^- a_{j}^+ \right) - \mu \sum_i N_i + \frac{U}{2} \sum_i N_i(N_i - 1).
```
This Hamiltonian is defined on the [Fock space](fock_space) associated to a chain of bosons, where the action bosonic creation, annihilation and number operators $a^+$, $a^-$ and $N = a^+ a^-$ in the local occupation number basis is given by
```{math}
:label: eq:bosonopmatel
\begin{align}
a^+ \ket{n} &= \sqrt{n + 1} \ket{n + 1} \\
a^- \ket{n} &= \sqrt{n} \ket{n - 1} \\
N \ket{n} &= n \ket{n}
\end{align}
```
Their bosonic nature can be summarized by the familiar the commutation relations
```{math}
\begin{align*}
\left[a_i^-, a_j^-\right] &= \left[a_i^+, a_j^+\right] = 0 \\
\left[a_i^-, a_j^+\right] &= \delta_{ij} \\
\left[N, a^+\right] &= a^+ \\
\left[N, a^-\right] &= -a^- \\
\end{align*}
```

This Hamiltonian is invariant under conjugation by the global particle number operator, $U H U^\dagger = H$, where
```{math}
U = \sum_i N_i
```
This invariance corresponds to a $\mathrm{U}(1)$ particle number symmetry, which can again be manifestly imposed when constructing the Hamiltonian terms as `TensorMap`s. From the representation theory of $\mathrm{U}(1)$ we know that it's irreps are all one-dimensional and can be labeled by integers, where the fusion of two irreps is given by addition.


### Directly Constructing the Hamiltonian Terms

We recall from our discussion on the $\mathbb{Z}_2$ symmetric Ising model that, in order to construct the Hamiltonian terms as symmetric tensors, we should work in the irrep basis where the symmetry transformation is block diagonal. In the current case, the symmetry operation is the particle number operator, which is already diagonal in the occupation number basis. Therefore, we don't need an additional local basis transformation this time, and can just observe that each local basis state can be identified with the $\mathrm{U}(1)$ irrep associated to the corresponding occupation number.

Following the same approach as before, we first write down the action of the Hamiltonian terms in the irrep basis:

```{math}
\begin{align*}
a_i^+ a_j^- \ket{n_i, n_j} &= \sqrt{(n_i + 1)n_j} \ket{n_i + 1, n_j - 1} \\
a_i^- a_j^+ \ket{n_i, n_j} &= \sqrt{n_i(n_j + 1)} \ket{n_i - 1, n_j + 1} \\
N \ket{n} &= n \ket{n}
\end{align*}
```

It is then a simple observation that these matrix elements are exactly captured by the following $\mathrm{U}(1)$ fusion trees with corresponding block values:

```{figure} ../_static/SymmetricTensors/U1_fusiontrees.svg
:scale: 12%
:name: U1_fusiontrees
```

This gives us all the information necessary to construct the corresponding `TensorMap`s. We follow the same steps as outlined in the previous example, starting with the construction of the physical space. This will now be a $\mathrm{U}(1)$ graded vector space `U1Space`, where each basis state $\ket{n}$ in the occupation number basis is represented by the corresponding $\mathrm{U}(1)$ irrep `U1Irrep(n)` with degeneracy 1. While this physical space is in principle infinite dimensional, we will impose a cutoff in occupation number at a maximum of 5 bosons per site, giving a 6-dimensional vector space:

```{code-cell} julia
cutoff = 5
V = U1Space(n => 1 for n in 0:cutoff)
```

We can now initialize the $a^+ a^-$, $a^- a^+$ and $N$ operators as empty `TensorMap`s with the correct domain and codomain vector spaces, and fill in the nonzero blocks associated to [the fusion trees shown above](U1_fusiontrees). To do this we need access to the integer label of the $\mathrm{U}(1)$ irreps in the fusion and splitting trees, which can be accessed through the `charge` field of the `U1Irrep` type.

```{code-cell} julia
:tags: [hide-output]
a⁺a⁻ = TensorMap(zeros, ComplexF64, V ⊗ V ← V ⊗ V)
for (s, f) in fusiontrees(a⁺a⁻)
    if s.uncoupled[1] == only(f.uncoupled[1] ⊗ U1Irrep(1)) && s.uncoupled[2] == only(f.uncoupled[2] ⊗ U1Irrep(-1))
        a⁺a⁻[s, f] .= sqrt(s.uncoupled[1].charge * f.uncoupled[2].charge)
    end
end
a⁺a⁻
```

```{code-cell} julia
:tags: [hide-output]
a⁻a⁺ = TensorMap(zeros, ComplexF64, V ⊗ V ← V ⊗ V)
for (s, f) in fusiontrees(a⁻a⁺)
    if s.uncoupled[1] == only(f.uncoupled[1] ⊗ U1Irrep(-1)) && s.uncoupled[2] == only(f.uncoupled[2] ⊗ U1Irrep(1))
        a⁻a⁺[s, f] .= sqrt(f.uncoupled[1].charge * s.uncoupled[2].charge)
    end
end
a⁻a⁺
```

```{code-cell} julia
:tags: [hide-output]
N = TensorMap(zeros, ComplexF64, V ← V)
for (s, f) in fusiontrees(N)
    N[s, f] .= f.uncoupled[1].charge
end
N
```


### Creation and Annihilation Operators as Symmetric Tensors

Just as in the $\mathbb{Z}_2$ case, it is obvious that we cannot directly construct the creation and annihilation operators as instances of a `TensorMap(..., V ← V)` since they are not invariant under conjugation by the symmetry operator. However, it is possible to construct them as `TensorMap`s using an *auxiliary vector space*, based on the following intuition. The creation operator $a^+$ violates particle number conservation by mapping the occupation number $n$ to $n + 1$. From the point of view of representation theory, this process can be thought af the *fusion* of an `U1Irrep(n)` with an `U1Irrep(1)`, naturally giving the fusion product `U1Irrep(n + 1)`. This means we can represent $a^+$ as a `TensorMap(..., V ← V ⊗ A)`, where the auxiliary vector space `A` is contains the $+1$ irrep with degeneracy 1,  `A = U1Space(1 => 1)`. Similarly, the decrease in occupation number when acting with $a^-$ can be thought of as the *splitting* of an `U1Irrep(n)` into an `U1Irrep(n - 1)` and an `U1Irrep(1)`, leading to a representation in terms of a `TensorMap(..., A ⊗ V ← V)`. Based on these observations, we can represent the matrix elements {eq}`eq:bosonopmatel` as blocks labeled by the $\mathrm{U}(1)$ fusion trees

```{figure} ../_static/SymmetricTensors/bosonops.svg
:scale: 12%
:name: bosonops
```

We can then combine these operators in to get the appropriate Hamiltonian terms,

```{figure} ../_static/SymmetricTensors/bosonham.svg
:scale: 12%
:name: bosonham
```

```{note}
Although we have made a suggestive distinction between the 'left' and 'right' versions of the operators $a_L^\pm$ and $a_R^\pm$, one can actually be obtained from the other by permuting the physical and auxiliary indices of the corresponding `TensorMap`s. This permutation has no effect on the actual array blocks of the tensors due to the [bosonic braiding style](https://jutho.github.io/TensorKit.jl/latest/lib/sectors/#TensorKit.BraidingStyle) of $\mathrm{U}(1)$ irreps, so the left and right operators can in essence by seen as the 'same' tensors. This is no longer the case when considering fermionic systems, where permuting indices can in fact change the array blocks as we will see next. As consequence, it is much less clear how to construct two-site symmetric operators in terms of local symmetric objects.
```

The explicit construction then looks something like

```{code-cell} julia
A = U1Space(1 => 1)
```

```{code-cell} julia
:tags: [hide-output]
a⁺ = TensorMap(zeros, ComplexF64, V ← V ⊗ A)
for (s, f) in fusiontrees(a⁺)
    a⁺[s, f] .= sqrt(f.uncoupled[1].charge+1)
end
a⁺
```

```{code-cell} julia
:tags: [hide-output]
a⁻ = TensorMap(zeros, ComplexF64, A ⊗ V ← V)
for (s, f) in fusiontrees(a⁻)
    a⁻[s, f] .= sqrt(f.uncoupled[1].charge)
end
a⁻
```

It is then simple to check that this is indeed what we expect.
```{code-cell} julia
@tensor a⁺a⁻_bis[-1 -2; -3 -4] := a⁺[-1; -3 1] * a⁻[1 -2; -4]
@tensor a⁻a⁺_bis[-1 -2; -3 -4] := a⁻[1 -1; -3] * a⁺[-2; -4 1]
@tensor N_bis[-1 ; -2] := a⁺[-1; 1 2] * a⁻[2 1; -2]

@test norm(a⁺a⁻_bis - a⁺a⁻) ≈ 0 atol=1e-14
@test norm(a⁻a⁺_bis - a⁻a⁺) ≈ 0 atol=1e-14
@test norm(N_bis - N) ≈ 0 atol=1e-14
```

```{note}
From the construction of the Hamiltonian operators [in terms of creation and annihilation operators](bosonham) we clearly see that they are invariant under a transformation $a^\pm \to e^{\pm i\theta} a^\pm$. More generally, any invertible transformation on the auxiliary space leaves the resulting contraction unchanged. This ambiguity in the definition clearly shows that one should really always think in terms of the fully symmetric procucts of $a^+$ and $a^-$ rather than in terms of these operators themselves. In particular, one can always decompose such a symmetric product into the [form above](bosonham) by means of an SVD. 
```


## Level 3: Fermions and the Kitaev Model

While we have already covered quite a lot of ground towards understanding symmetric tensors in terms of fusion trees and corresponding blocks, the symmetries considered so far have been quite 'simple' in the sense that sectors corresponding to irreps of $\mathbb{Z}_2$ and $\mathrm{U}(1)$ have [*Abelian fusion rules*](https://jutho.github.io/TensorKit.jl/latest/lib/sectors/#TensorKit.FusionStyle) and [*bosonic exchange statistics*](https://jutho.github.io/TensorKit.jl/latest/lib/sectors/#TensorKit.BraidingStyle). This means that the fusion of two irreps always gives a single irrep as the fusion product, and that exchanging two irreps in a tensor product is trivial. In practice, this implies that for tensors with these symmetries the fusion trees are completely fixed by the uncoupled charges, which uniquely define both the inner line and the coupled charge, and that tensor indices can be permuted freely without any 'strange' side effects.

In the following we will consider examples with fermionic and even anyonic exchange statistics, and non-Abelian fusion rules. In going through these examples it will become clear that the fusion trees labeling the blocks of a symmetric tensor imply more information than just a labeling.


### Fermion Parity Symmetry

As a simple example we will consider the Kitaev chain, which describes a chain of interacting spinless fermions with nearest-neighbor hopping and pairing terms. The Hamiltonian of this model is given by
```{math}
:label: eq:kitaev
H = \sum_{\langle i,j \rangle} \left(-\frac{t}{2}(c_i^+ c_j^- - c_i^- c_j^+) + \frac{\Delta}{2}(c_i^+ c_j^+ - c_i^- c_j^-) \right) - \mu \sum_{i} N_i
```
where $N_i = c_i^+ c_i^-$ is the local particle number operator. As opposed to the previous case, the fermionic creation and annihilation operators now satisfy the anticommutation relations
```{math}
\begin{align*}
\left\{c_i^-, c_j^-\right\} &= \left\{c_i^+, c_j^+\right\} = 0 \\
\left\{c_i^-, c_j^+\right\} &= \delta_{ij} .\\
\end{align*}
```
These commutation relations justify the choice of the relative minus sign in the hopping and pairing terms. Indeed, since fermionic operators on different sites always anticommute, these relative minus signs are needed to ensure that the Hamiltonian is Hermitian, since $\left( c_i^+ c_j^- \right)^\dagger = c_j^+ c_i^- = - c_i^- c_j^+$ and $\left( c_i^+ c_j^+ \right)^\dagger = c_j^- c_i^- = - c_i^- c_j^-$.
The anticommutation relations also naturally restrict the local occupation number to be 0 or 1, leading to a well-defined notion of *fermion-parity*. The local fermion-parity operator is related to the fermion number operator as $Q_i = (-1)^{n_i}$, and is diagonal in the occupation number basis. The Hamiltonian {eq}`eq:kitaev` is invariant under conjugation by the global fermion-parity operator, $Q H Q^\dagger = H$, where
```{math}
Q = \exp \left( i \pi \sum_i N_i \right) = (-1)^{\sum_i N_i}.
```
This fermion parity symmetry, which we will denote as $f\mathbb{Z}_2$, is a $\mathbb{Z}_2$-like symmetry in the sense that it has a trivial representation, which we call *even* and again denote by '0', and a sign representation which we call *odd* and denote by '1'. The fusion rules of these irreps are the same as for $\mathbb{Z}_2$. Similar to the previous case, the local symmetry operator $Q_i$ is already diagonal, so the occupation number basis coincides with the irrep basis and we don't need an additional basis transform. The important difference with a regular $\mathbb{Z}_2$ symmetry is that the irreps of $f\mathbb{Z}_2$ have fermionic braiding statistics, in the sense that exhanging two odd irreps gives rise to a minus sign.

In TensorKit.jl, an $f\mathbb{Z}_2$-graded vector spaces is represented as a `Vect[FermionParity]` space, where a given $f\mathbb{Z}_2$ irrep can be represented as a [`FermionParity`](https://jutho.github.io/TensorKit.jl/stable/lib/sectors/#TensorKit.FermionParity) sector instance. Using the simplest instance of a vector space containing a single even and odd irrep, we can already demonstrate the corresponding fermionic braiding behavior by [performing a permutation](https://jutho.github.io/TensorKit.jl/stable/lib/tensors/#TensorKit.permute-Union{Tuple{N%E2%82%82},%20Tuple{N%E2%82%81},%20Tuple{S},%20Tuple{AbstractTensorMap{S},%20Tuple{Tuple{Vararg{Int64,%20N%E2%82%81}},%20Tuple{Vararg{Int64,%20N%E2%82%82}}}}}%20where%20{S,%20N%E2%82%81,%20N%E2%82%82}) on a simple `TensorMap`.

```{code-cell} julia
V = Vect[FermionParity](0 => 1, 1 => 1)
t = TensorMap(ones, ComplexF64, V ← V ⊗ V)
```

```{code-cell} julia
permute(t, ((1,), (3, 2)))
```
In other words, when exchanging the two domain vector spaces, the block of the `TensorMap` for which both corresponding irreps are odd picks up a minus sign, exactly as we would expect for fermionic charges.


### Constructing the Hamiltonian

We can directly construct the Hamiltonian terms as symmetric `TensorMap`s using the same procedure as before starting from their matrix elements in the occupation number basis. However, in this case we should be a bit more careful about the precise definition of the basis states in composite systems. Indeed, the tensor product structure of fermionic systems is inherently tricky to deal with, and should ideally be treated in the context of [*super vector spaces*](https://en.wikipedia.org/wiki/Super_vector_space). For two sites, we can define the following basis states on top of the fermionic vacuuum $\ket{00}$:
```{math}
\begin{align*}
\ket{01} &= c_2^+ \ket{00} \\
\ket{10} &= c_1^+ \ket{00} \\
\ket{11} &= c_1^+ c_2^+ \ket{00} \\
\end{align*}
```
This definition in combination with the anticommutation relations above give rise to the nonzero matrix elements
```{math}
\begin{align*}
c_1^+ c_2^- \ket{0, 1} &= \ket{1, 0} \\
c_1^- c_2^+ \ket{1, 0} &= - \ket{0, 1} \\
c_1^+ c_2^+ \ket{0, 0} &= \ket{1, 1} \\
c_1^- c_2^- \ket{1, 1} &= - \ket{0, 0} \\
N \ket{n} &= n \ket{n}.
\end{align*}
```
While the signs in these expressions may seem a little unintuitive at first sight, they are essential to the fermionic nature of the system. Indeed, if we for example work out the matrix element of $c_1^- c_2^+$ we find
```{math}
\begin{align*}
c_1^- c_2^+ \ket{1, 0} = c_1^- c_2^+ c_1^+ \ket{0, 0} = - c_2^+ c_1^- c_1^+ \ket{0, 0} = - c_2^+ (\mathbb{1} - c_1^+ c_1^-) \ket{0, 0} = - c_2^+ \ket{0, 0} = - \ket{0, 1}. \\
\end{align*}
```

Once we have these matrix elements the hard part is done, and we can naively associate these to the following $f\mathbb{Z}_2$ fusion trees with corresponding block values,
```{figure} ../_static/SymmetricTensors/fZ2_fusiontrees.svg
:scale: 12%
:name: fZ2_fusiontrees
```
Given this information, we can go through the same procedure again to construct $c^+ c^-$, $c^- c^+$ and $N$ operators as `TensorMap`s over $f\mathbb{Z}_2$-graded vector spaces.

```{code-cell} julia
V = Vect[FermionParity](0 => 1, 1 => 1)
```

```{code-cell} julia
c⁺c⁻ = TensorMap(zeros, ComplexF64, V ⊗ V ← V ⊗ V)
odd = FermionParity(1)
for (s, f) in fusiontrees(c⁺c⁻)
    if s.uncoupled[1] == odd && f.uncoupled[2] == odd && f.coupled == odd
        c⁺c⁻[s, f] .= 1
    end
end
c⁺c⁻
```

```{code-cell} julia
c⁻c⁺ = TensorMap(zeros, ComplexF64, V ⊗ V ← V ⊗ V)
for (s, f) in fusiontrees(c⁻c⁺)
    if f.uncoupled[1] == odd && s.uncoupled[2] == odd && f.coupled == odd
        c⁻c⁺[s, f] .= -1
    end
end
c⁻c⁺
```

```{code-cell} julia
c⁺c⁺ = TensorMap(zeros, ComplexF64, V ⊗ V ← V ⊗ V)
odd = FermionParity(1)
for (s, f) in fusiontrees(c⁺c⁺)
    if s.uncoupled[1] == odd && f.uncoupled[1] != odd && f.coupled != odd
        c⁺c⁺[s, f] .= 1
    end
end
c⁺c⁺
```

```{code-cell} julia
c⁻c⁻ = TensorMap(zeros, ComplexF64, V ⊗ V ← V ⊗ V)
for (s, f) in fusiontrees(c⁻c⁻)
    if s.uncoupled[1] != odd && f.uncoupled[2] == odd && f.coupled != odd
        c⁻c⁻[s, f] .= -1
    end
end
c⁻c⁻
```

```{code-cell} julia
N = TensorMap(zeros, ComplexF64, V ← V)
for (s, f) in fusiontrees(N)
    N[s, f] .= f.coupled == odd ? 1 : 0
end
N
```

```{note}
Working with fermionic systems is inherently tricky, as can already be seen from something as simple as computing matrix elements of fermionic operators. Similarly, while constructing symmetric tensors that correspond to the symmetric Hamiltonian terms was still quite straightforward, it is far less clear in this case how to construct these terms as contractions of local symmetric tensors representing individual creation and annihilation operators in this case. While such a decomposition can always be in principle obtained using a (now explicitly fermionic) SVD, manually constructing such tensors as we did in the bosonic case is far from trivial. Trying this would be a good exercise in working with fermionic symmetries, but it is not something we will do here.
```


## Level 4: Non-Abelian Symmetries and the Quantum Heisenberg Model

We will now move on to systems which have more complicated *non-Abelian* symmetries. For a non-Abelian symmetry group $G$, the fact that its elements do not all commute has a profound impact on its representation theory. In particular, the irreps of such a group can be higher dimensional, and the fusion of two irreps can give rise to multiple different irreps. On the one hand this means that fusion trees of these irreps are no longer completely determined by the uncoupled charges. Indeed, in this case some of the [internal structure of the `FusionTree` type](fusion_trees) we have ignored before will become relevant (of which we will give an [example below](sun_heisenberg)). On the other hand, it follows that fusion trees of irreps now not only label blocks, but also encode a certain *nontrivial symmetry structure*. We will make this statement more precise in the following, but the fact that this is necessary is quite intuitive. If we recall our original statement that symmetric tensors consist of blocks associated to fusion trees which carry irrep labels, then for higher-dimensional irreps the corresponding fusion trees must encode some additional information that implicitly takes into account the internal structure of the representation spaces. In particular, this means that the conversion of an operator, given its matrix elements in the irrep basis, to the blocks of the corresponding symmetric `TensorMap` is less straightforward since it requires an understanding of exactly what this implied internal structure is. Therefore, we require some more discussion before we can actually move on to an example.

We'll start by discussing the general structure of a `TensorMap` which is symmetric under a non-Abelian group symmetry. We then given an example based on $\mathrm{SU}(2)$, where we construct the Heisenberg Hamiltonian using two different approaches. Finally, we show how the more intuitive approach can be used to obtain an elegant generalization to the $\mathrm{SU}(N)$-symmetric case.


### Block Sparsity Revisited: The Wigner-Eckart Theorem

Let us recall some basics of representation theory first. Consider a group $G$ and a corresponding representation space $V$, such that every element $g \in G$ can be realized as a unitary operator $U_g : V \to V$. Let $h$ be a `TensorMap` whose domain and codomain are given by the tensor product of two of these representation spaces. Recall that, by definition, the statement that '$h$ is symmetric under $G$' means that
```{figure} ../_static/SymmetricTensors/symmetric_tensor.svg
:scale: 12%
:name: symmetric_tensor
```
for every $g \in G$. If we label the irreducible representations of $G$ by $l$, then any representation space can be decomposed into a direct sum of irreducible representations, $V = \bigoplus_l V^{(l)}$, in such a way that $U_g$ is block-diagonal where each block is labeled by a particular irrep $l$. For each irrep space $V^{(l)}$ we can define an orthonormal basis labeled as $\ket{l, m}$, where the auxiliary label $m$ can take $\text{dim}\left( V^{(l)} \right)$ different values. Since we know that tensors are multilinear maps over tensor product spaces, it is natural to consider the tensor product of representation spaces in more detail.

[From the representation theory of groups](https://en.wikipedia.org/wiki/Tensor_product_of_representations#Clebsch%E2%80%93Gordan_theory), it is known that the product of two irreps can in turn be decomposed into a direct sum of irreps, $V^{(l_1)} \otimes V^{(l_2)} \cong \bigoplus_{k} V^{(k)}$. The precise nature of this decomposition, also refered to as the *Clebsch-Gordan problem*, is given by the so-called *Clebsch-Gordan coefficients*, which we will denote as $C^{k}_{l_1,l_2}$. This set of coefficients, which can be interpreted as a $\text{dim}\left( V^{(l_1)} \right) \times \text{dim}\left( V^{(l_2)} \right)  \times \text{dim}\left( V^{(l_3)} \right)$ array, that encodes how a basis state $\ket{k,n} \in V^{(k)}$ corresponding to some term in the direct sum can be decomposed into a linear combination of basis vectors $\ket{l_1,m_1} \otimes \ket{l_2,m_2}$ of the tensor product space:
```{math}
\ket{k,n} = \sum_{m_1, m_2} \left( C^{l_3}_{l_1,l_2} \right)^{n}_{m_1, m_2} \ket{l_1,m_1} \otimes \ket{l_2,m_2}.
```
These recoupling coefficients turn out to be essential to the structure of symmetric tensors, which can be best understood in the context of the [Wigner-Eckart theorem](https://en.wikipedia.org/wiki/Wigner%E2%80%93Eckart_theorem). This theorem implies that for any [`TensorMap` $h$ that is symmetric under $G$](symmetric_tensor), its matrix elements in the tensor product irrep basis are given by the product of Clebsch-Gordan coefficients which characterize the coupling of the basis states in the domain and codomain, and a so-called *reduced matrix element* which only depends on the irrep labels. Concretely, the matrix element $\bra{l_1,m_1} \otimes \bra{l_2,m_2} h \ket{l_3,m_3} \otimes \ket{l_4,m_4}$ is given by
```{figure} ../_static/SymmetricTensors/wignereckart.svg
:scale: 12%
:name: wignereckart
```
Here, the sum runs over all possible irreps $k$ in the fusion product $l_3 \otimes l_4$ and over all basis states $\ket{k,n}$ of $V^{(k)}$. The reduced matrix elements $h_{\text{red}}$ are independent of the basis state labels and only depend on the irrep labels themselves. Each reduced matrix element should be interpreted as being labeled by an irrep fusion tree,
```{figure} ../_static/SymmetricTensors/anotherfusiontree.svg
:scale: 12%
:name: anotherfusiontree
```
The fusion tree itself in turn implies the Clebsch-Gordan coefficients $C^{k}_{l_1,l_2}$ and conjugate coefficients ${C^{\dagger}}_{k}^{l_1,l_2}$ encode the splitting (decomposition) of the coupled basis state $\ket{k,n}$ to the codomain basis states $\ket{l_1,m_1} \otimes \ket{l_2,m_2}$ and the coupling of the domain basis states $\ket{l_3,m_3} \otimes \ket{l_4,m_4}$ to the coupled basis state $\ket{k,n}$ respectively.

The Wigner-Eckart theorem dictates that this structure in terms of Clebsch-Gordan coefficients is necessary to ensure that the corresponding is symmetric. It is precisely this structure that is inherently encoded into the fusion tree part of a symmetric `TensorMap`. In particular, **the array block value associated to each fusion tree in a symmetric tensor is precisely the reduced matrix element in the Clebsch-Gordan decomposition**.

As a small demonstration of this fact, we can make a simple $\mathrm{SU}(2)$-symmetric tensor with trivial block values and verify that its implied symmetry structure exactly corresponds to the expected Clebsch-Gordan coefficient. In TensorKit.jl, a $\mathrm{SU}(2)$-graded vector space is represented as an [`SU2Space`](https://jutho.github.io/TensorKit.jl/latest/lib/spaces/#TensorKit.SU2Space), where a given $\mathrm{SU}(2)$ irrep can be represented as an [`SU2Irrep`](https://jutho.github.io/TensorKit.jl/latest/lib/sectors/#TensorKit.SU2Irrep) instance of integer or halfinteger spin as encoded in its `j` field. If we construct a `TensorMap` whose symmetry structure corresponds to the coupling of two spin-$\frac{1}{2}$ irreps to a spin-$1$ irrep, we can then convert it to a plain array and compare it to the $\mathrm{SU}(2)$ Clebsch-Gordan coefficients exported by the [WignerSymbols.jl package](https://github.com/Jutho/WignerSymbols.jl).

```{code-cell} julia
V1 = SU2Space(1 => 1)
V2 = SU2Space(1//2 => 1)
t = TensorMap(ones, ComplexF64, V1 ← V2 ⊗ V2)
```

```{code-cell} julia
ta = convert(Array, t)
```
The conversion gives us a $3 \times 3 \times 2$ array, which exactly corresponds to the size of the $C_{1}^{\frac{1}{2},\frac{1}{2}}$ Clebsch-Gordan array. In order to explicitly compare whether the entries match we need to know the ordering of basis states assumed by TensorKit.jl when converting the tensor to its matrix elements in the irrep basis. For $\mathrm{SU}(2)$ the irrep basis is ordered in ascending magnetic quantum number $m$, which gives us a map $m = i - (l+1)$ for mapping an array index to a corresponding magnetic quantum number for the spin-$l$ irrep.
```{code-cell} julia
for i1 in 1:dim(V1), i2 in 1:dim(V2), i3 in 1:dim(V2)
    # map basis state index to magnetic quantum number
    m1 = i1 - (1 + 1)
    m2 = i2 - (1//2 + 1)
    m3 = i3 - (1//2 + 1)
    @test ta[i1, i2, i3] ≈ clebschgordan(1//2, m2, 1//2, m3, 1, m1)
end
```

Based on this discussion, we can quantify the aforementioned 'difficulties' in the inverse operation of what we just demonstrated, namely converting a given operator to a symmetric `TensorMap` given only its matrix elements in the irrep basis. Indeed, it is now clear that this precisely requires isolating the reduced matrix elements introduced above. Given the matrix elements of the operator in the irrep basis, this can in general be done by solving the system of equations implied by the [Clebsch-Gordan decomposition](wignereckart). A simpler way to achieve the same thing is to make use of the fact that the [Clebsch-Gordan tensors form a complete orthonormal basis](https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients#Orthogonality_relations) on the coupled space. Indeed, by projecting out the appropriate Clebsch-Gordan coefficients and using their orthogonality relations, we can construct a diagonal operator on each coupled irrep space $V^{(k)}$. Each of these diagonal operators is proportional to the identity, where the proportionality factor is precisely the reduced matrix element associated to the corresponding irrep fusion tree.
```{figure} ../_static/SymmetricTensors/none2symm.svg
:scale: 12%
:name: none2symm
```

This procedure works for any group symmetry, and all we need are matrix elements of the operator in the irrep basis and the Clebsch-Gordan coefficients. In the following we demonstrate this explicit procedure for the particular example of $G = \mathrm{SU}(2)$. However, it should be noted that for general groups the Clebsch-Gordan coefficients may not be as easy to compute (if it is even known how). In addition, the procedure for manually projecting out the reduced matrix elements requires being particularly careful about the correspondence between the basis states used to define the original matrix elements and those implied by the Clebsch-Gordan coefficients. Therefore, it is often easier to directly construct the symmetric tensor based on some representation theory, as we will see below.


### The 'Generic' Approach to the Spin-1 Heisenberg Model: Wigner-Eckart in Action

Consider the spin-1 Heisenberg model with Hamiltonian
```{math}
H = J \sum_{\langle i,j \rangle} \vec{S}_i \cdot \vec{S}_j
```
where $\vec{S} = (S^x, S^y, S^z)$ are the spin operators. The physical Hilbert space at each site is the three-dimensional spin-1 irrep of $\mathrm{SU}(2)$. Each two-site exchange operator $\vec{S}_i \cdot \vec{S}_j$ in the sum commutes with a global transformation $g \in \mathrm{SU}(2)$, so that it satisfies the [above symmetry condition](symmetric_tensor). Therefore, we can represent it as an $\mathrm{SU}(2)$-symmetric `TensorMap`, as long as we can isolate its reduced matrix elements.

In order to apply the above procedure, we first require the matrix elements in the irrep basis. These can be constructed as a $3 \times 3 \times 3 \times 3$ array `SS` using the [familiar representation of the $\mathrm{SU}(2)$ generators in the spin-1 representation](https://en.wikipedia.org/wiki/Spin_(physics)#Higher_spins), with respect to the $\{\ket{1,-1}, \ket{1,0}, \ket{1,1}\}$ basis.
```{code-cell} julia
:tags: [hide-output]
Sx = 1 / sqrt(2) * ComplexF64[0 1 0; 1 0 1; 0 1 0]
Sy = 1 / sqrt(2) * ComplexF64[0 1im 0; -1im 0 1im; 0 -1im 0]
Sz = ComplexF64[-1 0 0; 0 0 0; 0 0 1]

@tensor SS_arr[-1 -2; -3 -4] := Sx[-1; -3] * Sx[-2; -4] + Sy[-1; -3] * Sy[-2; -4] + Sz[-1; -3] * Sz[-2; -4]
```

The next step is to project out the reduced matrix elements by taking the overlap with the appropriate Clebsch-Gordan coefficients, according to [the Clebsch-Gordan decomposition given above](none2symm). In our current case of a spin-1 physical space we have $l_1 = l_2 = l_3 = l_4 = 1$, and the coupled irrep $k$ can therefore take the values $0, 1, 2$. The reduced matrix element for a given $k$ can then be implemented in the following way:
```{code-cell} julia
:tags: [hide-output]
function get_reduced_element(k)
    # construct Clebsch-Gordan coefficients for coupling 1 ⊗ 1 to k   
    CG = zeros(ComplexF64, 3, 3, 2*k + 1)
    for m1 in -k:k, m2 in -1:1, m3 in -1:1
        CG[m2 + 2, m3 + 2, m1 + k + 1] = clebschgordan(1, m2, 1, m3, k, m1)
    end

    # project out diagonal matrix on coupled irrep space
    @tensor reduced_matrix[-1; -2] := CG[1 2; -1] * SS_arr[1 2; 3 4] * conj(CG[3 4; -2])

    # check that it is proportianal to the identity
    @assert isapprox(reduced_matrix, reduced_matrix[1, 1] * I; atol=1e-12)

    # return the proportionality factor
    return reduced_matrix[1, 1]
end
```
If we use this to compute the reduced matrix elements for $k = 0, 1, 2$,
```{code-cell} julia
get_reduced_element(0)
```
```{code-cell} julia
get_reduced_element(1)
```
```{code-cell} julia
get_reduced_element(2)
```
we can read off the entries
```{math}
\renewcommand\thickspace{\kern .01ex}
\left[ (\vec{S}_i \cdot \vec{S}_j)_\text{red} \right] \,\,\!\!
\begin{smallmatrix}
    1,1\\
    0\\
    1,1
\end{smallmatrix} = -2, \quad
\left[ (\vec{S}_i \cdot \vec{S}_j)_\text{red} \right] \,\,\!\!
\begin{smallmatrix}
    1,1\\
    1\\
    1,1
\end{smallmatrix} = -1, \quad
\left[ (\vec{S}_i \cdot \vec{S}_j)_\text{red} \right] \,\,\!\!
\begin{smallmatrix}
    1,1\\
    2\\
    1,1
\end{smallmatrix} = 1, \quad
```
These can then be used to construct the symmetric `TensorMap` representing the exchange interaction:
```{code-cell} julia
V = SU2Space(1 => 1)
SS = TensorMap(zeros, ComplexF64, V ⊗ V ← V ⊗ V)
for (s, f) in fusiontrees(SS)
    k = Int(f.coupled.j)
    SS[s, f] .= get_reduced_element(k)
end
SS
```


### An 'Elegant' Approach to the Heisenberg Model

As noted above, the explicit procedure of projecting out the reduced matrix elements from the action of an operator in the irrep basis can be a bit cumbersome for more complicated groups. However, using some basic representation theory we can bypass this step altogether for the Heisenberg model. First, we rewrite the exchange interaction in the following way:
```{math}
:label: eq:casimir_decomp
\vec{S}_i \cdot \vec{S}_j = \frac{1}{2} \left( \left( \vec{S}_i + \vec{S}_j \right)^2 - \vec{S}_i^2 - \vec{S}_j^2 \right)
```
Here, $\vec{S}_i$ and $\vec{S}_j$ are spin operators on the physcial irrep, while total spin operator $\vec{S}_i + \vec{S}_j$ can be decomposed onto the different coupled irreps $k$. It is a well known fact that the quadratic sum of the generators of $\mathrm{SU}(2)$, often refered to as the [*quadratic Casimir*](https://en.wikipedia.org/wiki/Representation_theory_of_SU(2)#The_Casimir_element), commutes with all generators. By [Schur's lemma](https://en.wikipedia.org/wiki/Schur%27s_lemma), it must then act proportional to the identity on every irrep, where the corresponding eigenvalue is determined by the spin irrep label. In particular, we have for each irrep $l$
```{math}
\vec{S}^2 \ket{l,m} = l(l+1) \ket{l,m}.
```
It then follows from Eq. {eq}`eq:casimir_decomp` that the reduced matrix elements of the exchange interaction are completely determined by the eigenvalue of the quadratic Casimir on the uncoupled and coupled irreps. Indeed, to each fusion tree we can associate a well-defined value
```{figure} ../_static/SymmetricTensors/SU2_fusiontrees.svg
:scale: 12%
:name: SU2_fusiontrees
```
This gives us all we need to directly construct the exchange interaction as a symmetric `TensorMap`,
```{code-cell} julia
V = SU2Space(1 => 1)
SS = TensorMap(zeros, ComplexF64, V ⊗ V ← V ⊗ V)
for (s, f) in fusiontrees(SS)
    l3 = f.uncoupled[1].j
    l4 = f.uncoupled[2].j
    k = f.coupled.j
    SS[s, f] .= (k * (k + 1) - l3 * (l3 + 1) - l4 * (l4 + 1)) / 2
end
SS
```
which gives exactly the same result as the previous approach.

```{note}
This last construction for the exchange interaction immediatly generalizes to any value of the physical spin. All we need is to fill in the appropriate values for the uncoupled irreps $l_1$, $l_2$, $l_3$ and $l_4$.
```


(sun_heisenberg)=
### $\mathrm{SU}(N)$ generalization

We end this subsection with some comments on the generalization of the above discussion to $\mathrm{SU}(N)$. As foreshadowed above, the irreps of $\mathrm{SU}(N)$ in general have an even more complicated structure. In particular, they can admit so-called *fusion multiplicities*, where the fusion of two irreps can not just have multiple distinct outcomes, but they can even fuse to a given irrep in mutliple inequivalent ways. We can demonstrate this behavior for the adjoint representation of $\mathrm{SU}(3)$. For this we can use the the [SUNRepresentations.jl](https://github.com/maartenvd/SUNRepresentations.jl/tree/master) package which provides an interface for working with irreps of $\mathrm{SU}(N)$ and their Clebsch-Gordan coefficients. A particular representation is represented by an `SUNIrrep{N}` which can be used with TensorKit.jl. The eight-dimensional adjoint representation of $\mathrm{SU}(3)$ is given by
```{code-cell} julia
:tags: [remove-cell]
SUNRepresentations.display_mode("dimension")
```
```{code-cell} julia
l = SU3Irrep("8")
```
If we look at the possible outcomes of fusing two adjoint irreps, we find the by now familiar non-Abelian fusion behavior,
```{code-cell} julia
collect(l ⊗ l)
```
However, this particular fusion has multiplicities, since the adjoint irrep can actually fuse to itself in two distinct ways. The full decomposition of this fusion product is given by
```{math}
\mathbf{8} \otimes \mathbf{8} = \mathbf{1} \oplus \mathbf{3} \oplus 2 \cdot \mathbf{8} \oplus \mathbf{10} \oplus \mathbf{\overline{10}} \oplus \mathbf{27}
```
This fusion multiplicity can be detected by using [`Nsymbol`](https://jutho.github.io/TensorKit.jl/latest/lib/sectors/#TensorKit.Nsymbol) method from TensorKit.jl to inspect the number of times `l` appears in the fusion product `l ⊗ l`, 
```{code-cell} julia
Nsymbol(l, l, l)
```
When working with irreps with fusion multiplicities, each `FusionTree` carries additional `vertices` labels which label which of the distinct fusion vertices is being referred to. We will return to this at the end of this section.

Given the generators $T^k$ of $\mathrm{SU}(N)$, we can define a generalized Heisenberg model using a similar exchange interaction, giving the Hamiltonian
```{math}
H = J \sum_{\langle i,j \rangle} \vec{T}_i \cdot \vec{T}_j
```
For a particular choice of physical irrep, the exchange interaction can again be constructed as a symmetric `TensorMap` by first rewriting it as 
```{math}
\vec{T}_i \cdot \vec{T}_j = \frac{1}{2} \left( \left( \vec{T}_i + \vec{T}_j \right)^2 - \vec{T}_i^2 - \vec{T}_j^2 \right).
```
For any $N$, the quadratic Casimir
```{math}
C_1 = \sum_k T^k T^k
```
commutes with all $\mathrm{SU}(N)$ generators, meaning it has a well defined eigenvalue in each irrep. This observation then immediately given the reduced matrix elements of the exchange interaction as
```{figure} ../_static/SymmetricTensors/SUN_fusiontrees.svg
:scale: 12%
:name: SUN_fusiontrees
```
Using these to directly construct the corresponding symmetric `TensorMap` is much simpler than going through the explicit projection procedure using Clebsch-Gordan coefficients.

For the particular example of $\mathrm{SU}(3)$, the generators are given by $T^k = \frac{1}{2} \lambda^k$ , where $\lambda^k$ are the [Gell-Mann matrices](https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients_for_SU(3)#Generators_of_the_Lie_algebra). Each irrep can be labeled as $l = D(p,q)$ where $p$ and $q$ are refered to as the *Dynkin labels*. The eigenvalue of the quadratic Casimir for a given irrep is given by [Freudenthal's formula](https://en.wikipedia.org/wiki/Weyl_character_formula#Freudenthal's_formula),
```{math}
C_1(D(p,q)) = \frac{1}{3} (p^2 + q^2 + 3p + 3q + pq).
```
Using SUNRepresentations.jl, we can compute the Casimir as
```{code-cell} julia
:tags: [hide-output]
function casimir(l::SUNIrrep)
    p, q = dynkin_label(l)
    return (p^2 + q^2 + 3 * p + 3 * q + p * q) / 3
end
```
If we use the adjoint representation of $\mathrm{SU}(3)$ as physical space, the Heisenberg exchange interaction can then be constructed as
```{code-cell} julia
V = Vect[SUNIrrep{3}](SU3Irrep("8") => 1)
TT = TensorMap(zeros, ComplexF64, V ⊗ V ← V ⊗ V)
for (s, f) in fusiontrees(TT)
    l3 = f.uncoupled[1]
    l4 = f.uncoupled[2]
    k = f.coupled
    TT[s, f] .= (casimir(k) - casimir(l3) - casimir(l4)) / 2
end
TT
```
Circling back to our earlier remark, we clearly see that the fusion trees of this tensor indeed have non-trivial vertex labels.
```{code-cell} julia
f = collect(fusiontrees(TT))[3][2]
```
```{code-cell} julia
f.vertices
```

```{note}
While we have given an explicit example using $\mathrm{SU}(3)$ with the adoint irrep on the physical level, the same construction holds for the general $\mathrm{SU}(N)$ with arbitrary physical irreps. All we require is the expression for the eigenvalues of the quadratic Casimir in each irrep.
```


## Level 5: Anyonic Symmetries and the Golden Chain

While we have focussed exclusively on group-like symmetries in our discussion so far, the framework of symmetric actually extends beyond groups to so-called [*categorical symmetries*](https://jutho.github.io/TensorKit.jl/stable/man/sectors/#ss_representationtheory). These are quite exotic symmetries characterized in terms of [the topological data of a unitary fusion category](https://jutho.github.io/TensorKit.jl/stable/man/categories/#ss_topologicalfusion). While the precise details of all the terms in these statements fall beyond the scope of this tutorial, we can give a simple example of a Hamiltonion model with a categorical symmetry called [the golden chain](https://arxiv.org/abs/cond-mat/0612341).

This is a one-dimensional system defined as a spin chain, where each physical 'spin' corresponds to a so-called [Fibonacci anyon](https://arxiv.org/abs/0902.3275). There are two such Fibonacci anyons, which we will denote as $1$ and $\tau$. They obey the fusion rules
```{math}
1 \otimes 1 = 1, \quad 1 \otimes \tau = \tau, \quad \tau \otimes \tau = 1 \oplus \tau.
```
The Hilbert space of a chain of Fibonacci anyons is not a regular tensor product space, but rather a *constrained Hilbert space* where the only allowed basis states are labeled by valid Fibonacci fusion configurations.
In the golden chain model, we define a nearest-neighbor Hamiltonian on this Hilbert space by imposing an energy penalty when any two neighboring anyons fuse to a $\tau$ anyon.

Even just writing down an explicit expression for this interaction on such a constrained Hilbert space is not entirely straightforward. However, using the framework of symmetric tensors it can actually be explicitly constructed in a very straightforward way. Indeed, TensorKit.jl supports a dedicated [`FibonacciAnyon`](https://jutho.github.io/TensorKit.jl/latest/lib/sectors/#TensorKit.FibonacciAnyon) sector type which can be used to construct precisely such a constrained Fibonacci-graded vector space. A Hamiltonian
```{math}
H = \sum_{\langle i,j \rangle} h_{ij}
```
which favors neighboring anyons fusing to the vacuum can be constructed as a `TensorMap` on the product space of two Fibonacci-graded physical spaces
```{code-cell} julia
:tags: [hide-output]
V = Vect[FibonacciAnyon](:τ => 1)
```
and assigning the following nonzero block value to the two-site fusion trees
```{figure} ../_static/SymmetricTensors/Fib_fusiontrees.svg
:scale: 12%
:name: Fib_fusiontrees
```
This allows us to define this, at first sight, exotic and complicated Hamiltonian in a few simple lines of code,
```{code-cell} julia
h = TensorMap(ones, V ⊗ V ← V ⊗ V)
for (s, f) in fusiontrees(h)
    h[s, f] .= f.coupled == FibonacciAnyon(:I) ? -1 : 0
end
h
```

```{note}
In the previous section we have stressed the role of Clebsch-Gordan coefficients in the structure of symmetric tensors, and how they can be used to map between the representation of an operator in the irrep basis and its symmetric tensor representation. However, for categorical symmetries such as the Fibonacci anyons, there are no Clebsch-Gordan coefficients. Therefore, the 'matrix elements of the operator in the irrep basis' are not well-defined, meaning that a Fibonacci-symmetric tensor cannot actually be converted to a plain array in a meaningful way.
```


## To Be Added

- Add section on product symmetries and how to work with them.
    - Use the Hubbard model with $f\mathbb{Z}_2 \boxtimes \mathrm{U}(1) \boxtimes \mathrm{SU}(2)$ as an example.
- Add a section on classical $O(N)$ models to illustrate the ('Fourier') transformation from the group element to the irrep basis for continuous symmetries.

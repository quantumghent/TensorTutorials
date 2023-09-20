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

(tensor_network_states)=
# Tensor Network States

After our introduction on [quantum many body systems](many_body_intro) and [tensor networks](tensor_networks), we move on to considering how tensor networks can characterize many-body systems. We start with a constructive approach to approximating an arbitrary quantum state by a *tensor network state*. We then qualify in what settings such a representation is efficient, and introduce several classes of tensor network states used in different settings. We end this section by broadly commenting on how efficient manipulations of tensor network states can be used to simulate quantum systems

## Quantum States as Tensor Networks

Consider a quantum many body system which consists of physical spins with a local Hilbert space $
\mathbb H_i = \mathbb C^d $ of dimension $ d $, which we will call the *physical dimension*, are located at every site $ i $ of some lattice $ \Lambda $. This gives rise to a total Hilbert space of the system $ \mathbb H
= \bigotimes_{i = 1}^{N} \mathcal H_\lambda = \left( \mathbb C^d \right)^{\otimes
N}$ where $ N = |\Lambda| $ is the total number of sites in the lattice. A general quantum state in this [many-body Hilbert space](many_body) can be represented in terms a set of $d^N$ complex coefficients $C_{s_1,s_2,...,s_N} \in \mathbb C$, where $s_i\in \{0,...,d-1\}$, with respect to the computational basis as
```{math}
\ket{\psi} = \sum_{s_1,s_2,...,s_N} C_{s_1,s_2,...,s_N}\ket{s_1,s_2,...,s_N}.
```
The exponential increase in the number of coefficients with the system size means that it is entirely impossible to store the full state vector of a quantum system of any reasonable size in this way. For example, a system of $N=100$ spins with $d=2$ has $2^{100} \approx 10^{30}$ coefficients, which is far more than the number of atoms in the universe. 

Instead of directly storing this full state vector, we can alternatively parametrize it as a tensor network. Consider for example the case $N=4$. We can then represent the state vector as a tensor $C_{s_1,s_2,s_3,s_4}$ with four indices, where each index corresponds to a physical spin. The full state is then recovered as
```{figure} ../_static/TensorNetworkStates/full_state.svg
:scale: 12%
:name: full_state
:align: center
```
We can now split the full tensor $C$ into separate components by consecutively applying the SVD between pairs of physical indices. For example, splitting out the first index we can rewrite $C$ as
 ```{figure} ../_static/TensorNetworkStates/svd1.svg
:scale: 12%
:name: svd1
:align: center
```
In this expression we can interpret $L^{(1)}$ a a $d \times D$ matrix, $\lambda^{(1)}$ as a $D \times D$ matrix and $R^{(1)}$ as a $D$ by $d^{N-1}$ matrix. The horizontal edge in this diagram is called a *virtual bond* and the dimension $D$ of this bond is called the *bond dimension*. The bond dimension is a measure of the entanglement in the state, and in this case encodes the amount of entanglement between site 1 and the rest of the system. So far we have not actually done anything significant, since this decomposition of $C$ in general contains exactly the same amount of parameters. The key point is that we can reduce the number of parameters by *truncating* $\lambda^{(1)}$ to only keep the $D$ largest singular values. This results in a *low rank approximation* of the original state, where the quality of the approximation is controlled by the chosen final bond dimension $D$.

By repeatedly apply this procedure, grouping and splitting indices in the resulting diagrams and absorbing the bond tensors $\lambda^{(i)}$ into the site tensors we can decompose $C$ into a tensor network of any geometry. For example, we can approximate $C$ as the contraction of a square network to end up with a *tensor network state* of the form
```{figure} ../_static/TensorNetworkStates/tn_state.svg
:scale: 12%
:name: full_state
:align: center
```
In words, this expression means that for every basis state $ \ket{s_1,s_2,s_3,s_4} $ its corresponding coefficient in the superposition is obtained by indexing all of the *physical legs* pointing downward according to the corresponding physical basis state and contracting the resulting network.

We can therefore parametrize an arbitrary quantum state in terms of a set of local tensors $A^{(i)}$, where each of these tensors encodes a number of parameters that is polynomial in its physical dimension $d$ and bond dimensions $D$ (which can in principle be differnt for every virtual bond). For a general quantum state however, a good tensor network state approximation requires a bond dimension which scales exponentially with the system size, meaning that we have not actually gained anything in terms of efficiency. However, it turns out that for many physically relevant states the bond dimension can be bounded by a constant independent of the system size, in which case the tensor network representation then given an exponential reduction in the number of variational parameters.

(area_laws)=
## Area Laws and Tensor Network States

To see why this is the case, let us study the entanglement entropy of a tensor network state. Consider the following two-dimensional network, where all physical indices have a dimension $d$ and we assume all virtual bonds have the same dimension $D$,
```{figure} ../_static/TensorNetworkStates/peps.svg
:scale: 12%
:name: peps
:align: center
```
We now want to quantify the entanglement between the shaded region $ \mathcal A$ with the rest of the system for this specific state. To this end, we first recall the formula for the bipartite entanglement entropy Eq. {eq}`entanglement_entropy`, and note that the number of terms in this expression is determined by the number of nonzero Schmidt coefficients, the latter of which is referred to as the *Schmidt rank*. Looking back now at our initial decomposition of the full state tensor $C$ by splitting out its first index above, we see that the Shchmidt rank is precisely given by the bond dimension $D$ across this cut. From this, you should be able to convince yourself that the entanglement entropy across this cut is determined by the bond dimension as $S \sim \log(D)$. Extending this line of reasoning to our question of the entanglement between the region $ \mathcal A$ and the rest of the system, we see that each virtual leg connecting $\mathcal A$ to the rest of the system contributes a term $\log(D)$ to the entanglement entropy. Therefore we arrive at
```{math}
S(\mathcal A) \sim \log(D) \partial \; \mathcal A,
```
where $ \partial \mathcal A $ is the size boundary of $\mathcal A$ (which in this two-dimensional case is its circumference).

Clearly, this tensor network state then naturally obeys an area law for its entanglement entropy. In our discussion of the [low temparature properties of quantum many body systems](zero_temp) however, we have already seen that low-energy states of locally interacting Hamiltonians obey exactly such an area law. It is this fact that tensor network states inherently encode area law entanglement that makes them so well suited for representing low-energy states of quantum systems. They can only target a tiny corner of the full exponentially large Hilbert space, but this corner is precisely where the most relevant physics happens. This observation has given rise to a large family of tensor network states which allow for an efficient paramtrization of states with varying geometries.

<!-- TODO: figure with MPS, PEPS, tree, MERA, ... -->

```{note}
An equally important feature of tensor networks is that they, aside from providing an efficient parametrization of states, they also allow for efficient *manipulations* of these states. This means that they can be used to compute interesting features of quantum systems, and can be optimized to target states of specific interest such as ground states and low-lying excitations. For all of the network geometries depicted above there exist corresponding algorithms that put them to efficient use, some of which will be highlighted in future sections of this tutorial.
```

<!-- TODO: any more details on algorithms here? -->

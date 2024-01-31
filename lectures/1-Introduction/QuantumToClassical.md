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

(quantum_to_classical)=
# Quantum-to-Classical Mapping

In this final section, we introduce a general technique that essentially enables us to map
any quantum lattice system in $d$ dimensions to a classical partition function in $d+1$
dimensions, up to some caveats that we will return to at the end of this section.

## Suzuki-Trotter decomposition

Remember that thermal expectation values are given by

$$\braket{\hat{O}} = \mathrm{Tr}\left[\hat{O} \mathrm{e}^{-\beta \hat{H}}\right]/Z(\beta)= \mathrm{Tr}\left[\mathrm{e}^{-\beta \hat{H}/2} \hat{O} \mathrm{e}^{-\beta \hat{H}/2}\right]/Z(\beta)$$

with the thermal partition function $Z(\beta)$ given by

$$Z(\beta) = \mathrm{Tr} \mathrm{e}^{-\beta \hat{H}}.$$

The ground state physics is encoded in the limit $\beta \to \infty$. Note that, if the
system has a unique ground state, we can obtain the ground state $\ket{\Psi_0}$ of a quantum
system by starting from essentially a random state $\ket{\Phi}$ and evolving it in imaginary
time $\tau = -\mathrm{i} t$ for sufficiently long

$$\ket{\Psi_0} \sim \lim_{\tau \to \infty} \mathrm{e}^{-\tau \hat{H}} \ket{\Phi}$$

Expanding the initial state $\ket{\phi}$ in the energy eigenbasis of $\hat{H}$, we see that
the only condition is that it is not orthogonal to the ground state (subspace). In addition,
the ground state will be well approximated if $\tau \Delta E \gg 1$, with 
$\Delta E=E_1 - E_0$ the energy gap. This imaginary time evolution also forms the basic
ingredient of several numerical algorithms for approximating ground states of quantum many
body systems, often in combination with the Suzuki-Trotter decomposition which is introduced
below.

Using this approach, the following expression for the ground state expectation value of on
operator $\hat{O}$ is obtained

$$\braket{\hat{O}} = \lim_{\tau\to\infty} \braket{\phi\vert \mathrm{e}^{-\tau \hat{H}} \hat{O} \mathrm{e}^{-\tau \hat{H}}\vert \phi}/\braket{\phi\vert\mathrm{e}^{-2\tau \hat{H}} \vert\phi}$$

This expression can be compared to the thermal expectation value with $\beta = 2\tau$; the
only difference is in the boundary conditions.

For a quantum many body system, taking the exponential is as hard as determining the full
diagonalisation of the hamiltonian, which is impossible due to the exponentially large
Hilbert space. If the hamiltonian is a sum of local terms, each of these terms can be
exponentiated easily, but for arbitrary $\tau$ there is no relation ship between
$\exp(-\tau \sum_{i} \hat{h}_i)$ and the individual $\exp(-\tau \hat{h}_i)$, unless the
different $\hat{h}_i$ commute. However, for an infinitesimal time step $\epsilon$, we can
use to Baker-Campbell-Hausdorff formula (or better yet, the Zassenhaus formula) to obtain

$$\exp(-\epsilon \sum_{i} \hat{h}_i) = \prod_i \exp(-\epsilon \hat{h}_i) +
\mathcal{O}(\epsilon^2).$$

This then leads to the **Suzuki-Trotter decomposition**
 
$$\exp\left(-\tau \sum_{i} \hat{h}_i\right) = \lim_{M\to\infty} \left( \mathrm{e}^{-\frac{\tau}{M} \sum_i \hat{h}_i} \right)^M =\lim_{M\to\infty} \left(\prod_i \mathrm{e}^{-\frac{\tau}{M} \hat{h}_i} + \mathcal{O}(\tau^2/M^2)\right)^M = \lim_{M\to\infty} \left( \left[\prod_i \mathrm{e}^{-\frac{\tau}{M} \hat{h}_i}\right]^M + \mathcal{O}(\tau^2/M)\right)$$

The product in the final expressions requires chosing a specific order, exactly because the terms $\hat{h}_i$
and thus also the factors $\mathrm{e}^{-\frac{\tau}{M} \hat{h}_i}$ do not commute. The approximation
and error term are valid for arbitrary choices of ordering, but different orderings are not
equivalent. Particular choices can be more suitable for particular purposes. Furthermore,
note that splitting the time interval $[0,\tau]$ into small segments $\epsilon = \tau/M$ is
also the starting point for deriving a path integral representation of the quantum partition
function. The next step is to insert a resolution of the identity in between the $N$
different factors, where the labels of the basis will behave as classical degrees of
freedom. For obtaining a path integral, the basis should be labeled by a number of
continuous degrees of freedom, which can then become continuous functions of time in the
limit $\epsilon\to 0$. Here, instead, we will keep $\epsilon$ small but finite, and use a
discrete basis. 

## From quantum to statistical mechanics

Let's start with a quantum system in $d=0$, i.e. a small number of spins, or in particular, a single spin,
described by a hamiltonian

$$\hat{H} = - h_x \sigma^x - h_z \sigma^z$$

While we could in principle exponentiate $\hat{H}$ directly as it is a $2 \times 2$ matrix,
we will treat it using the Suzuki-Trotter decomposition. Throughout the remainder of this
section, we will use the $\sigma^z$ basis, which we denote as $\ket{1} = \ket{\uparrow}$ and
$\ket{-1} = \ket{\downarrow}$. Inserting resolutions of the identity, we write

$$Z(\beta) = \mathrm{Tr} \mathrm{e}^{-\beta \hat{H}} = \sum_{\{s_k\}=\pm 1}
\braket{s_1|\mathrm{e}^{-\epsilon \hat{H}}|s_2}\cdots \braket{s_M-1| \mathrm{e}^{-\epsilon
\hat{H}}|s_{M}} \cdots \braket{s_M| \mathrm{e}^{-\epsilon \hat{H}}|s_{M+1}}$$

with $s_{M+1} = s_1$, $M\epsilon =\beta$, and where

```{math}
\braket{s_{i}|\mathrm{e}^{-\epsilon \hat{H}}|s_{i+1}} &= \braket{s_{i} \vert \mathrm{e}^{-\epsilon H} \vert s_{i+1}} \approx \braket{s_{i} \vert \mathrm{e}^{\epsilon h_z \sigma^z}\mathrm{e}^{\epsilon h_x \sigma^x} \vert s_{i+1}}\\ 
&= \mathrm{e}^{\epsilon h_z s_i} \braket{s_{i} \vert \cosh(\epsilon h_x) \mathbb{1} + \sinh(\epsilon h_x) \sigma^x \vert s_{i+1}}\\
 &= \mathrm{e}^{K s_{i}s_{i+1} + h s_i + f_0}
```

where the parameters in the last line are given by
$K= -\frac{1}{2}\log \tanh(\epsilon h_x)$, $h = \epsilon h_z$ and
$f_0 =\frac{1}{2}\log[\cosh(\epsilon h_x)\sinh(\epsilon h_x)]$. 

We thus obtain

$$Z(\beta) = \sum_{s_k} \mathrm{e}^{\sum_{i=1}^{M}K s_i s_{i+1} + h s_i},$$

the partition function of the one-dimensional classical Ising model with periodic boundary
conditions. Indeed, $\braket{s_{i}|\mathrm{e}^{-\epsilon H}|s_{i+1}}$ does exactly
correspond to the transfer matrix, and diagonalising the transfer matrix is the most
straightforward approach to solving the one-dimensional classical Ising model.

## Higher dimensional generalisation

We now apply the same approach to the Ising model with both transverse and longitudinal
field in $d$ dimensions, on a hypercubic lattice. We separate the hamiltonian in two parts
according to

$$\hat{H} = \left(-J\sum_{\braket{i,j}} \sigma^z_i \sigma^z_j  - h_z \sum_{i} \sigma^z_i\right) + \left(- h_x \sum_{i} \sigma^x_i\right)= \hat{H}_1 + \hat{H}_2$$

Note that $\hat{H}_1$ and $\hat{H}_2$ in itself contain commuting terms, but of course don't
mutually commute. We follow the same strategy, and will in every (imaginary) time step
introduce a resolution of the identity using the tensor product $\sigma^z$ basis. We now denote
the basis at time step $k$ as $\ket{\{s_{i,k}\}}$, where $i$ labels a site in the $d$
dimensional lattice hosting the quantum degrees of freedom, and $k$ labels points along the
imaginary time axis, which emerges as a new dimension in the problem. We find

$$\exp(-\epsilon \hat{H}_1)\ket{\{s_{i,k}\}} = \exp(\epsilon J \sum_{\braket{i,j}} s_{i,k} s_{j,k}+\epsilon h \sum_i s_{i,k}) \ket{\{s_{i,k}\}}$$

as $\hat{H}_1$ is diagonal in this basis, and 

$$\braket{\{s_{i,k}\}|\exp(-\epsilon \hat{H}_2)|\{s_{i,k+1}\}} = \prod_i \braket{s_{i,k}|\mathrm{e}^{-\epsilon h_x \sigma^x_i} | s_{i,k+1}} \sim \exp(K_\perp\sum_{i} s_{i,k} s_{i,k+1})$$
 
with $K_\perp = \log \tanh(\epsilon h_x)$ as before. Here, we have now ignored an overall
proportionality factor, which is irrelevant when using the partition function to compute
expectation values. With this, we find

$$ Z(\beta) = \sum_{\{s_i,k\}} \exp\left(\sum_{k=1}^{M}\sum_{i} K_\perp s_{i,k} s_{i,k+1}+\sum_{k=1}^M\sum_{\braket{i,j}} K_\parallel s_{i,k} s_{j,k} + \sum_{k=1}^{M}\sum_{i} h s_{i,k} \right)$$

with $K_{\parallel} = \epsilon J$ and $h = \epsilon h_z$. We thus find the partition
function of the classical Ising model in $d+1$ dimensions with anisotropic interaction
strengths, periodic boundary condition in the imaginary time direction and a number of sites
in the time direction given by $M = \beta/\epsilon$. Hence, the ground state regime
$\beta\to \infty$ corresponds to the thermodynamic limit in this additional time direction
of the corresponding classical system, so that there are many similarities (or actually,
equivalences) between between quantum phenomena in $d$ dimensions and classical phenomena in
$D=d+1$ dimensions. On the other hand, when the quantum system is at finite temperature
$\beta$, the additional dimension is finite and never in the thermodynamic limit. In that
case, this extra dimension cannot cause new non-analyticities in the partition function and
finite temperature quantum systems in $d$ dimensions are very similar to classical systems
in $d$ dimensions.

This quantum to classical mapping can also be inverted. Taking a codimension $1$ slice out
of a $D$-dimensional classical partition function, one obtains a transfer matrix which can
be interpreted as the exponential of a quantum hamiltonian acting on the Hilbert space of a
$d=D-1$ dimensional quantum system. Hence, methods for targetting quantum ground states in
$d$ dimensions can also be used to study problems in $(d+1)$-dimensional classical
statistical mechanics.

It is clear that the quantum to classical mapping is not specific to the quantum Ising model
and can be applied to any hamiltonian. The path integral representation of the partition
function fits within the same scheme, and only differs in the fact that the limit
$\epsilon\to 0$ is taken such that the additional dimension becomes continuous. This is
particularly natural if also the spatial dimensions of the quantum dimension are continuous,
i.e. if we have a quantum field theory. In this case, a $D=d+1$ dimensional classical
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
be designed, as samples might annihilate each other. This is known as the **sign problem**.

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

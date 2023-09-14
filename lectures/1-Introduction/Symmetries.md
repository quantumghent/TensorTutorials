Symmetries in quantum many-body physics
=======================================

The goal of this section is to give a very gentle introduction to the concept of symmetries in quantum many-body physics. The general mathematical framework of symmetries in physics (or at least the framework we will restrict to) is that of group - and representation theory. Our goal is not to take this framework as given and illustrate it, but rather to first discuss a couple of important applications of symmetries in the context of some concrete models and gradually build up to the more general framework. We will finish our discussion with an outlook to generalizations of the framework presented here. It goes without saying that we will only scrape the surface of this vast topic. The interested reader is refered to the immense literature on this topic, or to a more specialized course.

## Examples and applications
### Symmetry breaking, order parameters and phases
Recall the one-dimensional transverse field Ising model defined above. It degrees of freedom are qubits ordered on a one-dimensional lattice, and its Hamiltonian reads
$$
    H = -\sum_{i} \sigma^z_i\sigma^z_{i+1} -h_x\sum_i\sigma^x_i.
$$
Let us not wory too much about boundary conditions consider the infinite chain for now. Besides the obvious translation symmetry, which we will discuss below, this model is also is invariant under flipping all spins simultanuously in the z-direction, ie. in the Pauli Z basis: $\ket{\uparrow}\leftrightarrow\ket{\downarrow}$. This is clear from the Hamiltonian as the energy of the first term only depends on the neighbouring spins being (anti-)aligned, which is clearly spin flip-invariant. The second spin is trivially invariant as this models an external magnetic field which is orthogonal to the Z-direction.

This spin flip is implemente, or more correctly represented, by the unitary operator $P=\bigotimes_i \sigma^x_i$. Notice that $P^2=1$ in accordance with our intuition that flipping all the spins twice is equivalent with leaving all spins untouched. The fact that this operator represents a symmetry of the model then translates to $[H,P]=0$, or equivalently $P^\dagger HP=H$.

Even though the Hamiltonian has the symmetry regardless of the value of the parameter $h_x$, you might know from a previous course that the ground state or states are not necessarily invariant under the symmetry, a phenomenom know as spontanuous symmetry breaking (SSB) or symmetry breaking for short. Let us investigate the ground state subspace of the trnasverse field Ising model in the extremal case of vanishing and infinite magnetic field.

- $h_x\rightarrow \infty$ In this case the model effectly reduces to a paramagnet. The uniqe ground state is the product state $\ket{\Psi_+}=\ket{+}^{\otimes N}$ where $\ket{+}=\frac{1}{\sqrt{2}}(\ket{\uparrow}+\ket{\downarrow})$ is the eigenvector of $\sigma^x$ with the highest eigenvalue, namely 1. Notice that this state is invariant under the symmetry operator $P$, $P\ket{\Psi_+}=\ket{\Psi_+}$. In other words, the ground state in this case is symmetric.
- $h_x=0$ In this case the energy is minimzed by aligning all the spins and behaves as a classical ferromagnet. Obviously, two distinct ground states are $\ket{\Psi_\uparrow}=\ket{\uparrow\uparrow...\uparrow}$ and $\ket{\Psi_\downarrow}=\ket{\downarrow\downarrow...\downarrow}$. Contrary as in the previous case they span a two-dimensional ground state subspace, and these states are not symmetric. In fact, under the action of $P$ they get mapped onto the other: $P\ket{\Psi_\uparrow}=\ket{\Psi_\downarrow}$ and vice versa. The ground state in this case is thus symmetry broken.

Since the ground state degeneracy is necessarily an integer, it is clear that it can not change smoothly from 2 to 1 when the magnetic field is slowly turned on from $h_x=0\rightarrow \infty$. Therefore the Ising model for small $h_x$ and large $h_x$ are said to belong to different phases, and for some finite value of $h_x$ a phase transition where the ground state degeneracy changes abruptly is expected to take place. As it turns out, this change happens for $h_x=1$, at which the Ising model becomes critical.

Inspired by the credo of symmetry we can introduce an opetor which is a local order parameter that is a witness of the phase in which we are in. In the case of the Ising model this order paramter is the local magnetisation on every site: $O=\sum_i\sigma^z_i$. It is clear that this order parameter anticommutes with the symmetry, $P^\dagger OP=-O$, from which it follows that in the symmetric phase the expectiation value of the order paramter vanishes, $\braket{\Psi_+|O|\Psi_+}=0$, while in the ferromagnetic phase $\braket{\Psi_\uparrow|O|\Psi_\uparrow}>0$ and $\braket{\Psi_\downarrow|O|\Psi_\downarrow}<0$. Notice however that for the latter we could also have chosen the ground state $\ket{\Psi_\uparrow}+\ket{\Psi_\downarrow}$ in which case the expectation value of $O$ becomes 0. This can be remedied by first adding a small symmetry breaking term $\lambda\sum_i\sigma^z_i$ to the Hamiltonian which depending on the sign of $\lambda$ selects one of the ground states $\ket{\Psi_{\uparrow/\downarrow}}$ after which the limit $\lambda\rightarrow 0$ is taken.

The synopsis of this example is thus the following. Symmetries in quantum many-body physics (but also in single-particle quantum mechanics) are represented by unitary operators which are closed under multiplication. Depending on the paramters in the Hamiltonian, part of these symmetries can be broken by the ground state subspace, and this pattern of symmetry breaking is a hallmark feature of different phases of the model. Different phases can be probed by a local order paramter which does not commute with the symmetries.

### Noether and conserved quantitites

You might remember Noether's theorem from a course on field theory. It states that every continuous symmetry of a system (often defined via its Lagrangian) gives rise to a conserved current. In the context of quantum physics Noether's theorem becomes almost trivial and states that the expectation value of every operator that commutes with the Hamiltonian has a conserved expectation value:
$$
    [H,O] = 0 \implies \frac{d}{dt}\braket{\Psi(t)|O|\Psi(t)} = 0.
$$
The proof is almost trivial and is left as a simple exercise.

The simplest example of this principle is obviously the Hamiltonain that trivially commutes with itself. The consequence is that the expectation value of the total energy is conserved.

Let us consider another non-trivial example to illustrate the implications of this theorem. Recall the spins s XXZ Heisenberg model whose Hamiltonian reads
$$
    H = -J\sum_i S^x_iS^x_{i+1}+S^y_iS^y_{i+1}+\Delta S^z_iS^z_{i+1}.
$$
The spin operators are $2s+1$-dimensional and satisfy the $\mathfrak{su}(2)$ commutation relations
$$
    [\sigma^a_i,\sigma^b_j]=i\delta_{i,j}\sum_c \varepsilon_{abc}S^c_i 
$$
Let us define the total spin
$$
    S^a = \sum_i S^a_i.
$$
From a direct computation it follows that in the case where $\Delta=1$, and the model thus reduces to the Heisenberg XXX model, H commutes with all $S^a$, $[H,S^a]=0$, $a=x,y,z$. However, when $\Delta\neq 1$ only the Z component $S^z$ commutes with H, $[H, S^z]=0$. Notice the difference with the Ising model where the same symmetry was present for all values of $h_x$.

This means that in the $\Delta=1$ case the Hamiltonian is symmetric under the full $SU(2)$ (half integer s) or $SO(3)$ (integer s) symmetry (see below), whereas when $\Delta\neq 1$ only an $SO(2)\neq U(1)$ symmetry generated by $S^z$ is retained.

According to Noether the Heisenberg model thus has conserved quantities. Regardless of $\Delta$ the Z component of the total spin is conserved, and for $\Delta$ all components of the total spin are conserved. In particular this means that the eigenvalue $M_z$ of $S^z$ and $S(S+1)$ of $\vec{S}\cdot\vec{S}$ are good quantum numbers to label the eigenstates of the Heisenberg Hamiltonian.

## Group - and representation theory
Motivated by the examples from above, we will gently introduce some notions of group - and representation theory that form the backbone of a general theory of symmetries.
### Group theory
First of all notice that a model can have a finite or infinite (discrete or continuous) number of symmetries. Clearly, the spin flip symmetry of the Ising model consists of only one non-trivial symmetry operation, namely flipping all spins. The operator carrying out this transformation is $\bigotimes_i\sigma^x_i$ The XXZ model however has a continuous symmetry, namely rotations around the Z axis, that is implemented via $\exp(i\theta S^z)$, $\theta\in[0,2\pi)$.

These symmetries can be composed or multiplied to form a new symmetry operation. Take for example flipping all the spins. Flipping all spins twice results in not flipping any spins at all, which is trivially also a symmetry of the Hamiltonian. Consider for example also the $U(1)$ symmetry of the XXZ model. First rotating over $\theta_2$ and then over $\theta_1$ gives a new rotation over $\theta_1+\theta_2$: $\exp(i\theta_1 S^z)\exp(i\theta_2 S^z)=\exp(i(\theta_1+\theta_2) S^z)$. This leads to the first part of the definition of what a group is.

1. A group $G$ is a set $G=\{g_1,g_2,...\}$ endowed with a multiplication $\cdot:G\times G\rightarrow G$. There exists an identity $1\in G$ for the multiplication such that $1\cdot g=g\cdot 1=g, \forall g\in G$.

Note that this multiplication is not necessarily abelian. A simple example is the full $SU(2)$ symmetry of the XXX model defined above. However, the composition of symmetries is still associative:

2. For all group elements $g,h,k$ we have that $g(hk)=(gh)k$.

A property we also would like to formalize is the fact that every symmetry can be undone. Take for example a $U(1)$ rotation $\exp(i\theta S^z)$, if we compose it with the opposite rotation $\exp(i(2\pi-\theta) S^z)$ we get the identity.

3. Every group element $g$ has a unique inverse $g^{-1}$: $gg^{-1}=g^{-1}g=1$.

Together 1. 2. and 3. constitute the definition of a group. Before mentioning some examples let us also introduce the concept of a subgroup. As the name suggests, a subgroup is a subset of a group which itself constitutes a group.

The concept of subgroups lies at the heart of symmetry breaking. Recall that in the ferromagnetic limit, the Ising model break the spin flip symmetry. In Landau's paradigm we say that the pattern of symmetry breaking is $\mathbb{Z}_2\rightarrow \{1\}$ (see below for an explanation of the notation). In other words, the full symmetry is broken in the ferromagnetic phase. More generally, a theory with a $G$ symmetry can undergo a pattern of symmetry breaking $G\rightarrow H$ where $H$ is a subgroup of $G$. The meaning of this symbolic expression is that in this particular case the ground states keep an H symmetry, and the ground state degeneracy is $|G|/|H|$.
#### Examples

- $\mathbb{Z}_N$ is the additive group integers modulo $N$. The group elements are the integers $\{0,1,...,N-1\}$ and hence it is clearly a finite group. The group multiplication is addition modulo $N$. In particular, the spin flip symmetry from above corresponds to the group $\mathbb{Z}_2$. Notice that for all $N$ $\mathbb{Z}_N$ is abelian.
- Another abelian group is $U(1)$. $U(1)=\left\{z\in\mathbb{C}:|z|^2 = 1\right\}$, with group multiplication the multiplication of complex numbers. Note we encountered this group in the XXZ model as being the rotations around the Z axis.
- $SU(2)$ is the group of unimodular unitary $2\times 2$ matrices:
    $$
        SU(2) := \left\{U\in\mathbb{C}^{2\times 2}|\det U=1, UU^\dagger=U^\dagger U=\mathbb{I}\right\}.
    $$
    The group multiplication is given by group multiplication. Similarly, one defines $SU(N),N\geq 2$. Note that none of these groups are abelian.

- The 3D rotation group or special orthogonal group $SO(3)$ is the group of real $3\times 3$ orthogonal matrices with unit determinant:
	$$
		SO(3) := \left\{M\in\mathbb{R}^{3\times 3}|MM^T=M^TM=\mathbb{I},\det M=1\right\}.
	$$ 
    Similarly, one defines $SO(N),N\geq 2$. Note that only $SO(2)$ is abelian.
### Representation theory
#### Examples
## Outlook and generalizations
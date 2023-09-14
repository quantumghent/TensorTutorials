Matrix Product States
=====================

1. Matrix Product States
	- MPS as a polynomial Ansatz
	- entanglement and area laws
	- gauge freedom
	- expectation values
	- correlators

# Introduction
Consider a quantum system of $N$ particles, where each particle can be in in $d$ different states (we call this a system of qudits). The total Hilbert space is of this system is spanned by $d^N$ basis states. This means that if we want to track all possible states of our system, we need to track an exponential amount of states. This is called **the dimensionality curse** and is one of the prime reasons why studying systems with many particles is difficult. To circumvent this problem we can only consider a small corner of the full Hilbert space and hope that the essential physics of our system is represented by this small corner. One such small corner is parametrized by **matrix product states (MPS)**, which has been succesful in describing the physics of many different interacting quantum systems.

# Matrix product states through succesive SVD
We can represent every basis state of our $N$ qudit Hilbert space as $\ket{s_1s_2...s_N}$, where $s_i\in \{0,...,d-1\}$. To specify a general quantum state can be expressed as a linear combination of these basis states. We can write
```{math}
\ket{\psi} = \sum_{s_1,s_2,...,s_N} C_{s_1,s_2,...,s_N}\ket{s_1s_2...s_N}.
```
So for a general state we need the values of all the different $C_{s_1,s_2,...,s_N}$, hence we need to track $d^N$ complex numbers. Let us perform some tricks to handle this exponential number of coefficients. First we are going to isolate the first index by interpreting the coefficients as a $d$ by $d^(N-1)$ matrix. Next we can write this matrix as
```{math}
C_{s_1,s_2,...,s_N} = \sum_{k\in{0,D-1}} U_{s_1,k}C_{k,s_2,...,s_N}',
```
where $U$ is a $d$ by $D$ matrix and $C'$ is a $D$ by $d^(N-1)$ matrix. One can find such a matrix by performing an **Singular-Value-Decomposition (SVD)** on the original matrix $C$. A SVD decomposes the original matrix $C$ as $U\Sigma V^\dagger$ with $U$ and $V$ unitaries and $\Sigma$ a diagonal matrix containing the singular values of $C$. $\Sigma$ will be a $d$ by $d^(N-1)$ matrix with rank($C$) non-zero elements. One can then define $C'$ as $\Sigma V^\dagger$.

So far it seems that we made things worse, we now need to store $d^2$ elements of $U$ and $d^N$ elements of $C'$. We can make things better by making a **low rank approximation** of the state instead of representing the full state. This low rank approximation can be made by **truncating** $\Sigma$ by only retaining the largest $D$ singular values. This truncates $C'$ to a $D$ by $d^(N-1)$ matrix. $D$ is called the bond dimension and increasing $D$ will give you a better approximation of the original state. 

By isolating the first index of our original tensor $C$ and by performing an SVD and truncating, we have gained a low rank approximation of the original tensor $C$. Now we can use the same trick on index 2. We can do this by writing
```{math}
C_{s_1,s_2,...,s_N} = \sum_{k,l\in{0,D-1}} U^1_{s_1,k}U^2_{k,s_2,l}C_{l,s_3,...,s_N}'.
``` 
The above equation is an equality if $D=rank(C)$, but again we are free to choose a lower bond dimension to make an approximation of $C$. By performing succesive SVD's we decompose/approximate our original tensor $C$ as a "**product of matrices**". We hence write
```{math}
C_{s_1,s_2,...,s_N} = \sum_{k,l,...,z\in{0,D-1}} U^1_{s_1,k}U^2_{k,s_2,l}...U^N_{z,s_2,l}.
``` 
Notice that our original state is now given by a contraction of $N$ three leg tensors. A state that can be written in this way is called a **matrix product state**. In principle any state can be written as a MPS, but one would need a large bond dimensions for an arbitray state. The true power lies in using MPS as a variational ansatz. One can fix the bond dimension $D$ and find the optimal MPS within this subspace. Furthermore many systems are translationally invariant meaning that one can use the same tensor $U$ on every site, which gives a huge compression of information. 

# Matrix product states through projected entangled pairs

Another approach to construct a MPS is by using projected entangled pairs. Consider a chain with $N$ sites. On every site we place 2 qudits, giving us a local Hilbert space of $\mathbb{C}^D\otimes\mathbb{C}^D$ on each site. Now we maximally entangle each qudit with the nearest qubit on a neighbouring site. We have now a chain of entangled pairs, where the full quantum state of the chain is given by 
```{math}
\ket{\psi} = \bigotimes_i=1^{N}\sum_{j=0}^{D-1}\ket{jj}.
```
Next we construct a map $f:\mathbb{C}^D\otimes\mathbb{C}^D\xrightarrow \mathbb{C}^d$ acting on a local site. This map projects the $D^2$-dimensional "virtual" Hilbert space onto a "physical" $d$-dimensional Hilbert space. Acting with $f$ on each site gives a $d$-dimensional quantum chain with entanglement shared between neighbouring sites. Since the state is fully determined by the map $f$, which can be represented as a $D$ by $d$ by $D$ three leg tensor it is clear that this construction indeed gives rise to a MPS. One advantage of this construction is that it allows for easy generalization to higher dimensions.

# Entanglement of Matrix product states
One important question is which states can be represented by an MPS and what are their properties. To answer this question we will look at the entanglement structure of a MPS. Consider a MPS with a fixed bond dimension $D$. We can choose any bond of the MPS to make a bipartition of the MPS, where all sites to the left of the bond are in one partition and all sites to the right of the partition are in the other partition. One can now calculate the entanglement entropy $S$ between these two partitions. To calculate this we will use the SVD decomposition. If we perform an SVD along the bond, then the singular values will completely determine the entanglement structure between the two partitions. The entanglement entropy is then given by
```{math}
S = -\sum_i\abs{\sigma_i}^2\log{(\abs{\sigma_i}^2)}.
```
Each bond we cut to make the bipartition will contribute to the total entanglement entropy. To bipartition a MPS we only need to cross one bond giving us a constant $S$. In higher dimensional systems the amount of bonds that you need to cross to bipartition the system will scale with the area of the partition, hence the entanglement entropy also scales with the area. This is scaling relation is called an area law. This leads us to the following statement: "Any state with an area law for the entanglement entropy can be arbitrarily good approximated by an MPS wiht **poly**$(N)$ bond dimension.". The small corner that MPS parametrizes is thus given by the weakly entangled states that are governed by an area law. This seems very restrictive however it appears that many of the relevant states for different physical systems fall within this class.
# Gauge freedom of MPS

# Expectation values of Matrix product states

# Correlation functions of Matrix product states.
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
C_{s_1,s_2,...,s_N} = \sum_{k\in{0,\chi-1}} U_{s_1,k}C_{k,s_2,...,s_N}',
```
where $U$ is a $d$ by $\chi$ matrix and $C'$ is a $\chi$ by $d^(N-1)$ matrix. One can find such a matrix by performing an **Singular-Value-Decomposition (SVD)** on the original matrix $C$. A SVD decomposes the original matrix $C$ as $U\Sigma V^\dagger$ with $U$ and $V$ unitaries and $\Sigma$ a diagonal matrix containing the singular values of $C$. $\Sigma$ will be a $d$ by $d^(N-1)$ matrix with rank($C$) non-zero elements. One can then define $C'$ as $\Sigma V^\dagger$.

So far it seems that we made things worse, we now need to store $d^2$ elements of $U$ and $d^N$ elements of $C'$. We can make things better by making a **low rank approximation** of the state instead of representing the full state. This low rank approximation can be made by **truncating** $\Sigma$ by only retaining the largest $\chi$ singular values. This truncates $C'$ to a $\chi$ by $d^(N-1)$ matrix. $\chi$ is called the bond dimension and increasing $\chi$ will give you a better approximation of the original state. 

By isolating the first index of our original tensor $C$ and by performing an SVD and truncating, we have gained a low rank approximation of the original tensor $C$. Now we can use the same trick on index 2. We can do this by writing
```{math}
C_{s_1,s_2,...,s_N} = \sum_{k,l\in{0,\chi-1}} U^1_{s_1,k}U^2_{k,s_2,l}C_{l,s_3,...,s_N}'.
``` 
The above equation is an equality if $\chi=rank(C)$, but again we are free to choose a lower bond dimension to make an approximation of $C$. By performing succesive SVD's we decompose/approximate our original tensor $C$ as a **product of matrices**. We hence write
```{math}
C_{s_1,s_2,...,s_N} = \sum_{k,l,...,z\in{0,\chi-1}} U^1_{s_1,k}U^2_{k,s_2,l}...U^N_{z,s_2,l}.
``` 
Notice that our original state is now given by a contraction of $N$ three leg tensors. A state that can be written in this way is called a **matrix product state**. In principle any state can be written as a MPS, but one would need an large bond dimensions for an arbitray state. The true power lies in using MPS as a variational ansatz. One can fix the bond dimension $\chi$ and find the optimal MPS within this subspace. Furthermore many systems are translationally invariant meaning that one can use the same tensor $U$ on every site, which gives a huge compression of information.

# Matrix product states through projected entangled pairs

# Entanglement of Matrix product states

# Expectation values of Matrix product states

# Correlation functions of Matrix product states.
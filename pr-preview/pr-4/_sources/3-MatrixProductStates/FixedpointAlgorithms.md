# Fixed-Point algorithms

In this section we introduce two algorithms for approximating the ground state of local gapped Hamiltonians using matrix product state techniques. Approximating ground states in a variational manner boils down to minimizing
```{math}
    \min_{\ket{\psi}\in D} \frac{\braket{\psi|H|\psi}}{\braket{\psi|\psi}},
```
over a restricted class of states $D$. For simplicity, we will restrict ourselves to Hamiltonians that consist of a sum of local nearest-neighbor terms, i.e.
```{math}
    H = \sum_i h_{i,i+1}.
```
It is understood that $h_{i,i+1}$ only acts non-trivially on sites $i$ and $i+1$. In the algorithms discussed below we optimize over matrix product states of a fixed finite bond dimension. In the first algorithm known as DMRG (density matrix renormalization group) the states we consider are finite MPS, whereas the second algorithm VUMPS (variational uniform matrix product state algorithm), as the name suggests, optimizes over uniform MPS. Hence, VUMPS enables direct optimization in the thermodynamic limit, without breaking translation invariance.

Our exposition of DMRG closes follows the one in {cite}`bridgeman2017handwaving`, and that of VUMPS closely follows the excellent set of lecture notes {cite}`vanderstraeten2019tangentspace`.

## DMRG
Starting from a random MPS ansatz, DMRG tries to approximate the ground state by sequentially optimizing over all the MPS tensors one by one and sweeping through the chain, until convergence is reached. Let us discuss this algorithm in a bit more detail step by step.

### Algorithm

Let us consider a random ansatz, by taking random tensors $\{A_1,A_2,...,A_L\}$, $L$ being the number of sites. Fixing all tensors but the one at site $i$, the local tensor $A_i$ is updated according to

(5.4, Hand waving)

Though seemingly daunting we can turn this problem in a simple eigenvalue problem by making full use of the mixed gauge. By bringing all tensors on the right of $A_i$ in the right canonical form and those to the left in left canonical form, the denominator simply becomes $\braket{A_i|A_i}$.

The numerator then reduces to 

(Numerator 5.7, Hand waving),

So that in the end the problem we have to solve is to find the minimum eigenvector of $\mathcal{H}_i$, and this repeatedly for every site, sweeping back and forth through the chain. Notice however that DMRG manifestly breaks translation invariance by updating one tensor at a time. As we will see, VUMPS does not suffer from this artefact.

From this brief explanation it should be clear that DMRG is a surprisingly simple algorithm. Nevertheless DMRG has proven itself time and time again, and is the most successful algorithm for variationally approximating the ground state of local gapped (1+1)d Hamiltonians. DMRG is implemented in MPSKit and can be called by `DMRG()`.

### Example

Let us illustrate the use of DMRG in MPSKit by approximating the ground state of the transverse field Ising model. The Ising model is implemented in MPSKitModels as follows
```{math}
    H = -J\left(\sum_{<i,j>} Z_i Z_j + \sum_i h_x X_i + h_z Z_i\right),
```
where we are free to choose the parameters $J$, $h_x$ and $h_z$, and $X$ and $Z$ are the generators of $\mathfrak{su}(2)$, and thus differ from the usual Pauli matrices by a factor of $\frac{1}{2}$.

Let us consider 16 lattice sites, bond dimension 12, open boundary conditions and let's stick to the default critical values of $h_x=0.5$ and $h_z=0$. Finding the ground state using DMRG then only takes a handful of iterations!

```{code-cell} julia
    L = 16; # Length spin chain
    D = 12; # Bond dimension

    H = transverse_field_ising();

    algorithm = DMRG(); # Summon DMRG
    Ψ = FiniteMPS(L,ℂ^2,ℂ^D); # Random MPS ansatz with bond dimension D
    Ψ₀,_ = find_groundstate(Ψ,H,algorithm);
```

## VUMPS

As mentioned above, VUMPS optimizes uniform MPS directly in the thermodynamic limit. Since the total energy becomes unbounded in this limit, our objective should be to rather minimize the energy density. Pictorially:

(105, Tangent space methods).

The VUMPS algorithm offers the advantage of global optimalization by design, since the algorithm, contrary to DMRG, does not rely on updates of local tensors.

Given a local Hamiltonian of the form mentioned above and an intial random uniform MPS defined by $\{A_L, A_R,C\}$, VUMPS approximates the ground state by finding an approximate solution to the fixed-point equations
```{math}
    A_C' = H_{A_C}(A_C), \\
    C' = H_C(C), \\
    A_C = A_LC = CA_R.
```
A detailed derivation that these equations characterize the variational minimum in the manifold of uniform MPS is beyond the scope of these notes, but see {cite}`vanderstraeten2019tangentspace`.

In these equations the effective Hamiltonian $H_{A_C}$ acting on $A_C$ is given by

(131, Tangent space methods).

Where $h$ is the local Hamiltonian term, and $H_{L/R}$ are environments that capture the contributions of all local terms acting on the left/right. The other effective Hamiltonian $H_C$ acting on the bond tensor $C$ can then be represented by the network

(132, Tangent space methods).

The last equation then simply states that $C$ intertwines the left - and right-orthonormal form of the tensor $A$.


### Algorithm

Let us now explain step-by-step how VUMPS finds an approximate solution to the fixed-point equations in an iterative way.


1. We initialize the algorithm with the random guess $\{A_L, A_R,C\}$, and chose a tolerance $\eta$.

2. We first solve the first two eigenvalue equations
    ```{math}
        A_C = H_{A_C}(A_C), \\
        C = H_C(C),
    ```
    using for example an Arnoldi algorithm with the previous approximations of $A_C$ and $C$ as initial guess. This yields two tensors $\tilde A_C$ and $\tilde C$.

3. From $\tilde A_C$ and $\tilde C$ we compute $\tilde A_L$ and $\tilde A_R$ that minimize following two-norms
    
    ```{math}
        \epsilon_L = \min_{A_L^\dagger A_L=1} ||\tilde A_C-\tilde A_L\tilde C||_2, \\
        \epsilon_R = \min_{A_R A_R^\dagger=1} ||\tilde A_C-\tilde C\tilde A_R||_2,
    ```
    and thus approximately solve the last equation. Note that the minimum is taken over respectively left - and right isometric matrices. We comment below on the analytic soltuion of these equations and how this analytic solution can be approximated efficiently.

4. Update $A_L\leftarrow\tilde A_L$, $A_R\leftarrow\tilde A_R$ and $C\leftarrow\tilde C$.

5. Evaluate $\epsilon=\max(\epsilon_L,\epsilon_R)$ and repeat until $\epsilon$ is below the tolerance $\eta$.

Let us finally comment on solving the minimization problem to approximate $\tilde A_{L/R}$.

A beautiful result in linear algebra states that the minimum is exactly given by $\tilde A_L=U_lV_l^\dagger$ where $U_l$ and $V_l$ are the isometries arising from the singular value decomposition of $\tilde A_C\tilde C^\dagger=U_l\Sigma_lV_l^\dagger$, and similarly $\tilde A_R=U_rV_r^\dagger$, where $\tilde C^\dagger\tilde  A_C=U_r\Sigma_rV_r^\dagger$. Even though this approach will work well for the first iteration steps, this might not be the best solution close to convergence. When approaching the exact solution $A^s_C=A^s_LC=CA^s_R$ the singular values in $\Sigma_{l/r}$ become really small so that in finite precision arithmetic the singular vector in the isometries $U_{l/r}$ and $V_{l/r}$ are poor approximations of the exact singular vectors. A robust and close to optimal solution turns out to be
```{math}
    \tilde A_L = U^l_{A_C}(U^l_C)^\dagger,\qquad \tilde A_R = (U^r_C)^\dagger U^r_{A_C},
```
where the $U$'s are the unitaries appearing in the polar decomposition of
```{math}
    \tilde A_C = U^l_{A_C}P^l_{A_C},\qquad \tilde C = U^l_CP^l_C,\\
    \tilde A^r_C = P^r_{A_C}U^r_{A_C},\qquad \tilde C = P^r_CU^r_C.
```

### Example

Let us demonstrate the algorithm using MPSKit by estimating the ground state energy density of the spin 1 XXX model. The VUMPS algorithm is called in the same way as we called DMRG. We initialize a random initial MPS with bond dimension 12 and physical dimension 3 (because the spin 1 representation of SU(2) is $2\cdot1+1=3$-dimensional). Obviously we don't have to specify a system size because we work directly in the thermodynamic limit.

```{code-cell} julia
    H = xyz();

    Ψ = InfiniteMPS(ℂ^3,ℂ^D);
    algorithm = VUMPS();
    Ψ₀, envs = find_groundstate(Ψ,H,algorithm);
```

It takes about 30 iterations and a second or two to reach convergence. Let us gauge how well the ground state energy density was approximated by calling

```{code-cell} julia
    expectation_value(Ψ₀,H)
```

The value we obtain here is to be compared with the quasi-exact value -1.401 484 038 971 2(2) obtained in {cite}`haegeman2011time`. As you can see, even with such a small bond dimension we can easily approximate the ground state energy up to 3 decimals.
# Fixed point algorithms

In this section we introduce two algorithms for approximating the ground state of local gapped Hamiltonians using matrix product state techniques. Approximating ground states in a variational manner boils down to minimizing
$$
    \min_{\ket{\psi}\in D} \frac{\braket{\psi|H|\psi}}{\braket{\psi|\psi}},
$$
over a restricted class of states $D$. In the algorithms discussed below we optimize over matrix product states of a fixed finite bond dimension. In the first algorithm going by the name of DMRG (density matrix renormalization group) the states we consider are finite MPS, whereas the second algorithm VUMPS (variational uniform matrix product state algorithm), as the name suggests, optimizes over uniform MPS.

## DMRG
Starting from a random MPS ansatz, DMRG tries to approximate the ground state by sequentially optimizing over all the MPS tensors one by one and sweeping through the chain, until convergence is reached. Let us discuss this algorithm in a bit more detail step by step.

Fixing all tensors but the one at site $i$, the local tensor $A_i$ is updated according to

(5.4)

Though seemingly daunting we can turn this problem in a simple eigenvalue problem by making full use of the center gauge. By bringing all tensors on the right of $A_i$ in the right canonical form and those to the left in left canonical form, the denominator simply becomes $\braket{A_i|A_i}$.

The numerator then reduces to 

(5.7),

So that in the end the problem we have to solve is to find the minimum eigenvector of $\mathcal{H}_i$, and this repeatedly for every site, sweeping back and forth through the chain

From this brief explanation it should be clear that DMRG is a surprisingly simple algorithm. Nevertheless DMRG has proven itself time and time again, and is the most succesful algorithm for variationally approximating the ground state of local gapped (1+1)d Hamiltonians. DMRG is implemented in MPSKit and can be called by `DMRG()`.

Let us illustrate the use of DMRG in MPSKit by approximating the ground state of the transverse field Ising model. The Ising model is implemented in MPSKitModels as follows
$$
    H = -J\left(\sum_{<i,j>} Z_i Z_j + \sum_i h_x X_i + h_z Z_i\right),
$$
where we are free to choose the parameters $J$, $h_x$ and $h_z$.

Let us consider 16 lattice sites, bond dimension 6, open boundary conditions and let's stick to the default critical values of $h_x=.5$ and $h_z=0$. Finding the ground state using DMRG then only takes a handful of iterations!

```{code-cell} julia
    L = 16; # Length spin chain
    D = 6; # Bond dimension

    H = transverse_field_ising();

    algorithm = DMRG(verbose=false); # Summon DMRG
    Ψ = FiniteMPS(L,ℂ^2,ℂ^D); # Random MPS ansatz with bond dimension D
    Ψ₀,_ = find_groundstate(Ψ,H,algorithm);
```
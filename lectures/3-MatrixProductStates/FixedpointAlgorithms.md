# Fixed-Point algorithms

In this section we introduce two algorithms for approximating the ground state of local gapped Hamiltonians using matrix product state techniques. Approximating ground states in a variational manner boils down to minimizing
$$
    \min_{\ket{\psi}\in D} \frac{\braket{\psi|H|\psi}}{\braket{\psi|\psi}},
$$
over a restricted class of states $D$. In the algorithms discussed below we optimize over matrix product states of a fixed finite bond dimension. In the first algorithm known as DMRG (density matrix renormalization group) the states we consider are finite MPS, whereas the second algorithm VUMPS (variational uniform matrix product state algorithm), as the name suggests, optimizes over uniform MPS. Hence, VUMPS allows for direct optimization in the thermodynamic limit and does not break translation invariance.

## DMRG
Starting from a random MPS ansatz, DMRG tries to approximate the ground state by sequentially optimizing over all the MPS tensors one by one and sweeping through the chain, until convergence is reached. Let us discuss this algorithm in a bit more detail step by step.

Fixing all tensors but the one at site $i$, the local tensor $A_i$ is updated according to

(5.4)

Though seemingly daunting we can turn this problem in a simple eigenvalue problem by making full use of the center gauge. By bringing all tensors on the right of $A_i$ in the right canonical form and those to the left in left canonical form, the denominator simply becomes $\braket{A_i|A_i}$.

The numerator then reduces to 

(5.7),

So that in the end the problem we have to solve is to find the minimum eigenvector of $\mathcal{H}_i$, and this repeatedly for every site, sweeping back and forth through the chain.

Notice however that DMRG manifestly breaks translation invariance by updating one tensor at a time. VUMPS does not suffer from this artefact.

From this brief explanation it should be clear that DMRG is a surprisingly simple algorithm. Nevertheless DMRG has proven itself time and time again, and is the most succesful algorithm for variationally approximating the ground state of local gapped (1+1)d Hamiltonians. DMRG is implemented in MPSKit and can be called by `DMRG()`.

Let us illustrate the use of DMRG in MPSKit by approximating the ground state of the transverse field Ising model. The Ising model is implemented in MPSKitModels as follows
$$
    H = -J\left(\sum_{<i,j>} Z_i Z_j + \sum_i h_x X_i + h_z Z_i\right),
$$
where we are free to choose the parameters $J$, $h_x$ and $h_z$.

Let us consider 16 lattice sites, bond dimension 12, open boundary conditions and let's stick to the default critical values of $h_x=.5$ and $h_z=0$. Finding the ground state using DMRG then only takes a handful of iterations!

```{code-cell} julia
    L = 16; # Length spin chain
    D = 12; # Bond dimension

    H = transverse_field_ising();

    algorithm = DMRG(); # Summon DMRG
    Ψ = FiniteMPS(L,ℂ^2,ℂ^D); # Random MPS ansatz with bond dimension D
    Ψ₀,_ = find_groundstate(Ψ,H,algorithm);
```

## VUMPS

As mentioned above, VUMPS optimizes uniform MPS directly in the thermodynamic limit. By construction VUMPS respects translational invariance. As will become clear below, VUMPS can be thought of as a generalization of DMRG and we will consider this point of view in these notes.

Let us demonstrate the algorithm by estimating the ground state energy density of the spin 1 XXX model. The VUMPS algorithm is called in the same way as we called DMRG. We initialize a random initial MPS with bond dimension 12 and physical dimension 3 (because the spin 1 representation of SU(2) is $2\cdot1+1=3$-dimensional). Obviously we don't have to specify a sytsem size because we work directly in the thermodynamic limit.

```{code-cell} julia
    H = xyz();

    Ψ = InfiniteMPS(ℂ^3,ℂ^D);
    algorithm = VUMPS();
    Ψ₀, envs = find_groundstate(Ψ,H,algorithm);
    expectation_value(Ψ₀,H);
```

It takes about 30 iterations and a couple of seconds to reach convergence. Let us gauge how well the ground state energy density was approximated by calling

```{code-cell} julia
    expectation_value(Ψ₀,H);
```

The value we obtain here is to be compared with the quasi-exact value -1.401 484 038 971 2(2) obtained in [TDVP Jutho].
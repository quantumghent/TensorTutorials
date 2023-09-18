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



## Quantum Mechanics in 1+1D

One of the most important techniques for MPOs is the ability to write a sum of local
operators in MPO-form. The resulting operator has a very specific structure, and is often
referred to as a "Jordan block" MPO.

### Jordan Block MPOs

For example, if we consider the transverse field Ising model,

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

### Expectation Values

 

### Jordan MPOs in the Thermodynamic Limit






## Quasi-1D Systems



2. Matrix Product Operators
	- Statistical mechanics as an MPO
	- Jordan-block form of MPO
	- expectation values
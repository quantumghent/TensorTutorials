https://quantumghent.github.io/TensorTutorials/

# Tutorial outline
## I. Prerequisites
1. General Introduction
	- who is who?
	- practicalities during the year
	- Julia
	- Git(Hub)
2. Introduction to Quantum Many-Body theory
	- context, history, purpose, relevance
	- second quantization
	- tensor product
	- quantum <-> classical
3. Tensor Networks
	- context, history, purpose, relevance
	- multi-linear algebra
	- graphical notation
	- computational complexity
4. *Symmetries*
	- group theory
	- fermions
	- fusion categories

## II. Matrix Product States
1. Matrix Product States
	- MPS as a polynomial Ansatz
	- entanglement and area laws
	- gauge freedom
	- expectation values
	- correlators
2. Matrix Product Operators
	- Statistical mechanics as an MPO
	- Jordan-block form of MPO
	- expectation values
3. Infinite Matrix Product States
	- gauging revisited
	- expectation values
	- structure factors
4. Basic Algorithms
	- *Simple Update*
	- Trotterization
	- TEBD
5. Applications
	- Transverse-Field Ising Model
	- Classical Ising Model

## III. Tensor Network Algorithms
1. Fixedpoint Algorithms
	- DMRG
	- VUMPS
	- Krylov methods
2. Time-Evolution Algorithms
	- TDVP
	- MPO Compression
	- *Finite Temperature*
3. *Excitations*
	- Finite Excitations
	- Quasiparticle ansatz
4. *Gradient Methods*
	- Quasi-Newton methods
	- Grassmann optimization
	- AD
5. *PEPS*
6. Dynamic Bonddimensions
	- Multi-Site Algorithms
	- MPS Truncation
	- MPO-MPS approximation

---

# Local development:

Run the `build.sh` script to build the pages locally. Optionally, you can pass the `-d` flag to update the dependencies.

```bash
./build.sh [-d]
```
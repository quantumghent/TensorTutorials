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

# Getting Started with Numerics

On this page, there are some links with relevant information for getting started with
numerical computing. We point to some references for the Julia programming language, as well
as some resources for learning about version control software.

## Version Control

Version control software is a tool used in software development (and sometimes in other
fields) to manage and track changes made to a project's source code, documents, or any other
set of files. It allows multiple contributors to work collaboratively on a project, keeping
a history of changes, and facilitating the organization and synchronization of different
versions of the project. The most popular version control system is
[git](https://git-scm.com/), which is a free tool developed by Linus Torvalds in 2005, and
has become the de facto standard in the software development industry.

Again, multiple resources are available for learning about git. From the official website,
the book [Pro Git}(https://git-scm.com/book/en/v2) is a good place to start.

## Julia

In order to get started with Julia, there are many resources already available. The
[official documentation](https://docs.julialang.org/en/v1/) is a good place to start, and a
full _getting started_ exposition can be found for example
[here](https://julia.quantecon.org/intro.html). There is also a
[learning page](https://julialang.org/learning/) that has tutorials on different topics, a
list of books, and much more.

Additionally, there is an active [forum](https://discourse.julialang.org/) for asking
questions, as well as a [slack channel](https://julialang.org/slack/) and a stack overflow page.

## Julia Packages

Julia has a very active open-source community, and many packages are available for different
purposes. These typically have their own documentation, and are hosted on GitHub. An
(incomplete) list of packages that are relevant for this course are given below:

- [TensorKit.jl](https://github.com/Jutho/TensorKit.jl)
- [TensorOperations.jl](https://github.com/Jutho/TensorOperations.jl)
- [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl)
- [OptimKit.jl](https://github.com/Jutho/OptimKit.jl)
- [MPSKit.jl](https://github.com/maartenvd/MPSKit.jl)
- [PEPSKit.jl](https://github.com/quantumghent/PEPSKit.jl)

## Noteworthy Tensor Network Software

There are many additional software libraries available for tensor network computations, or
more generally for quantum physics research. Below you can find an incomplete list of some
of these.

- [ITensor](https://itensor.org/) Julia/C++ library for tensor network calculations
- [TenPy](https://tenpy.readthedocs.io/en/latest/) Python library for tensor network calculations
- [QUIMB](https://quimb.readthedocs.io/en/latest/) Python library for quantum information many-body calculations

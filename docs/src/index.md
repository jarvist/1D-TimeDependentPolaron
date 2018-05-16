These codes simulate the time-evolution of a polaron state. 
Currently they are limited to 1D only. 
The classical response of the lattice is treated as a dielectric material being
polarised by the electron density. 
The electron density is found by forming a tight-binding Hamiltonian for the
lattice, with the site energy modified by the (back-reaction) polarisation
field of the lattice. 

The (classical) lattice degree of freedom and (quantum) electron degree of
freedom are time evolved by direct evaluation of the exponential equations of
motion. 
Initial work has been done to add surface-hopping between different potential
energy surfaces. 

## Installation

These codes require Julia >0.6 . They are structured as a full Julia package,
but are not yet registered with the central METADATA package repository. 

To install, type the following at the Julia REPL:

```
julia> Pkg.clone("git://github.com/Jarvist/TheDancer.jl.git")
```

## Community guidelines

Contributions to the code (extending that which is calculated), or additional
physical systems / examples, are very welcome. 

If you have questions about the software, scientific questions, or find errors,
please create a [GitHub issue](https://github.com/jarvist/TheDancer.jl/issues). 

If you find this package (or snippets, such as the entered and tested
free-energy expressions) useful for your work, please cite the paper as and
when we publish something! 


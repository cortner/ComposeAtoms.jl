# ComposeAtoms

[![Build Status](https://travis-ci.org/cortner/ComposeAtoms.jl.svg?branch=master)](https://travis-ci.org/cortner/ComposeAtoms.jl)

[![Coverage Status](https://coveralls.io/repos/cortner/ComposeAtoms.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cortner/ComposeAtoms.jl?branch=master)

[![codecov.io](http://codecov.io/github/cortner/ComposeAtoms.jl/coverage.svg?branch=master)](http://codecov.io/github/cortner/ComposeAtoms.jl?branch=master)

Visualisation of atomistic configurations, mostly for materials science
(rather than molecules). At the moment, this collection of codes is experimental,
but hopefully it will grow into a useful module.

Requires [JuLIP.jl](https://github.com/libAtoms/JuLIP.jl)

## Current capabilities:

* 2D and quasi-2D configurations via [Compose.jl](https://github.com/GiovineItalia/Compose.jl)


## Examples

* Graphene sheet
```julia
using ComposeAtoms, JuLIP
at = JuLIP.ASE.graphene_nanoribbon(4, 10)
JuLIP.swapxy!(JuLIP.swapyz!(at))
display(at, atradius=0.35)
```

* Edge Dislocation in Bulk Silicon (110 plane)

(the core is not equilibrated!)
```julia
using ComposeAtoms, JuLIP, MaterialsScienceTools
at = MaterialsScienceTools.Silicon.edge_dipole_110(D = 10, L = (20, 10))
set_pbc!(at, (false, false, true))
a0 = rnn("Si")
display_layers(at, [1.9, 0.0], rcut = 1.4 * a0, axis = [0, 16*a0, 6*a0, 19*a0])
```

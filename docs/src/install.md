# Installation

`TwoDimensional` is an [offical Julia package](https://pkg.julialang.org/)
which can be installed in several ways as explained here.


## Using the package manager

At the [REPL of
Julia](https://docs.julialang.org/en/stable/manual/interacting-with-julia/),
hit the `]` key to switch to the package manager REPL (you should get a
`... pkg>` prompt) and type:

```julia
pkg> add TwoDimensional
```

where `pkg>` represents the package manager prompt.

To check whether the `TwoDimensional` package works correctly, type:

```julia
pkg> test TwoDimensional
```

Later, to update to the last version (and run tests), you can type:

```julia
pkg> update TwoDimensional
pkg> build TwoDimensional
pkg> test TwoDimensional
```

If something goes wrong, it may be because you already have an old version of
`TwoDimensional`. Uninstall `TwoDimensional` as follows:

```julia
pkg> rm TwoDimensional
pkg> gc
pkg> add TwoDimensional
```

before re-installing.

To revert to Julia's REPL, hit the `Backspace` key at the `... pkg>` prompt.


## Installation in scripts

To install `TwoDimensional` in a Julia script, write:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/emmt/TwoDimensional.jl", rev="master"));
```

or with `url="git@github.com:emmt/TwoDimensional.jl"` if you want to use `ssh`.

This also works from the Julia REPL.

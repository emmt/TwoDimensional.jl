# TwoDimensional

[![License][license-img]][license-url]
[![Stable][doc-stable-img]][doc-stable-url]
[![Dev][doc-dev-img]][doc-dev-url]
[![Build Status][github-ci-img]][github-ci-url]
[![Build Status][appveyor-img]][appveyor-url]
[![Coverage][codecov-img]][codecov-url]

**TwoDimensional** is a [Julia][julia-url] package which provides useful types
and methods to define and manipulate 2-dimensional points, bounding-boxes and
affine coordinate transforms.

Other related packages:
- [CoordinateTransformations](https://github.com/FugroRoames/CoordinateTransformations.jl)


## Usage

```julia
using TwoDimensional
```

gives you types `AffineTransform{T}`, `Point{T}` and `BoundingBox{T}`
parameterized by the type `T` of their components (`T` must be floating point
for `AffineTransform{T}`).

To avoid conflicts with other packages, you may use/import
`TwoDimensional.Suffixed` which gives you types `AffineTransform2D{T}`,
`Point2D{T}` and `BoundingBox2D{T}` instead, that is with suffix `2D`.

You can also fine tune what you want.  For instance:

```julia
using TwoDimensional: AffineTransform, Point2D
```

## Documentation

Latest documentation is
[here](https://emmt.github.io/TwoDimensional.jl/latest).


## Installation

`TwoDimensional` is an [official Julia package][julia-pkgs-url] so you can
install it from Julia's package manager.  In an interactive Julia session, hit
the `]` key to switch to the package manager REPL (you should get a `... pkg>`
prompt) and type:

```julia
pkg> add TwoDimensional
```

You can also execute the following statements (in a Julia script or from Julia
REPL):

```julia
using Pkg
Pkg.add("TwoDimensional")
```

[doc-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[doc-stable-url]: https://emmt.github.io/TwoDimensional.jl/stable

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/TwoDimensional.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[github-ci-img]: https://github.com/emmt/TwoDimensional.jl/actions/workflows/CI.yml/badge.svg?branch=main
[github-ci-url]: https://github.com/emmt/TwoDimensional.jl/actions/workflows/CI.yml?query=branch%3Amain

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/TwoDimensional.jl?svg=true
[appveyor-url]: https://ci.appveyor.com/project/emmt/TwoDimensional-jl

[codecov-img]: https://codecov.io/gh/emmt/TwoDimensional.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/emmt/TwoDimensional.jl

[julia-url]: https://julialang.org/
[julia-pkgs-url]: https://pkg.julialang.org/

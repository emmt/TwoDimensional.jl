# TwoDimensional

[![License][license-img]][license-url]
[![Stable][doc-stable-img]][doc-stable-url]
[![Dev][doc-dev-img]][doc-dev-url]
[![Build Status][github-ci-img]][github-ci-url]
[![Build Status][appveyor-img]][appveyor-url]
[![codecov](https://codecov.io/gh/emmt/TwoDimensional.jl/graph/badge.svg?token=m7V2JI1HA6)](https://codecov.io/gh/emmt/TwoDimensional.jl)

`TwoDimensional` is a [Julia][julia-url] package which provides useful types
and methods to define and manipulate 2-dimensional objects (points, rectangles,
circles, polygons, and bounding-boxes) and affine coordinate transforms. This
package also offers methods to build masks from the composition of elementary
shapes.

Other related packages:
- [CoordinateTransformations](https://github.com/FugroRoames/CoordinateTransformations.jl)
  for coordinate transformations;
- [GeometryBasics](https://github.com/JuliaGeometry/GeometryBasics.jl) for
  basic geometric types;
- [Graphics](https://github.com/JuliaGraphics/Graphics.jl) for basic graphical
  objects and methods;

## Usage

```julia
using TwoDimensional
```

gives you types `AffineTransform`, `Point` and `BoundingBox`.

To avoid conflicts with other packages, you may specifically use/import aliases
to these types with suffixes `2D` like `AffineTransform2D`, `Point2D`,
`BoundingBox2D`, etc. For example:

```julia
using TwoDimensional: AffineTransform2D, Point2D, BoundingBox2D
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

[github-ci-img]: https://github.com/emmt/TwoDimensional.jl/actions/workflows/CI.yml/badge.svg?branch=master
[github-ci-url]: https://github.com/emmt/TwoDimensional.jl/actions/workflows/CI.yml?query=branch%3Amaster

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/TwoDimensional.jl?svg=true
[appveyor-url]: https://ci.appveyor.com/project/emmt/TwoDimensional-jl

[julia-url]: https://julialang.org/
[julia-pkgs-url]: https://pkg.julialang.org/

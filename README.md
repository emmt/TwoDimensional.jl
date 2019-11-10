# TwoDimensional

| **Documentation**               | **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][doc-dev-img]][doc-dev-url] | [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

**TwoDimensional** is a [Julia][julia-url] package which provides useful types
and methods to define and manipulate 2-dimensional points, bounding boxes and
affine coordinate transforms.

Other related packages:
- [CoordinateTransformations](https://github.com/FugroRoames/CoordinateTransformations.jl)

## Usage

```julia
using TwoDimensional
```

gives you types `AffineTransform{T}`, `Point{T}` and `BoundingBox{T}` parameterized by the
type `T` of their components (`T` must be floating point for `AffineTransform{T}`).

To avoid conflicts with other packages, you may use/import
`TwoDimensional.Suffixed` which gives you types `AffineTransform2D{T}`,
`Point2D{T}` and `BoundingBox2D{T}` instead, that is with suffix `2D`.


## Installation

`TwoDimensional` is not yet an [official Julia package][julia-pkgs-url] so you
have to clone the repository.  In Julia, hit the `]` key to switch to the
package manager REPL (you should get a `... pkg>` prompt) and type:

```julia
pkg> add https://github.com/emmt/TwoDimensional.jl.git
```

if you use HTTPS, or:

```julia
pkg> add git@github.com:emmt/TwoDimensional.jl.git
```

if you use SSH.

[doc-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[doc-stable-url]: https://emmt.github.io/TwoDimensional.jl/stable

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/TwoDimensional.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[travis-img]: https://travis-ci.org/emmt/TwoDimensional.jl.svg?branch=master
[travis-url]: https://travis-ci.org/emmt/TwoDimensional.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/TwoDimensional.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/TwoDimensional-jl/branch/master

[coveralls-img]: https://coveralls.io/repos/emmt/TwoDimensional.jl/badge.svg?branch=master&service=github
[coveralls-url]: https://coveralls.io/github/emmt/TwoDimensional.jl?branch=master

[codecov-img]: http://codecov.io/github/emmt/TwoDimensional.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/emmt/TwoDimensional.jl?branch=master

[julia-url]: https://julialang.org/
[julia-pkgs-url]: https://pkg.julialang.org/

# Introduction

`TwoDimensional` is a [Julia](https://julialang.org/) package which provides
useful types and methods to define and manipulate 2-dimensional objects
(points, bounding-boxes) and affine coordinate transforms.

Other related packages:
- [CoordinateTransformations](https://github.com/FugroRoames/CoordinateTransformations.jl)
  for coordinate transformations;
- [GeometryBasics](https://github.com/JuliaGeometry/GeometryBasics.jl) for
  basic geometric types;
- [Graphics](https://github.com/JuliaGraphics/Graphics.jl) for basic graphical
  objects and methods;

The source code of `TwoDimensional` is available on
[GitHub](https://github.com/emmt/TwoDimensional.jl).


## Exported Types

```julia
using TwoDimensional
```

gives you types [`AffineTransform{T}`](@ref AffineTransform), [`Point{T}`](@ref
Point) and [`BoundingBox{T}`](@ref BoundingBox) parameterized by the type `T`
of their components (`T` must be floating point for [`AffineTransform{T}`](@ref
AffineTransform)).

To avoid conflicts with other packages, you may use/import
`TwoDimensional.Suffixed` which gives you types [`AffineTransform2D{T}`](@ref
AffineTransform2D), [`Point2D{T}`](@ref Point2D) and [`BoundingBox2D{T}`](@ref
BoundingBox2D) instead, that is with suffix `2D`.

You can also fine tune what you want.  For instance:

```julia
using TwoDimensional: AffineTransform, Point2D
```


## Table of contents

```@contents
Pages = ["install.md", "AffineTransforms.md",
         "Points.md", "BoundingBoxes.md", "reference.md"]
```

## Index

```@index
```

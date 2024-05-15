# Introduction

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

The source code of `TwoDimensional` is available on
[GitHub](https://github.com/emmt/TwoDimensional.jl).


## Exported Types

```julia
using TwoDimensional
```

gives you types [`AffineTransform`](@ref TwoDimensional.AffineTransform),
[`Point`](@ref TwoDimensional.Point), [`Rectangle`](@ref
TwoDimensional.Rectangle), [`Circle`](@ref TwoDimensional.Circle),
[`Polygon`](@ref TwoDimensional.Polygon), and [`BoundingBox`](@ref
TwoDimensional.BoundingBox).

To avoid conflicts with other packages, you may specifically use/import aliases
to these types with suffixes `2D` like [`AffineTransform2D`](@ref
TwoDimensional.AffineTransform2D), [`Point2D`](@ref TwoDimensional.Point2D),
[`BoundingBox2D`](@ref TwoDimensional.BoundingBox2D), etc. For example:

```julia
using TwoDimensional: AffineTransform2D, Point2D, Polygon2D, BoundingBox2D
```

## Table of contents

```@contents
Pages = [
    "install.md",
    "points.md",
    "boxes.md",
    "xforms.md",
    "reference.md",
]
```

## Index

```@index
```

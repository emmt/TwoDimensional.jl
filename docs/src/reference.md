# Reference

The following summarizes the documentation of types and methods provided by the
`TwoDimensional` package.  This information is also available from the REPL by
typing `?` followed by the name of a method or a type.


## Affine 2D Coordinate Transforms

```@docs
AffineTransform
scale
rotate
translate
jacobian
intercept
compose
```

## Points

```@docs
AbstractPoint
Point
WeightedPoint
distance
```

## Bounding-Boxes

```@docs
BoundingBox
area
center
interior
exterior
```


## Methods

```@docs
round(::Union{Point,BoundingBox})
floor(::Point)
ceil(::Point)
```

## Aliases

```@docs
AbstractPoint2D
Point2D
WeightedPoint2D
BoundingBox2D
AffineTransform2D
```

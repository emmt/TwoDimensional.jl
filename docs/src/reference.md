# Reference

The following summarizes the documentation of types and methods provided by the
`TwoDimensional` package.  This information is also available from the REPL by
typing `?` followed by the name of a method or a type.


## Points

```@docs
TwoDimensional.AbstractPoint
TwoDimensional.AbstractPoint2D
TwoDimensional.Point
TwoDimensional.Point2D
TwoDimensional.distance
```

## Bounding-Boxes

```@docs
TwoDimensional.BoundingBox
TwoDimensional.BoundingBox2D
TwoDimensional.interior
TwoDimensional.exterior
TwoDimensional.grow
TwoDimensional.shrink
```

## Rectangles

```@docs
TwoDimensional.Rectangle
TwoDimensional.Rectangle2D
TwoDimensional.area
```

## Coordinate Type

```@docs
TwoDimensional.coord_type
TwoDimensional.convert_coord_type
TwoDimensional.promote_coord_type
```

## Affine 2D Coordinate Transforms

```@docs
TwoDimensional.AffineTransform
TwoDimensional.AffineTransform2D
TwoDimensional.compose
TwoDimensional.intercept
TwoDimensional.jacobian
TwoDimensional.rotate
TwoDimensional.factors_type
TwoDimensional.offsets_type
```

## Masks

```@docs
TwoDimensional.MaskElement
TwoDimensional.MaskElement2D
TwoDimensional.Overlap
TwoDimensional.apply_mask
TwoDimensional.apply_mask!
TwoDimensional.circular_aperture
TwoDimensional.circular_obscuration
TwoDimensional.forge_mask
TwoDimensional.forge_mask!
TwoDimensional.grid_step
TwoDimensional.interpolate
TwoDimensional.is_opaque
TwoDimensional.is_transparent
TwoDimensional.polygonal_aperture
TwoDimensional.polygonal_obscuration
TwoDimensional.rectangular_aperture
TwoDimensional.rectangular_obscuration
```

## Other Public Methods

```@docs
TwoDimensional.center
TwoDimensional.scale
TwoDimensional.translate
Base.floor(::TwoDimensional.Point)
Base.ceil(::TwoDimensional.Point)
Base.round(::TwoDimensional.Point)
Base.vec(::TwoDimensional.Polygon)
```

## Internal Methods

```@docs
TwoDimensional.apply
TwoDimensional.geometric_properties
TwoDimensional.is_convex
TwoDimensional.is_nothing
TwoDimensional.is_something
TwoDimensional.parts
TwoDimensional.shape
TwoDimensional.vertices
```

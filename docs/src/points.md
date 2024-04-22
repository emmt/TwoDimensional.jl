# Points

Any object whose type is derived from [`AbstractPoint{T}`](@ref AbstractPoint)
has 2-D coordinates: its abscissa and ordinate respectively named `x` and `y`.
The parameter `T` is the type of the coordinates and can be retrieved by the
`eltype` method.


## Aliases

Call:

```julia
using TwoDimensional: AbstractPoint2D, Point2D
```

instead of:

```julia
using TwoDimensional
```

to have `Point2D`, `AbstractPoint2D` provided as respective aliases to
`TwoDimensional.Point` and `TwoDimensional.AbstractPoint`.


## Construction

The most simple concrete type is [`Point{T}`](@ref Point) constructed by:

```julia
Point(x,y)
```

where `(x,y)` are the coordinates of the point.


Coordinates and weights can also be specified by keywords:

```julia
Point(x=xval, y=yval)
```

There are no default values for keywords `x` and `y` so they must both be
specified.


## Fields

The fields of a `Point`, say `pnt`, can be retrieved in different ways:

```julia
x = pnt.x
y = pnt.y
```

or:

```julia
x = pnt[1]
y = pnt[2]
```

or:

```julia
x, y = pnt
```


## Conversion

Simple points can be constructed from a 2-tuple of coordinates or from an
instance of 2-dimensional `CartesianIndex`:

```julia
v = (x,y)
I = CartesianIndex(x,y)
Point(v)    # yields Point(x,y)
Point(I)    # yields Point(x,y)
```

and reciprocally, assuming `P = Point(x,y)`:

```julia
Tuple(P)          # yields the 2-tuple (P.x, P.y)
CartesianIndex(P) # yields CartesianIndex(P.x, P.y)
```

Coordinate type conversion, say to type `T`, for a point `P = Point(x,y)` is
done by:

```julia
Point{T}(P)
convert(Point{T}, P)
T.(P)
```

The latter form involves broadcasting rules and may be a bit slower.


## Operations on Points

The addition (resp. subtraction) of two points adds (resp. subtracts) their
coordinates:

```julia
Point(x1,y1) + Point(x2,y2)   # yields Point(x1+x2,y1+y2)
Point(x1,y1) - Point(x2,y2)   # yields Point(x1-x2,y1-y2)
```

Unary minus of a point negates its coordinates:

```julia
-Point(x,y)   # yields Point(-x,-y)
```

A point may be multiplied or divided by a scalar, say `α`, to scale its
coordinates:

```julia
α*Point(x,y)  # yields Point(α*x,α*y)
Point(x,y)*α  # yields Point(α*x,α*y)
α\Point(x,y)  # yields Point((1/α)*x,(1/α)*y)
Point(x,y)/α  # yields Point((1/α)*x,(1/α)*y)
```

Taking the hypothenuse or the arctangent of a point `P = Point(x,y)` yield its
distance to the origin `O = Point(0,0)` and the angle between `OP` and
the abscissa axis:

```julia
hypot(Point(x,y))  # yields hypot(x, y)
atan(Point(x,y))   # yields atan(y, x)
```

The distance between two points is given by the [`distance`](@ref) method:

```julia
distance(Point(x1,y1),Point(x2,y2))  # yields hypot(x1-x2,y1-y2)
```

The nearest point to an instance `obj` of [`Point`](@ref) is given by:

```julia
round([T,] obj, [r::RoundingMode])
```

which rounds the coordinates of `obj` to the nearest integral value. Optional
argument `T` is to specify the type of the result or the type of the
coordinates of the result. Optional argument `r` is the rounding mode, the
default is the same as `round` for a scalar value.

Similarly, `floor([T,],P)` and `ceil([T,],P)` yield the point with integer
coordinates immediately (inclusively) before and after [`Point`](@ref) `P`.


A point can be clamped within the limits of a bounding-box:

```julia
clamp(Point(-1.1, 6.3), BoundingBox(1:4,1:5)) # yields Point(1.0,5.0)
```

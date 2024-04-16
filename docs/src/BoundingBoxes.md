# Bounding-Boxes

2-D bounding-boxes are built by:

```julia
BoundingBox(xmin,xmax,ymin,ymax)
```

to represent a 2D rectangular box whose sides are aligned with the coordinate
axes and containing points of coordinates `(x,y)` such that `xmin ≤ x ≤ xmax`
and `ymin ≤ y ≤ ymax`.

The type of the bounds, say `T`, can be explicitly specified:

```julia
BoundingBox{T}(xmin,xmax,ymin,ymax)
```

If unspecified, it is assumed to be he types of the bounds promoted to a common
type.  The type of the bounds and can be retrieved by the `eltype` method.


## Aliases

Call:

```julia
using TwoDimensional.Suffixed
```

instead of:

```julia
using TwoDimensional
```

to have [`BoundingBox2D`](@ref BoundingBox2D) provided as an alias to
[`TwoDimensional.BoundingBox`](@ref BoundingBox).


## Construction

The coordinates of a bounding-box can be specified by keywords:

```julia
BoundingBox(xmin=x0, ymin=y0, xmax=x1, ymax=y1)
```

There are no default values for keywords `xmin`, `xmax`, `ymin` and `ymax` so
all must be specified.

A bounding-box can be constructed from a 4-tuple of coordinates and conversely:

```julia
BoundingBox((x0,x1,y0,y1))      # yields BoundingBox(x0,x1,y0,y1)
Tuple(BoundingBox(x0,x1,y0,y1)) # yields (x0,x1,y0,y1)
```

A bounding-box can be constructed from its first and last points (*i.e.* at the
lower-left and upper right opposite corners) specified as instances of `Point`,
of `CartesianIndex{2}` or of `Tuple{Real,Real}`:

```julia
BoundingBox(Point(x0,y0), Point(x1,y1))
BoundingBox(CartesianIndex(x0,y0), CartesianIndex(x1,y1))
BoundingBox((x0,y0), (x1,y1))
```

which all yield the same result:

```julia
BoundingBox(x0,x1,y0,y1)
```

Conversely, methods `first(B)` and `last(B)` respectively yield the lower left
and upper right corners of the bounding-box `B` (as a [`Point`](@ref)
instance):

```julia
first(BoundingBox(x0,x1,y0,y1)) # yields Point(x0,y0)
last(BoundingBox(x0,x1,y0,y1))  # yields Point(x1,y1)
```

Integer-valued unit-ranges can be specified to define a bounding-box.  For
example:

```julia
BoundingBox(x0:x1, y0:y1)    # 2 unit-range
BoundingBox((x0:x1, y0:y1))  # a 2-tuple of unit range
```

This makes possible writing:

```julia
BoundingBox(axes(A))
```

to get the bounding-box corresponding to all indices of array `A`.  Conversely:

```julia
axes(BoundingBox(x0,x1,y0,y1))
```

yields the axes of a bounding-box with integer coordinates, that is
`(x0:x1,y0:y1)`.  To get the `k`-th axis of a bounding-box `B`, call
`axes(B,k)`.

To loop over the Cartesian indices defined by a bounding-box `B` with integer
coordinates, you can just write:

```julia
for I in CartesianIndices(B)
   ...
end
```

A bounding-box may also be constructed by applying a predicate function to the
elements of a 2-dimensional array:

```julia
BoundingBox(f, A)
```

yields the bounding-box of all integer coordinates `(x,y)` such that
`f(A[x,y])` yields `true`.  If the elements of `A` are booleans (of type
`Bool`), then `BoundingBox(A)` is equivalent to `BoundingBox(identity,A)`.


## Fields

The fields of a `BoundingBox`, say `box`, can be retrieved in different ways:

```julia
xmin = box.xmin
xmax = box.xmax
ymin = box.ymin
ymax = box.ymax
```

or:

```julia
xmin = box[1]
xmax = box[2]
ymin = box[3]
ymax = box[4]
```

or:

```julia
xmin, xmax, ymin, ymax = box
```

## Conversion

Coordinate type conversion, say to type `T`, is done by:

```julia
B = BoundingBox(x0,x1,y0,y1)
BoundingBox{T}(B)
convert(BoundingBox{T}, B)
T.(B)
```

The latter form involves broadcasting rules and may be a bit slower.


## Union and Intersection of Bounding-Boxes

The union of bounding-boxes `b1`, `b2`, ... is given by one of:

```julia
B1 ∪ B2 ∪ ...
union(B1, B2, ...)
```

which both yield the largest bounding-box contained into the bounding-boxes
`B1`, `B2`, ...

The intersection of bounding-boxes `B1`, `B2`, ... is given by one of:

```julia
B1 ∩ B2 ∩ ...
intersect(B1, B2, ...)
```

which both yield the smallest bounding-box containing the bounding-boxes `B1`,
`B2`, ...

The maximal or minimal bounding-box with coordinates of type `T` that can be
constructed are respectively given by `typemax(BoundingBox{T})` and
`typemin(BoundingBox{T})`. These can be useful to initiate a shrinking or a
growing bounding-box. The call:

```julia
BoundingBox{T}(nothing)
```

yields the same result as `typemin(BoundingBox{T})`.

The method `isempty(B)` yields whether a bounding-box `B` is empty or not.


## Interior, Exterior, Nearest, etc.

Given the bounding-box `B`, [`interior(B)`](@ref interior) and
[`exterior(B)`](@ref exterior) respectively yield the largest interior and
smallest exterior bounding-boxes with integer bounds.

[`round(B)`](@ref round), [`round(T,B)`](@ref round), or [`round(T,B,r)`](@ref
round) yield a bounding-box whose limits are those of the bounding-box `B`
rounded to the nearest integral values of type `T` with rounding-mode `r`.
Default type `T` is `eltype(B)` and default rounding-mode `r` is
`RoundingNearest`.

[`center(B)`](@ref center) yields the `Point` whose coordinates are the
geometrical center of the bounding-box `B`.

[`area(B)`](@ref area) yields the area of a bounding-box `B`.


## Arithmetic and Basic Methods

Adding or subtracting a scalar `δ` to a bounding-box `B` adds or removes a
margin `δ` to the bounding-box `B`:

```julia
BoundingBox((x0,y0),(x1,y1)) + δ -> BoundingBox((x0-δ,y0-δ),(x1+δ,y1+δ))
BoundingBox((x0,y0),(x1,y1)) - δ -> BoundingBox((x0+δ,y0+δ),(x1-δ,y1-δ))
```

Adding or subtracting a point `P` to a bounding-box `B` shifts the limits of
the bounding-box `B`:

```julia
BoundingBox((x0,y0),(x1,y1)) + Point(x,y) -> BoundingBox((x0+x,y0+y),(x1+x,y1+y))
```

A bounding-box `B` can be negated:
```julia
-BoundingBox((x0,y0),(x1,y1)) -> BoundingBox((-x1,-y1),(-x0,-y0))
```

`eltype(B)` yields the type of the coordinates of a bounding-box `B`.

Basic methods `size(B[,k])` and `axes(B[,k])` can be applied to an
**integer-valued** bounding-box `B`.  These two methods are type-stable:
`size(B)` yields a 2-tuple of `Int`, `size(B,k)` yields an `Int`, `axes(B)`
yields a 2-tuple of `UnitRange{Int}` and `axes(B,k)` yields a `UnitRange{Int}`.

Method `in` and operator `∈`, obtained by `\in`-tab, yield whether a point
`pnt` is inside a bounding-box `box`:

```julia
pnt ∈ box
```

is a shortcut for:

```julia
(box.xmin ≤ pnt.x ≤ box.xmax) & (box.ymin ≤ pnt.y ≤ box.ymax)
```

Method `issubset` or operator `⊆`, obtained by `\subseteq`-tab, yield whether a
bounding-box, say `A`, is inside another one, say `B`. That is `issubset(A, B)`
and `A ⊆ B` are shortcuts for:

```julia
(isempty(A) | ((A.xmin ≥ B.xmin) & (A.xmax ≤ B.xmax) &
               (A.ymin ≥ B.ymin) & (A.ymax ≤ B.ymax)))
```

# User visible changes in the `TwoDimensional` package

## Unreleased

### Added

- Use `@public` macro from [`TypeUtils`](https://github.com/emmt/TypeUtils.jl) package for
  all documented non-exported symbols.

- Extend and export `TypeUtils.get_precision` and `TypeUtils.adapt_precision` for geometric
  objects.

- Array with offsets can be automatically built by `forge_mask` when the `OffsetArrays`
  package is loaded.

- `forge_mask` can create output array given its element type and dimensions.

### Deprecated

- Computing the cross-product of the coordinates of points `a` and `b` by `a*b` is
  deprecated in favor of `cross(a,b)`.

### Fixed

- `coord_type` calls `TypeUtils.to_same_concrete_type`, not `promote_type`, to ensure that
  the result is a common concrete type.


## Version 0.5.0

This version of `TwoDimensional` introduces new geometric object types (`Rectangle`,
`Polygon`, and `Circle`), the possibility for coordinates to have units, more consistent
arithmetic rules, and methods to build complex masks by composing elementary geometric
objects.

### General changes

- **Coordinate type:** Coordinates of geometric objects may be of any type implementing
  basic arithmetic operations and may thus have units.

  - Method `coord_type` yields the type of the coordinates of a geometric object/type or,
    with a list of geometric objects/types, yields the promoted coordinates type of these
    objects/types.

  - Method `convert_coord_type` converts the coordinate type of a geometric object/type.

  - Method `promote_coord_type` converts all its arguments to a common coordinate type and
    return them as a tuple.

  - Methods `convert_bare_type`, `convert_real_type`, and `convert_floating_point_type` of
    the `TypeUtils` package have been extended to convert the numerical type of the
    coordinates of graphical objects.

  - Method `float` may be called to convert the coordinate type of a geometric object/type
    to a floating-point type.

- For an object `obj` having homogeneous components (points are like 2-tuple of
  coordinates, rectangles and boxes are like 2-tuple of points, polygons are like vectors
  of points), `length(obj)` and `eltype(obj)` yield the number and the type of these
  components.

- New geometric object types `Rectangle`, `Circle`, and `Polygon`. Rectangles are similar
  to bounding-boxes except that they can never be empty (a rectangle contains at least a
  single point) and that they are converted into polygons when transformed by an affine
  coordinate transform.

- Non-exported method `TwoDimensional.apply(f, obj)` applies the function `f` to each
  component of the geometric object `obj` and rebuilds an object of the same kind with the
  result. This method is only applicable to geometric objects having homogeneous parts:
  points, rectangles, bounding-boxes, and polygons. For bounding-boxes, keyword `swap`
  specifies whether to swap the bounds of the box, `box.start` and `box.stop`.

- `zero(obj)` and `one(obj)` yield the additive and multiplicative identities for the type
  of the geometric object `obj`. That is such that `zero(obj) + obj == obj + zero(obj) ==
  obj` and `one(obj)*obj == obj*one(obj) == obj` hold. These rules implies that `one(obj)`
  is the scalar `one(coord_type(obj))`; while, `zero(obj)` is a point with the same
  coordinate type as `obj` and with coordinates equal to zero.

- `Base.convert` method has been specialized to implement most allowed conversion from
   point-like objects (2-tuple of coordinates, abstract points, and Cartesian indices) to
   points, from bounding-box-like objects (2-tuple of points, of Cartesian indices, of
   unit-ranges, of 2-tuple of coordinates, etc.) to bounding-boxes, and similarly for
   rectangles and other geometric object types. Hence, the `TypeUtils.as` method no longer
   needs to be extended to perform those conversions.

### Points and bounding-boxes

- Use `A ⊆ B` to check whether bounding-box `A` is contained in bounding-box `B` (the
  syntax `A ∈ B` is no longer implemented for that) and `P ∈ B` to check whether point `P`
  is inside bounding-box `B`.

- Points can be defined by their polar coordinates with `pnt = Point(r=..., θ=...)`.
  `hypot(pnt)` and `abs(pnt)` yield the distance `r` of `pnt` to the origin, `atan(pnt)`
  yields its polar angle `θ`, and `abs2(pnt)` yields its squared distance `r^2` to the
  origin.

- Points are considered as 2-dimensional vectors for some operations: `norm(a)` yields the
  Euclidean norm of the point `a`; `dot(a, b)` or `a ⋅ b` compute the dot product (a.k.a.
  scalar product or inner product) of the two points `a` and `b`; while `cross(a, b)` or
  `a*b` compute their cross product (a.k.a. vector product or directed area product).

- `BoundingBox(obj)` yields the bounding-box of the geometric object `obj`.
  `BoundingBox(f, A)` yields the bounding-box of the entries of the 2-dimensional array
  `A` such that `f(A[i,j])` is true. If `A` is an array of Booleans, `f` is assumed to be
  the identity if not specified.

- Addition and subtraction of bounding-boxes, say `C = A ± B`, yields the bounding-box `C`
  for all points `c = a ± b` whatever `a ∈ A` and `b ∈ B`.

- `BoundingBox{T}()` yields an empty bounding-box for coordinate type `T`.

- Taking the center of an empty bounding-box throws an error.

### Affine transforms

- Unary minus is implemented for affine transforms.

- The type of the result of applying an affine transform follows Julia type promotion
  rules.

- Since the factors and the offsets in an affine transform `A` may have different units,
  `eltype(A)` is no longer applicable; `bare_type(A)`, `real_type(A)`, or
  `floating_point_type(A)` may be used to retrieve the floating-point type of the
  coefficients of `A` (these methods require `using TypeUtils`). Non-exported methods
  `TwoDimensional.factors_type(A)` and `TwoDimensional.offsets_type(A)` yield the types of
  the factors ( the coefficients `A.xx`, `A,xy`, `A.yx`, and `A.yy`) and of the offsets
  (the coefficients `A.x` and `A.y`) of `A`. For the same reasons, `T(A)` and `T.(A)` with
  `T` a floating-point type is no longer supported; call any of `convert_bare_type(T,A)`,
  `convert_real_type(T,A)`, or `convert_floating_point_type(T,A)` to convert the
  floating-point type of the coefficients of `A` (these methods require `using
  TypeUtils`).

- Method `TwoDimensional.compose` is no longer exported. Explicitly use/import it or use
  `*`, or `∘` (`\circ<tab>`) to compose affine transforms.

- Due to the definition of *"intercept"* in mathematics, `TwoDimensional.intercept(A)`,
  with `A` an affine transform, has been replaced by `TwoDimensional.ldiv(A,b=(0,0))`
  which yields `c` such that `A*c = b` where `b` is `(0,0)` by default. If `b` is
  specified as a `Point`, `c` is returned as a `Point` as well. With the same (slight)
  abuse of notation, than `A*c -> A(c)`, the `\\` operator is overloaded so that `A\\b ->
  solve(A,b)`.

- Non-exported methods `TwoDimensional.leftdivide` and `TwoDimensional.rightdivide` have
  been fixed and renamed `TwoDimensional.ldiv` and `TwoDimensional.rdiv`.

### Masks

It is now possible to define and apply arbitrary masks by combining elementary
*obscurations* (opaque elementary geometric objects) and *apertures* (transparent
elementary geometric objects). The `Mask` constructor takes any number of elementary mask
objects and combine them into a composite mask. Methods `forge_mask` and `forge_mask!`
compute the transmission of a mask at each cell of an array using antialiasing rules and
multi-threaded computations.

### Things no longer supported

- **Weighted points**, of non-exported type `TwoDimensional.WeightedPoint`, have been
  removed due to coordinates possibly having units. They may come back but with 2 type
  parameters: one for the weight, the other for the coordinates.

- Four arguments constructor `BoundingBox(xmin,xmax,ymin,xmax)` is no longer supported as
  it was a source of confusion and of many errors. Use
  `BoundingBox((xmin,ymin),(xmax,xmax))` or `BoundingBox(xmin:xmax,ymin:xmax)` (for
  integer coordinates), or use keywords `BoundingBox(xmin=..., xmax=..., ymin=...,
  xmax=...)`.

- Operator `⋅` (`\cdot<tab>`) is reserved for the dot product. Applying an affine
  transform `A` to an object `obj` is done by `A(obj)` or `A*obj`. Composing affine
  transforms `A` and `B` is done by `A*B` or `A∘B` (with `∘` given by `\circ<tab>` in the
  REPL)

- Julia versions older than 1.0 are no longer supported.

- Sub-module `TwoDimensional.Suffixed` has been suppressed, aliases must specifically be
  used/imported. For example with: `using TwoDimensional: AffineTransform2D, Point2D,
  BoundingBox2D`.

- `AffineTransforms` is no longer a sub-module of `TwoDimensional`. If only 2-dimensional
  affine transforms are needed, do `using TwoDimensional: AffineTransform` or `using
  TwoDimensional: AffineTransform2D`. Most related operations are available by symbolic
  operators; otherwise, they can also be used/imported. For example: `using
  TwoDimensional: AffineTransform2D, rotate`.

- Deprecated `BoundingBox(A::AbstractArray) -> BoundingBox(axes(A))` has been removed.

- Due to the definition of *"Jacobian"* in mathematics, the `jacobian` method has been
  removed, call `abs(det(A))`.

## Version 0.4.1

- Tests run on Julia 1.11.

- Fix documentation.

- Add missing specialization for `zero`, `one`, and `oneunit` for `Point` instances.

- Non-exported aliases `TwoDimensional.PointLike` and `TwoDimensional.BoundingBoxLike` for
  unions of types that can be used to specify a `Point` or a `BoundingBox`.

- Out of bounds indices in point and bounding-box instances throw `BoundError` (was
  `ErrorException`).

- `size(box)` for bounding-box `box` may yield non-zero dimension even though `box` is
  empty (i.e. same behavior as with arrays).

## Version 0.4.0

- Extend `map` for `Point` and `BoundingBox` objects.

- `round` for a `Point` or a `BoundingBox` has an optional `r::RoundingMode` argument.

- Extend `ceil`, and `floor` methods for `BoundingBox` objects.

- Unify extends for `Point` and `BoundingBox` of the `round`, `ceil`, and `floor` methods
  and make their rules more consistent with the corresponding method for scalar values.

## Version 0.3.1

- Julia 0.7 is supported.

## Version 0.3.0

- `A ⊆ B` and `issubset(A, B)` yield whether bounding-box `A` is inside the bounding-box
  `B`. For bounding-boxes `A` and `B`, `A ∈ B` has been deprecated in favor of `A ⊆ B`.

## Version 0.2.1

- Instances of `Point`, `WeightedPoint` or `BoundingBox` are iterable. Hence `x, y = pnt`
  can be used to extract the coordinates of `pnt` an instance of `Point`.

- `A ∈ B` yields whether `A`, a point or a bounding-box, is inside the bounding-box `B`.

- Unary plus can be applied to a point or to a bounding-box.

## Version 0.2.0

### New functionalities and improvements

- The bounding-box algorithm has been largely improved. Except for very small matrices, it
  is much faster than the former naive implementation.

- A bounding-box can be negated and can be built from two 2-tuples of coordinates.

- A point can be clamped within the limits of a bounding-box.

- `zero` and `one` methods have been extended for type `Point`.

### Changes of behavior

- `BoundingBox(A)` is no longer equivalent to `BoundingBox(axes(A))` for a 2-dimensional
  array `A`. `BoundingBox(A)` is equivalent to `BoundingBox(identity,A)` if the elements
  of `A` are Booleans (of type `Bool`).

## Version 0.1.0

First official version of `TwoDimensional`.

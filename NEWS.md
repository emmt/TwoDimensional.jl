# User visible changes in the `TwoDimensional` package

## Version 0.5.0

### General changes

- Points, bounding-boxes, and affine transforms coefficients may have units.

- Julia versions older than 1.0 are no longer supported.

- Sub-module `TwoDimensional.Suffixed` has been suppressed, aliases must
  specifically used/imported. For example with: `using TwoDimensional:
  AffineTransform2D, Point2D, BoundingBox2D`.

- `AffineTransforms` is no longer a sub-module of `TwoDimensional`. If only
  2-dimensional affine transforms are needed, do `using TwoDimensional:
  AffineTransform` or `using TwoDimensional: AffineTransform2D`. Most related
  operations are available by symbolic operators; otherwise, they can also be
  used/imported. For example: `using TwoDimensional: AffineTransform2D,
  rotate`.

- Exported method `coord_type` yields the type of the coordinates for a
  geometric object or for its type. Method `promote_coord_type` yields the
  promoted coordinates type from a list of geometric objects/types.

### Points and bounding-boxes

- `map` and broadcasting rules have been extended to be more consistent for
  points and bounding-boxes:
  - for a point: `f.(pnt) = map(f,pnt) = Point(map(f,Tuple(pnt)))`
  - for a bounding-box: `f.(box) = map(f,box) = BoundingBox(map(f,Tuple(box)))`
    and `map` takes a `swap` keyword (`false` by default) to specify whether to
    swap the inferior and superior bounds along each dimension which is needed
    when, for example, multiplying a bounding-box by a negative factor.

- `zero(obj)` and `one(one)` yield the additive and multiplicative identities
  for the type of object `obj`, a point or a bounding-box. That is such that
  `zero(obj) + obj == obj + zero(obj) == obj` and `one(obj)*obj == obj*one(obj)
  == obj` hold. These rules implies that `one(obj)` is the scalar
  `one(eltype(obj))`; while, `zero(obj)` is an object of the same type as `obj`
  with all values set to zero.

- Use `A ⊆ B` to check whether bounding-box `A` is contained in bounding-box
  `B` (the syntax `A ∈ B` is no longer implemented for that) and `P ∈ B` to
  check whether point `P` is inside bounding-box `B`.

- Points can be defined by their polar coordinates and `hypot`, `abs`, and
  `norm` yield the distance of a point distance to the origin, while `atan`
  yields its polar angle and `abs2` yields its squared distance to the origin.
  `TwoDimensional.inner` and `TwoDimensional.outer` compute the inner and outer
  products of two points.

- `BoundingBox(obj)` yields the bounding-box of the geometric object `obj`.
  `BoundingBox(f,arr)` yields the bounding-box of the entries of the
  2-dimensional array `A` such that `f(A[i,j])` is true. If `A` is an array of
  Booleans, `f` is assumed to the identity if not specified.

- Weighted points, of non-exported type `TwoDimensional.WeightedPoint`, have
  been removed due to coordinates possibly having units. They may come back but
  with 2 type parameters: one for the weight, the other for the coordinates.

### Bounding-boxes

- Addition and subtraction of bounding-boxes, say `C = A ± B`, yields the
  bounding-box `C` for all points `c = a ± b` whatever `a ∈ A` and `b ∈ B`.

- Deprecated `BoundingBox(A::AbstractArray) -> BoundingBox(axes(A))` has been
  removed.

- Method `TypeUtils.as` is extended to convert points and bounding-boxes
  to/from tuples.

- Taking the center of an empty bounding-box throws an error.

### Affine transforms

- Unary minus is implemented for affine transforms.

- The type of the result of applying an affine transform follows Julia
  type promotion rules.

- Since the factors and the offsets in an affine transform `A` may have
  different units, `eltype(A)` is no longer applicable; `bare_type(A)`,
  `real_type(A)`, or `floating_point_type(A)` may be used to retrieve the
  floating-point type of the coefficients of `A` (these methods require `using
  Unitless`). Non-exported methods `TwoDimensional.factors_type(A)` and
  `TwoDimensional.offsets_type(A)` yield the types of the factors ( the
  coefficients `A.xx`, `A,xy`, `A.yx`, and `A.yy`) and of the offsets (the
  coefficients `A.x` and `A.y`) of `A`. For the same reasons, `T(A)` and
  `T.(A)` with `T` a floating-point type is no longer supported; call any of
  `convert_bare_type(T,A)`, `convert_real_type(T,A)`, or
  `convert_floating_point_type(T,A)` to convert the floating-point type of the
  coefficients of `A` (these methods require `using Unitless`).

- Method `TwoDimensional.compose` is no longer exported. Use `*`, `⋅`
  (`\cdot<tab>`), or `∘` (`\circ<tab>`) to compose affine transforms.

## Version 0.4.1

- Tests run on Julia 1.11.

- Fix documentation.

- Add missing specialization for `zero`, `one`, and `oneunit` for `Point`
  instances.

- Non-exported aliases `TwoDimensional.PointLike` and
  `TwoDimensional.BoundingBoxLike` for unions of types that can be used to
  specify a `Point` or a `BoundingBox`.

- Out of bounds indices in point and bounding-box instances throw `BoundError`
  (was `ErrorException`).

- `size(box)` for bounding-box `box` may yield non-zero dimension even though
  `box` is empty (i.e. same behavior as with arrays).

## Version 0.4.0

- Extend `map` for `Point` and `BoundingBox` objects.

- `round` for a `Point` or a `BoundingBox` has an optional `r::RoundingMode`
  argument.

- Extend `ceil`, and `floor` methods for `BoundingBox` objects.

- Unify extends for `Point` and `BoundingBox` of the `round`, `ceil`, and
  `floor` methods and make their rules more consistent with the corresponding
  method for scalar values.

## Version 0.3.1

- Julia 0.7 is supported.

## Version 0.3.0

- `A ⊆ B` and `issubset(A, B)` yield whether bounding-box `A` is inside the
  bounding-box `B`. For bounding-boxes `A` and `B`, `A ∈ B` has been deprecated
  in favor of `A ⊆ B`.

## Version 0.2.1

- Instances of `Point`, `WeightedPoint` or `BoundingBox` are iterable. Hence
  `x, y = pnt` can be used to extract the coordinates of `pnt` an instance of
  `Point`.

- `A ∈ B` yields whether `A`, a point or a bounding-box, is inside the
  bounding-box `B`.

- Unary plus can be applied to a point or to a bounding-box.

## Version 0.2.0

### New functionalities and improvements

- The bounding-box algorithm has been largely improved. Except for very small
  matrices, it is much faster than the former naive implementation.

- A bounding-box can be negated and can be built from two 2-tuples of
  coordinates.

- A point can be clamped within the limits of a bounding-box.

- `zero` and `one` methods have been extended for type `Point`.

### Changes of behavior

- `BoundingBox(A)` is no longer equivalent to `BoundingBox(axes(A))` for a
  2-dimensional array `A`. `BoundingBox(A)` is equivalent to
  `BoundingBox(identity,A)` if the elements of `A` are Booleans (of type
  `Bool`).

## Version 0.1.0

First official version of TwoDimensional.

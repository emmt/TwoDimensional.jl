# User visible changes in TwoDimensional package

## Version 0.5.0

- Points, bounding-boxes, and affine transforms coefficients may have units.


- `map` and broadcasting rules have been extended to be more consistent for
  points and bounding-boxes:
  - for a point: `f.(pnt) = map(f,pnt) = Point(map(f,Tuple(pnt)))`
  - for a bounding-box: `f.(box) = map(f,box) = BoundingBox(map(f,Tuple(box)))`
    and `map` takes a `swap` keyword (`false` by default) to specify whether to
    swap the inferior and superior bounds along each dimension which is needed
    when, for example, multiplying a bounding-box by a negative factor.
- Unary minus is implemented for affine transforms.

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
  `convert_floating_point_type(T,A)` may be used to convert the floating-point
  type of the coefficients of `A` (these methods require `using Unitless`).

- The type of the result of applying an affine transform follows Julia
  type promotion rules.

- Deprecated `BoundingBox(A::AbstractArray) -> BoundingBox(axes(A))` has been
  removed.

- Method `TwoDimensional.compose` is no longer exported. Use `*`, `⋅`
  (`\cdot<tab>`), or `∘` (`\circ<tab>`) to compose affine transforms.

- Method `TypeUtils.as` is extended to convert points and bounding-boxes
  to/from tuples.

- Julia versions older than 1.0 are no longer supported.

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

# User visible changes in TwoDimensional package

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

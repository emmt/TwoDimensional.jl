# User visible changes in TwoDimensional package

## Version 0.3.2

- `round` has a `r::RoundingMode` optional argument.

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

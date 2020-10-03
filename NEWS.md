# User visible changes in TwoDimensional package

## Version master

- `A ∈ B` yields whether `A`, a point or a bounding-box, is inside the
  bounding-box `B`.

- Unary plus can be applied to a point or to a bounding-box.


## Version 0.2.0

### New functionalities and improvements

- The bounding-box algorithm has been largely improved.  Except for very small
  matrices, it is much faster than the former naive implementation.

- A bounding-box can be negated and can be built from two 2-tuples of
  coordinates.

- A point can be clamped whithin the limits of a bounding-box.

- `zero` and `one` methods have been extended for type `Point`.


### Changes of behavior

- `BoundingBox(A)` is no longer equivalent to `BoundingBox(axes(A))` for a
  2-dimensional array `A`.  `BoundingBox(A)` is equivalent to
  `BoundingBox(identity,A)` if the elements of `A` are booleans (of type
  `Bool`).


## Version 0.1.0

First official version of TwoDimensional.

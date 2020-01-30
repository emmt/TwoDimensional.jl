# Changes in TwoDimensional package

## Version 0.2.0

### Changes of behavior

- `BoundingBox(A)` is no longer equivalent to `BoundingBox(axes(A))` for a
  2-dimensional array `A`.  `BoundingBox(A)` is equivalent to
  `BoundingBox(identity,A)` if the elements of `A` are booleans (of type
  `Bool`).


## Version 0.1.0

First official version of TwoDimensional.

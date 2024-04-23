- Use `IntervalSets` for bounding-boxes.

- Add masks.

- Implement `Ellipse` and/or `TransformedCircle` to be able to apply an affine
  transform to circles (and ellipses).

- Optimization: Add a special argument in `Rectangle` constructor to skip
  re-ordering coordinates.

- Unify a bit more the API for rectangles and bounding-boxes.

- In current code, points can only be `INSIDE` or `OUTSIDE` a polygon. Fix this
  so that a point may `OVERLAP` to boundaries of a polygon. Also check that the
  algorithm for a rectangular cell overlapping a polygon is correct for a cell
  having zero width or height.

- Implement [Welzl's
  algorithm](https://en.wikipedia.org/wiki/Smallest-circle_problem) for
  determining the diameter of a polygon.

- Extending `convert` for other purpose of type conversion is potentially
  dangerous as this method is called by type constructors to set fields:
  ```julia
  Base.convert(::Type{CartesianIndex}, P::Point{<:Integer}) = CartesianIndex(P)
  Base.convert(::Type{CartesianIndex{2}}, P::Point{<:Integer}) = CartesianIndex(P)
  Base.convert(::Type{T}, P::Point) where {T<:Tuple} = convert(T, Tuple(P))

  Base.convert(::Type{T}, B::BoundingBox) where {T<:Tuple} = convert(T, Tuple(B))
  Base.convert(::Type{CartesianIndices}, B::BoundingBox{<:Integer}) =
      CartesianIndices(B)
  Base.convert(::Type{CartesianIndices{2}}, B::BoundingBox{<:Integer}) =
      CartesianIndices(B)
  ```

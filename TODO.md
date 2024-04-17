- Remove/deprecate `TwoDimensional.Suffixed`.

- Use `IntervalSets`.

- `BoundingBox(pnt::Point) = BoundingBox(pnt, pnt)` and
  `BoundingBox{T}(pnt::Point) where {T} = BoundingBox{T}(pnt, pnt)`

- Add geometrical shape objects: `Rectangle`, `Polygon`, `Circle`, and `Ellipse`.

- Add masks.

- Iterating and indexing bounding-boxes yield `Point(xmin,ymin)` and
  `Point(xmax,ymax)`.

- Deprecate/discard weighted points or weight may have different units than
  coordinates.

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

- `interior` and `exterior` should follow the same rules as `nearest`.
- Make `inv(A)` faster (and a bit less accurate) by avoiding division.
- Doc. for points and bounding boxes.
- Improve bounding box algorithm with a predicate.

- Extending `convert` for other prurpose of type conversion is potentially
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

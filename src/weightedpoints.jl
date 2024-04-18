# Other constructors of weighted points.
WeightedPoint{T}(w, x, y) where {T<:AbstractFloat} = WeightedPoint{T}((w, x, y))
WeightedPoint(w::Real, x::Real, y::Real) = WeightedPoint(map(float, promote(w, x, y)))
WeightedPoint(w::T, x::T, y::T) where {T<:AbstractFloat} = WeightedPoint{T}(w, x, y)

WeightedPoint(; w::Real, x::Real, y::Real) = WeightedPoint(w, x, y)
WeightedPoint{T}(; w::Real, x::Real, y::Real) where {T<:AbstractFloat} = WeightedPoint{T}(w, x, y)

WeightedPoint(pnt::WeightedPoint) = pnt
WeightedPoint{T}(pnt::WeightedPoint{T}) where {T<:AbstractFloat} = pnt
WeightedPoint{T}(pnt::WeightedPoint) where {T<:AbstractFloat} = WeightedPoint{T}(Tuple(pnt))

WeightedPoint(pnt::Point) = WeightedPoint(oneunit(eltype(pnt)), pnt...)
WeightedPoint{T}(pnt::Point) where {T<:AbstractFloat} = WeightedPoint(oneunit(T), pnt...)

# Properties.
Base.propertynames(::WeightedPoint) = (:w, :x, :y)
Base.getproperty(pnt::WeightedPoint, key::Symbol) =
    key === :w ? pnt[1] :
    key === :x ? pnt[2] :
    key === :y ? pnt[3] : throw(KeyError(key))

function Base.show(io::IO, pnt::WeightedPoint{T}) where {T}
    print(io, "WeightedPoint{")
    show(io, T)
    print(io, "}(w = "); show(io, pnt.w)
    print(io, ", x = "); show(io, pnt.x)
    print(io, ", y = "); show(io, pnt.y)
    print(io, ")")
end

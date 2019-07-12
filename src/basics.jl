#
# basics.jl --
#
# Utilities for Julia interface to TAO.
#
#-------------------------------------------------------------------------------
#
# This file if part of the TAO software (https://github.com/emmt/TAO) licensed
# under the MIT license.
#
# Copyright (C) 2019, Éric Thiébaut.
#

"""

`dimensions(A)` yields the list of dimensions of array `A`.  The result is the
same as `size(A)` but it is also checked that `A` has standard indices
(starting at 1).

"""
function dimensions(A::AbstractArray{<:Any,N}) where {N}
    # The fastest code is LazyAlgebra.has_standard_indexing(A)
    LazyAlgebra.has_standard_indexing(A) ||
        error("array has non-standard indexing")
    return size(A)
end

#------------------------------------------------------------------------------
# WEIGHTED ARRAYS

"""
```julia
weights(A)
```

yields the array of weights of the weighted array `A`.

Also see [`WeightedArray`](@ref), [`values`](@ref).

"""
weights(A::WeightedArray) = A.wgt

"""
```julia
values(A)
```

yields the array of values of the weighted array `A`.

Also see [`WeightedArray`](@ref), [`weights`](@ref).

"""
values(A::WeightedArray) = A.dat

#------------------------------------------------------------------------------
# POINTS AND BOUNDING BOXES

# Constructors of points and conversion to/from a Cartesian index.
Point(P::Point) = P
Point{T}(P::Point{T}) where {T} = P
Point{T}(P::Point) where {T} = Point{T}(P.x, P.y)
Point(x::Tx, y::Ty) where {Tx<:Real,Ty<:Real} =
    Point{promote_type{Tx,Ty}}(x, y)
Point(I::CartesianIndex{2}) = Point(I[1], I[2])
Base.CartesianIndex(P::Point{<:Integer}) = CartesianIndex(P.x, P.y)
Base.convert(::Type{CartesianIndex}, P::Point{<:Integer}) = CartesianIndex(P)
Base.convert(::Type{CartesianIndex{2}}, P::Point{<:Integer}) = CartesianIndex(P)

# Conversion of points and rounding to nearest integer coordinates.
Base.convert(::Type{T}, P::Point) where {T<:Point} = T(P)
Base.round(::Type{Point{T}}, P::Point) where {T<:Real} =
    Point(round(T, P.x), round(T, P.y))

# Constructors of bounding boxes and conversion.
function BoundingBox(xmin::Txmin, xmax::Txmax,
                     ymin::Tymin, ymax::Tymax) where {Txmin,Txmax,Tymin,Tymax}
    T = promote_type(Txmin, Txmax, Tymin, Tymax)
    return BoundingBox{T}(xmin, xmax, ymin, ymax)
end
BoundingBox(B::BoundingBox) = B
BoundingBox{T}(B::BoundingBox{T}) where {T} = B
BoundingBox{T}(B::BoundingBox) where {T} =
    BoundingBox{T}(B.xmin, B.xmax, B.ymin, B.ymax)
BoundingBox(I0::CartesianIndex{2}, I1::CartesianIndex{2}) =
    BoundingBox(I0[1], I1[1], I0[2], I1[2])
BoundingBox(P0::Point, P1::Point) =
    BoundingBox(P0.x, P1.x, P0.y, P1.y)

# Empty bounding and unlimited boxes.
BoundingBox{T}(::Nothing) where {T<:Real} = typemin(BoundingBox{T})
Base.typemin(::Type{BoundingBox{T}}) where {T<:Real} =
    BoundingBox(typemax(T), typemin(T), typemax(T), typemin(T))
Base.typemax(::Type{BoundingBox{T}}) where {T<:Real} =
    BoundingBox(typemin(T), typemax(T), typemin(T), typemax(T))

# Conversion of bounding boxes to/from Cartesian indices.
Base.CartesianIndices(B::BoundingBox{<:Integer}) =
    CartesianIndices((Int(B.xmin):Int(B.xmax), Int(B.ymin):Int(B.ymax)))
Base.convert(::Type{CartesianIndices}, B::BoundingBox{<:Integer}) =
    CartesianIndices(B)
Base.convert(::Type{CartesianIndices{2}}, B::BoundingBox{<:Integer}) =
    CartesianIndices(B)

# Conversion of bounding boxes and rounding to nearest integer coordinates.
Base.convert(::Type{T}, B::BoundingBox) where {T<:BoundingBox} = T(B)
Base.round(::Type{BoundingBox{T}}, B::BoundingBox) where {T<:Real} =
    BoundingBox(round(T, B.xmin), round(T, B.xmax),
                round(T, B.ymin), round(T, B.ymax))

# Lower left and upper right corners of a bounding box.
Base.first(B::BoundingBox) = Point(B.xmin, B.ymin)
Base.last(B::BoundingBox) = Point(B.xmax, B.ymax)

Base.size(B::BoundingBox{Int}) = (max(B.xmax - B.xmin + 1, 0),
                                  max(B.ymax - B.ymin + 1, 0))
Base.size(B::BoundingBox{Int}, k::Integer) =
    (k == 1 ? max(B.xmax - B.xmin + 1, 0) :
     k == 2 ? max(B.ymax - B.ymin + 1, 0) :
     k > 2 ? 1 : throw_bad_dimension_index())

Base.size(B::BoundingBox{<:Integer}) = (max(Int(B.xmax) - Int(B.xmin) + 1, 0),
                                        max(Int(B.ymax) - Int(B.ymin) + 1, 0))
Base.size(B::BoundingBox{<:Integer}, k::Integer) =
    (k == 1 ? max(Int(B.xmax) - Int(B.xmin) + 1, 0) :
     k == 2 ? max(Int(B.ymax) - Int(B.ymin) + 1, 0) :
     k > 2 ? 1 : throw_bad_dimension_index())

Base.axes(B::BoundingBox{Int}) = (B.xmin:B.xmax, B.ymin:B.ymax)
Base.axes(B::BoundingBox{Int}, k::Integer) =
    (k == 1 ? (B.xmin:B.xmax) :
     k == 2 ? (B.ymin:B.ymax) :
     k > 2 ? Base.OneTo(1) : throw_bad_dimension_index())

Base.axes(B::BoundingBox{<:Integer}) = (Int(B.xmin):Int(B.xmax),
                                        Int(B.ymin):Int(B.ymax))
Base.axes(B::BoundingBox{<:Integer}, k::Integer) =
    (k == 1 ? (Int(B.xmin):Int(B.xmax)) :
     k == 2 ? (Int(B.ymin):Int(B.ymax)) :
     k > 2 ? Base.OneTo(1) : throw_bad_dimension_index())

@noinline throw_bad_dimension_index() =
    error("invalid dimension index")

# Union of bounding boxes:
Base.:(∪)(A::BoundingBox, B::BoundingBox) =
    BoundingBox(min(A.xmin, B.xmin), max(A.xmax, B.xmax),
                min(A.ymin, B.ymin), max(A.ymax, B.ymax))

# Intersection of bounding boxes:
Base.:(∩)(A::BoundingBox, B::BoundingBox) =
    BoundingBox(max(A.xmin, B.xmin), min(A.xmax, B.xmax),
                max(A.ymin, B.ymin), min(A.ymax, B.ymax))

# Use bounding boxes to extract a sub-array or a view.
@propagate_inbounds Base.getindex(A::AbstractMatrix, B::BoundingBox{<:Integer}) =
    A[B.xmin:B.xmax, B.ymin:B.ymax]
@propagate_inbounds Base.getindex(A::WeightedMatrix, B::BoundingBox{<:Integer}) =
    WeightedMatrix(weights(A)[B], values(A)[B])

Base.view(A::AbstractMatrix, B::BoundingBox{<:Integer}) =
    view(A, B.xmin:B.xmax, B.ymin:B.ymax)
Base.view(A::WeightedMatrix, B::BoundingBox{<:Integer}) =
    WeightedMatrix(view(weights(A), B), view(values(A), B))

"""

`interior([BoundingBox{T},] B::BoundingBox)` yields the largest bounding box
with integer valued bounds which is contained by box `B`.  Optional first
argument is to specify the type of the result which is that of `B` by default.

"""
interior(B::BoundingBox{<:Integer}) = B
interior(B::BoundingBox{<:AbstractFloat}) =
    BoundingBox(ceil(B.xmin), floor(B.xmax),
                ceil(B.ymin), floor(B.ymax))
interior(::Type{BoundingBox{T}}, B::BoundingBox{<:AbstractFloat}) where {T<:Integer} =
    BoundingBox(ceil(T, B.xmin), floor(T, B.xmax),
                ceil(T, B.ymin), floor(T, B.ymax))
interior(::Type{BoundingBox{T}}, B::BoundingBox) where {T} =
    BoundingBox{T}(interior(B))

"""

`exterior([BoundingBox{T},] B::BoundingBox)` yields the smallest boundingbox
box with integer valued bounds which contains box `B`.  Optional first argument
is to specify the type of the result which is that of `B` by default.

"""
exterior(B::BoundingBox{<:Integer}) = B
exterior(B::BoundingBox{<:AbstractFloat}) =
    BoundingBox(floor(B.xmin), ceil(B.xmax),
                floor(B.ymin), ceil(B.ymax))
exterior(::Type{BoundingBox{T}}, B::BoundingBox{<:AbstractFloat}) where {T<:Integer} =
    BoundingBox(floor(T, B.xmin), ceil(T, B.xmax),
                floor(T, B.ymin), ceil(T, B.ymax))
exterior(::Type{BoundingBox{T}}, B::BoundingBox) where {T} =
    BoundingBox{T}(exterior(B))

"""

```julia
center(B::BoundingBox) -> c::Point
```

yields the central point of the bounding box `B`.

"""
center(B::BoundingBox{<:Integer}) =
    Point(0.5*(B.xmin + B.xmax), 0.5*(B.ymin + B.ymax))
center(B::BoundingBox{T}) where {T<:AbstractFloat} =
    Point(half(T)*(B.xmin + B.xmax), half(T)*(B.ymin + B.ymax))

half(::Type{T}) where {T<:AbstractFloat} = one(T)/convert(T, 2)

"""

`area(B)` yields the area of the bounding box `B`.

"""
area(B::BoundingBox{T}) where {T<:Real} =
    max(B.xmax - B.xmin, zero(T))*max(B.ymax - B.ymin, zero(T))

# Scaling of a point coordinates.
Base.:(*)(P::Point, α::Real) = α*P
Base.:(*)(α::Real, P::Point) = Point(α*P.x, α*P.y)
Base.:(/)(P::Point, α::Real) = Point(P.x/α, P.y/α)
Base.:(\)(α::Real, P::Point) = P/α

# Unary minus applied to point negate coordinates.
Base.:(-)(P::Point{<:Union{AbstractFloat,Signed,Irrational}}) =
    Point(-P.x, -P.y)

# Addition and subtraction of point coordinates.
Base.:(+)(A::Point, B::Point) = Point(A.x + B.x, A.y + B.y)
Base.:(-)(A::Point, B::Point) = Point(A.x - B.x, A.y - B.y)

# Scaling of bounding box bounds (e.g. to change units).
Base.:(*)(B::BoundingBox, α::Real) = α*B
Base.:(*)(α::Real, B::BoundingBox) = BoundingBox(α*B.xmin, α*B.xmax,
                                                 α*B.ymin, α*B.ymax)
Base.:(/)(B::BoundingBox{T}, α::Real) where {T} = (one(T)/α)*B
Base.:(\)(α::Real, B::BoundingBox) = B/α

# Add or remove a margin δ to a bounding box B.
Base.:(+)(B::BoundingBox, δ::Real) =
    BoundingBox(B.xmin - δ, B.xmax + δ,
                B.ymin - δ, B.ymax + δ)
Base.:(-)(B::BoundingBox, δ::Real) =
    BoundingBox(B.xmin + δ, B.xmax - δ,
                B.ymin + δ, B.ymax - δ)

# Translate a bounding box.
Base.:(+)(B::BoundingBox, P::Point) =
    BoundingBox(B.xmin + P.x, B.xmax + P.x,
                B.ymin + P.y, B.ymax + P.y)
Base.:(-)(B::BoundingBox, P::Point) =
    BoundingBox(B.xmin - P.x, B.xmax - P.x,
                B.ymin - P.y, B.ymax - P.y)

#------------------------------------------------------------------------------

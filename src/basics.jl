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
# POINTS AND BOUNDING BOXES

# Constructors of points and conversion to/from a Cartesian index.
Point(I::CartesianIndex{2}) = Point(I[1], I[2])
Base.CartesianIndex(P::Point{<:Integer}) = CartesianIndex(P.x, P.y)
Base.convert(::Type{CartesianIndex}, P::Point{<:Integer}) = CartesianIndex(P)
Base.convert(::Type{CartesianIndex{2}}, P::Point{<:Integer}) = CartesianIndex(P)

# Conversion of points and rounding to nearest integer coordinates.
Base.convert(::Type{Point{T}}, p::Point{T}) where {T<:Real} = p
Base.convert(::Type{Point{T}}, p::Point) where {T<:Real} =
    Point(T(p.x), T(p.y))
Base.round(::Type{Point{T}}, p::Point) where {T<:Real} =
    Point(round(T, p.x), round(T, p.y))

# Constructors of bounding boxes and conversion to/from Cartesian indices.
BoundingBox(I0::CartesianIndex{2}, I1::CartesianIndex{2}) =
    BoundingBox(I0[1], I1[1], I0[2], I1[2])
BoundingBox(P0::Point, P1::Point) =
    BoundingBox(P0.x, P1.x, P0.y, P1.y)
Base.CartesianIndices(B::BoundingBox{<:Integer}) =
    CartesianIndices((Int(B.xmin):Int(B.xmax), Int(B.ymin):Int(B.ymax)))
Base.convert(::Type{CartesianIndices}, B::BoundingBox{<:Integer}) =
    CartesianIndices(B)
Base.convert(::Type{CartesianIndices{2}}, B::BoundingBox{<:Integer}) =
    CartesianIndices(B)

# Conversion of bounding boxes and rounding to nearest integer coordinates.
Base.convert(::Type{BoundingBox{T}}, b::BoundingBox{T}) where {T<:Real} = b
Base.convert(::Type{BoundingBox{T}}, b::BoundingBox) where {T<:Real} =
    BoundingBox(T(b.xmin), T(b.ymin), T(b.xmax), T(b.ymax))
Base.round(::Type{BoundingBox{T}}, b::BoundingBox) where {T<:Real} =
    BoundingBox(round(T, b.xmin), round(T, b.xmax),
                round(T, b.ymin), round(T, b.ymax))

# Lower left and upper right corners of a bounding box.
Base.first(b::BoundingBox) = Point(b.xmin, b.ymin)
Base.last(b::BoundingBox) = Point(b.xmax, b.ymax)

Base.size(b::BoundingBox{Int}) = (max(b.xmax - b.xmin + 1, 0),
                                  max(b.ymax - b.ymin + 1, 0))
Base.size(b::BoundingBox{Int}, k::Integer) =
    (k == 1 ? max(b.xmax - b.xmin + 1, 0) :
     k == 2 ? max(b.ymax - b.ymin + 1, 0) :
     k > 2 ? 1 : error("invalid dimension index"))

Base.size(b::BoundingBox{<:Integer}) = (max(Int(b.xmax) - Int(b.xmin) + 1, 0),
                                        max(Int(b.ymax) - Int(b.ymin) + 1, 0))
Base.size(b::BoundingBox{<:Integer}, k::Integer) =
    (k == 1 ? max(Int(b.xmax) - Int(b.xmin) + 1, 0) :
     k == 2 ? max(Int(b.ymax) - Int(b.ymin) + 1, 0) :
     k > 2 ? 1 : error("invalid dimension index"))

# Union of bounding boxes:
Base.:(∪)(a::BoundingBox, b::BoundingBox) =
    BoundingBox(min(a.xmin, b.xmin), max(a.xmax, b.xmax),
                min(a.ymin, b.ymin), max(a.ymax, b.ymax))

# Intersection of bounding boxes:
Base.:(∩)(a::BoundingBox, b::BoundingBox) =
    BoundingBox(max(a.xmin, b.xmin), min(a.xmax, b.xmax),
                max(a.ymin, b.ymin), min(a.ymax, b.ymax))

"""

`interior([BoundingBox{T},] b::BoundingBox)` yields the largest bounding box
with integer valued bounds which is contained by box `b`.  Optional first
argument is to specify the type of the result which is that of `b` by default.

"""
interior(b::BoundingBox{<:Integer}) = b
interior(b::BoundingBox{T}) where {T<:Real} =
    BoundingBox(ceil(T, b.xmin), floor(T, b.xmax),
                ceil(T, b.ymin), floor(T, b.ymax))
interior(::Type{BoundingBox{T}}, b::BoundingBox) where {T<:Real} =
    BoundingBox(ceil(T, b.xmin), floor(T, b.xmax),
                ceil(T, b.ymin), floor(T, b.ymax))
interior(::Type{BoundingBox{T}}, b::BoundingBox{<:Integer}) where {T<:Real} =
    convert(BoundingBox{T}, b)

"""

`exterior([BoundingBox{T},] b::BoundingBox)` yields the smallest boundingbox
box with integer valued bounds which contains box `b`.  Optional first argument
is to specify the type of the result which is that of `b` by default.

"""
exterior(b::BoundingBox{<:Integer}) = b
exterior(b::BoundingBox{T}) where {T<:Real} =
    BoundingBox(floor(T, b.xmin), ceil(T, b.xmax),
                floor(T, b.ymin), ceil(T, b.ymax))
exterior(::Type{BoundingBox{T}}, b::BoundingBox) where {T<:Real} =
    BoundingBox(floor(T, b.xmin), ceil(T, b.xmax),
                floor(T, b.ymin), ceil(T, b.ymax))
exterior(::Type{BoundingBox{T}}, b::BoundingBox{<:Integer}) where {T<:Real} =
    convert(BoundingBox{T}, b)

"""

```julia
center(b::BoundingBox) -> c::Point
```

yields the central point of the bounding box `b`.

"""
center(b::BoundingBox{<:Integer}) =
    Point(0.5*(b.xmin + b.xmax), 0.5*(b.ymin + b.ymax))
center(b::BoundingBox{T}) where {T<:AbstractFloat} =
    Point(half(T)*(b.xmin + b.xmax), half(T)*(b.ymin + b.ymax))

half(::Type{T}) where {T<:AbstractFloat} = one(T)/convert(T, 2)

"""

`area(b)` yields the area of the bounding box `b`.

"""
area(b::BoundingBox{T}) where {T<:Real} =
    max(b.xmax - b.xmin, zero(T))*max(b.ymax - b.ymin, zero(T))

# Scaling of a point coordinates.
Base.:(*)(p::Point, α::Real) = α*p
Base.:(*)(α::Real, p::Point) = Point(α*p.x, α*p.y)
Base.:(/)(p::Point, α::Real) = Point(p.x/α, p.y/α)
Base.:(\)(α::Real, p::Point) = p/α

# Unary minus applied to point negate coordinates.
Base.:(-)(P::Point{<:Union{AbstractFloat,Signed,Irrational}}) =
    Point(-P.x, -P.y)

# Addition and subtraction of point coordinates.
Base.:(+)(a::Point, b::Point) = Point(a.x + b.x, a.y + b.y)
Base.:(-)(a::Point, b::Point) = Point(a.x - b.x, a.y - b.y)

# Add or remove a margin δ to a bounding box b.
Base.:(+)(b::BoundingBox, δ::Real) =
    BoundingBox(b.xmin - δ, b.xmax + δ,
                b.ymin - δ, b.ymax + δ)
Base.:(-)(b::BoundingBox, δ::Real) =
    BoundingBox(b.xmin + δ, b.xmax - δ,
                b.ymin + δ, b.ymax - δ)

# Translate a bounding box.
Base.:(+)(b::BoundingBox, p::Point) =
    BoundingBox(b.xmin + p.x, b.xmax + p.x,
                b.ymin + p.y, b.ymax + p.y)
Base.:(-)(b::BoundingBox, p::Point) =
    BoundingBox(b.xmin - p.x, b.xmax - p.x,
                b.ymin - p.y, b.ymax - p.y)

#------------------------------------------------------------------------------

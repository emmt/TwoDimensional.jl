#
# basics.jl --
#
# Basic types.
#
#-------------------------------------------------------------------------------
#
# This file is part of the TwoDimensional Julia package licensed under the MIT
# license (https://github.com/emmt/TwoDimensional.jl).
#
# Copyright (C) 2019, Éric Thiébaut.
#

using Base: @propagate_inbounds

"""

Any object whose type is derived from `AbstractPoint{T}` has at least 2 fields:
`x` its abscissa and `y` its ordinate, both of type `T`.

See also: [`Point`](@ref), [`WeightedPoint`](@ref).

"""
abstract type AbstractPoint{T<:Real} end

"""

```julia
Point(x,y)
```

yields an instance of a 2D point of coordinates `(x,y)`.

A point may be multiplied or divided by a scalar to scale its coordinates.  The
addition (resp. subtraction) of two points adds (resp. subtracts) their
coordinates.

Coordinates can be specified by keywords:

```julia
Point(x=xval, y=yval)
```

There are no default values for keywords `x` and `y` so both must be specified.

See also: [`WeightedPoint`](@ref), [`AbstractPoint`](@ref).

""" Point

struct Point{T} <: AbstractPoint{T}
    x::T
    y::T
end

"""

A `WeightedPoint{T}` has just 3 fields: `w` its weight, `x` its abscissa and
`y` its ordinate, all of type `T`.  By convention `w ≥ 0` but this is not
checked for efficiency reasons.

See also: [`Point`](@ref), [`AbstractPoint`](@ref).

""" WeightedPoint

struct WeightedPoint{T<:AbstractFloat}  <: AbstractPoint{T}
    w::T # weight
    x::T # abscissa
    y::T # ordinate
end

"""

`BoundingBox(xmin,xmax,ymin,ymax)` yields an instance of a 2D rectangular
bounding-box whose sides are aligned with the coordinate axes and containing
points of coordinates `(x,y)` such that `xmin ≤ x ≤ xmax` and `ymin ≤ y ≤
ymax`.  The box is *empty* if `xmin > xmax` or `ymin > ymax`.

A bounding-box can be constructed from the first and last points (i.e. at the
lower-left and upper right opposite corners) of the box:

```julia
BoundingBox(P0::Point, P1::Point)
BoundingBox(I0::CartesianIndex{2}, I1::CartesianIndex{2})
```

Coordinates can be specified by keywords:

```julia
BoundingBox(xmin=x0, ymin=y0, xmax=x1, ymax=y1)
```

There are no default values for keywords `xmin`, `xmax`, `ymin` and `ymax` so
all must be specified.

See also [`Point`](@ref), [`interior`](@ref), [`exterior`](@ref).

""" BoundingBox

struct BoundingBox{T<:Real}
    xmin::T
    xmax::T
    ymin::T
    ymax::T
end

# Constructors of points and conversion to/from a Cartesian index.
Point(P::Point) = P
Point{T}(P::Point{T}) where {T} = P
Point{T}(P::Point) where {T} = Point{T}(P.x, P.y)
Point(x::Tx, y::Ty) where {Tx<:Real,Ty<:Real} = Point{promote_type(Tx,Ty)}(x, y)
Point(; x::Real, y::Real) = Point(x, y)
Point(I::CartesianIndex{2}) = Point(I[1], I[2])
Point(P::NTuple{2,Real}) = Point(P[1], P[2])
Base.CartesianIndex(P::Point{<:Integer}) = CartesianIndex(P.x, P.y)
Base.Tuple(P::Point) = (P.x, P.y)
Base.eltype(::AbstractPoint{T}) where {T} = T
Broadcast.broadcasted(::Type{T}, obj::Point) where {T<:Real} =
    Point{T}(obj)
Base.promote_type(::Type{Point{T}}, ::Type{Point{U}}) where {T,U} =
    Point{promote_type(T,U)}
Base.promote_type(::Type{Point{T}}, ::Type{Point{T}}) where {T} =
    Point{T}

# Constructors of weighted points.
WeightedPoint(P::WeightedPoint) = P
WeightedPoint{T}(P::WeightedPoint{T}) where {T<:AbstractFloat} = P
WeightedPoint{T}(P::WeightedPoint) where {T} = WeightedPoint{T}(P.w, P.x, P.y)
WeightedPoint(w::Tw, x::Tx, y::Ty) where {Tw<:Real,Tx<:Real,Ty<:Real} =
    WeightedPoint{float(promote_type(Tw,Tx,Ty))}(w, x, y)
WeightedPoint(; w::Real, x::Real, y::Real) = WeightedPoint(w, x, y)
WeightedPoint(P::Point{T}) where {T} = WeightedPoint(one(T), P.x, P.y)
WeightedPoint(P::NTuple{3,Real}) = WeightedPoint(P[1], P[2], P[3])
Base.Tuple(P::WeightedPoint) = (P.w, P.x, P.y)
Broadcast.broadcasted(::Type{T}, obj::WeightedPoint) where {T<:AbstractFloat} =
    WeightedPoint{T}(obj)
Base.promote_type(::Type{WeightedPoint{T}}, ::Type{WeightedPoint{U}}) where {T,U} =
    WeightedPoint{promote_type(T,U)}
Base.promote_type(::Type{WeightedPoint{T}}, ::Type{WeightedPoint{T}}) where {T} =
    WeightedPoint{T}

# Conversion of points and rounding to nearest integer coordinates.
Base.convert(::Type{T}, obj::AbstractPoint) where {T<:AbstractPoint} = T(obj)

"""

```julia
round([T,] obj)
```

yields the object that is the nearest to `obj` by rounding its coordinates to
the nearest integer.  Argument `T` can be the type of the result (a point or a
bounding-box) or the type of the coordinates of the result.

For points, see also: [`floor`](@ref), [`ceil`](@ref).

For bounding-boxes, see also: [`interior`](@ref), [`exterior`](@ref).

"""
Base.round(obj::Point{T}) where {T} = round(T, obj)
Base.round(::Type{Point{T}}, obj::Point) where {T} = round(T, obj)
Base.round(::Type{T}, obj::Point{T}) where {T<:Integer} = obj
Base.round(::Type{T}, obj::Point{T}) where {T<:Real} =
    Point(round(obj.x, RoundNearest),
          round(obj.y, RoundNearest))
Base.round(::Type{T}, obj::Point{<:Integer}) where {T<:Integer} =
    Point{T}(obj)
Base.round(::Type{T}, obj::Point{<:Real}) where {T<:Integer} =
    Point(round(T, obj.x, RoundNearest),
          round(T, obj.y, RoundNearest))
Base.round(::Type{T}, obj::Point{<:Integer}) where {T<:Real} =
    Point{T}(obj)
Base.round(::Type{T}, obj::Point{U}) where {T<:Real,U<:Real} =
    Point{T}(round(U, obj))

# Extend round for bounding-boxes.
Base.round(obj::BoundingBox{T}) where {T} = round(T, obj)
Base.round(::Type{BoundingBox{T}}, obj::BoundingBox) where {T} = round(T, obj)
Base.round(::Type{T}, obj::BoundingBox{T}) where {T<:Integer} = obj
Base.round(::Type{T}, obj::BoundingBox{T}) where {T<:Real} =
    BoundingBox(round(obj.xmin, RoundNearest),
                round(obj.xmax, RoundNearest),
                round(obj.ymin, RoundNearest),
                round(obj.ymax, RoundNearest))
Base.round(::Type{T}, obj::BoundingBox{<:Integer}) where {T<:Integer} =
    BoundingBox{T}(obj)
Base.round(::Type{T}, obj::BoundingBox{<:Real}) where {T<:Integer} =
    BoundingBox(round(T, obj.xmin, RoundNearest),
                round(T, obj.xmax, RoundNearest),
                round(T, obj.ymin, RoundNearest),
                round(T, obj.ymax, RoundNearest))
Base.round(::Type{T}, obj::BoundingBox{<:Integer}) where {T<:Real} =
    BoundingBox{T}(obj)
Base.round(::Type{T}, obj::BoundingBox{U}) where {T<:Real,U<:Real} =
    BoundingBox{T}(round(U, obj))

"""

```julia
floor([T,] P)
```

yields the point with the largest integer coordinates smaller or equal those of
the point `P`.  Argument `T` can be the type of the result or the type of the
coordinates of the result.

See also: [`round`](@ref), [`ceil`](@ref).

"""
Base.floor(obj::Point{T}) where {T} = floor(T, obj)
Base.floor(::Type{Point{T}}, obj::Point) where {T} = floor(T, obj)
Base.floor(::Type{T}, obj::Point{T}) where {T<:Integer} = obj
Base.floor(::Type{T}, obj::Point{T}) where {T<:Real} =
    Point(floor(obj.x), floor(obj.y))
Base.floor(::Type{T}, obj::Point{<:Integer}) where {T<:Integer} =
    Point{T}(obj)
Base.floor(::Type{T}, obj::Point{<:Real}) where {T<:Integer} =
    Point(floor(T, obj.x), floor(T, obj.y))
Base.floor(::Type{T}, obj::Point{<:Integer}) where {T<:Real} =
    Point{T}(obj)
Base.floor(::Type{T}, obj::Point{U}) where {T<:Real,U<:Real} =
    Point{T}(floor(U, obj))

"""

```julia
ceil([T,] P)
```

yields the point with the smallest integer coordinates larger or equal those of
the point `P`.  Argument `T` can be the type of the result or the type of the
coordinates of the result.

See also: [`round`](@ref), [`floor`](@ref).

"""
Base.ceil(obj::Point{T}) where {T} = ceil(T, obj)
Base.ceil(::Type{Point{T}}, obj::Point) where {T} = ceil(T, obj)
Base.ceil(::Type{T}, obj::Point{T}) where {T<:Integer} = obj
Base.ceil(::Type{T}, obj::Point{T}) where {T<:Real} =
    Point(ceil(obj.x), ceil(obj.y))
Base.ceil(::Type{T}, obj::Point{<:Integer}) where {T<:Integer} =
    Point{T}(obj)
Base.ceil(::Type{T}, obj::Point{<:Real}) where {T<:Integer} =
    Point(ceil(T, obj.x), ceil(T, obj.y))
Base.ceil(::Type{T}, obj::Point{<:Integer}) where {T<:Real} =
    Point{T}(obj)
Base.ceil(::Type{T}, obj::Point{U}) where {T<:Real,U<:Real} =
    Point{T}(ceil(U, obj))

# Methods hypot() and atan() yield the polar coordinates of a point.
Base.hypot(P::Point) = hypot(P.x, P.y)
Base.atan(P::Point) = atan(P.y, P.x)

"""

```julia
distance(A, B)
```

yields the Euclidean distance between the 2 points `A` and `B`.

"""
distance(A::Point, B::Point) =
    distance(promote(A, B)...)
distance(A::Point{T}, B::Point{T}) where {T<:Unsigned} =
    hypot(ifelse(A.x > B.x, A.x - B.x, B.x - A.x),
          ifelse(A.y > B.y, A.y - B.y, B.y - A.y))
distance(A::Point{T}, B::Point{T}) where {T<:Real} =
    hypot(B.x - A.x, B.y - A.y)

# Constructors of bounding-boxes and conversion.
function BoundingBox(xmin::Txmin, xmax::Txmax,
                     ymin::Tymin, ymax::Tymax) where {Txmin,Txmax,Tymin,Tymax}
    T = promote_type(Txmin, Txmax, Tymin, Tymax)
    return BoundingBox{T}(xmin, xmax, ymin, ymax)
end
BoundingBox(; xmin::Real, xmax::Real, ymin::Real, ymax::Real) =
    BoundingBox(xmin, xmax, ymin, ymax)
BoundingBox(B::BoundingBox) = B
BoundingBox{T}(B::BoundingBox{T}) where {T<:Real} = B
BoundingBox{T}(B::BoundingBox) where {T<:Real} =
    BoundingBox{T}(B.xmin, B.xmax, B.ymin, B.ymax)
BoundingBox(I0::CartesianIndex{2}, I1::CartesianIndex{2}) =
    BoundingBox(I0[1], I1[1], I0[2], I1[2])
BoundingBox(P0::AbstractPoint, P1::AbstractPoint) =
    BoundingBox(P0.x, P1.x, P0.y, P1.y)
BoundingBox(P::NTuple{4,Real}) = BoundingBox(P[1], P[2], P[3], P[4])
Base.eltype(::BoundingBox{T}) where {T} = T
Base.Tuple(B::BoundingBox) = (B.xmin, B.xmax, B.ymin, B.ymax)
Broadcast.broadcasted(::Type{T}, obj::BoundingBox) where {T<:Real} =
    BoundingBox{T}(obj)
Base.promote_type(::Type{BoundingBox{T}}, ::Type{BoundingBox{U}}) where {T,U} =
    BoundingBox{promote_type(T,U)}
Base.promote_type(::Type{BoundingBox{T}}, ::Type{BoundingBox{T}}) where {T} =
    BoundingBox{T}

BoundingBox(A::AbstractMatrix) = BoundingBox(axes(A))
BoundingBox(inds::NTuple{2,AbstractUnitRange{<:Integer}}) =
    BoundingBox(inds[1], inds[2])
BoundingBox(X::AbstractUnitRange{<:Integer}, Y::AbstractUnitRange{<:Integer}) =
    BoundingBox(Int(first(X)), Int(last(X)), Int(first(Y)), Int(last(Y)))

# FIXME: speed-up stupid algorithm!
function BoundingBox(f::Function, A::AbstractMatrix)
    inds = axes(A)
    imin = typemax(Int)
    imax = typemin(Int)
    jmin = typemax(Int)
    jmax = typemin(Int)
    for j in inds[2], i in inds[1]
        if f(A[i,j])
            imin = min(imin, i)
            imax = max(imax, i)
            jmin = min(jmin, j)
            jmax = max(jmax, j)
        end
    end
    return BoundingBox(imin, imax, jmin, jmax)
end

# Empty bounding and unlimited boxes.
BoundingBox{T}(::Nothing) where {T<:Real} = typemin(BoundingBox{T})
Base.typemin(::Type{BoundingBox{T}}) where {T<:Real} =
    BoundingBox(typemax(T), typemin(T), typemax(T), typemin(T))
Base.typemax(::Type{BoundingBox{T}}) where {T<:Real} =
    BoundingBox(typemin(T), typemax(T), typemin(T), typemax(T))

# Conversion of bounding-boxes to/from Cartesian indices.
Base.CartesianIndices(B::BoundingBox{<:Integer}) =
    CartesianIndices((Int(B.xmin):Int(B.xmax), Int(B.ymin):Int(B.ymax)))

# Conversion of bounding-boxes and rounding to nearest integer coordinates.
Base.convert(::Type{T}, obj::T) where {T<:BoundingBox} = obj
Base.convert(::Type{T}, obj::BoundingBox) where {T<:BoundingBox} = T(obj)

# Lower left and upper right corners of a bounding-box.
Base.first(B::BoundingBox) = Point(B.xmin, B.ymin)
Base.last(B::BoundingBox) = Point(B.xmax, B.ymax)

Base.isempty(B::BoundingBox) = ((B.xmin > B.xmax)|(B.ymin > B.ymax))

Base.size(B::BoundingBox{Int}) =
    (isempty(B) ? (0, 0) : (B.xmax - B.xmin + 1,
                            B.ymax - B.ymin + 1))
Base.size(B::BoundingBox{Int}, k::Integer) =
    (isempty(B) ? (k ≥ 1 ? 0 : throw_bad_dimension_index()) :
     k == 1 ? B.xmax - B.xmin + 1 :
     k == 2 ? B.ymax - B.ymin + 1 :
     k > 2 ? 1 : throw_bad_dimension_index())

Base.size(B::BoundingBox{<:Integer}) =
    (isempty(B) ? (0, 0) : (Int(B.xmax) - Int(B.xmin) + 1,
                            Int(B.ymax) - Int(B.ymin) + 1))

Base.size(B::BoundingBox{<:Integer}, k::Integer) =
    (isempty(B) ? (k ≥ 1 ? 0 : throw_bad_dimension_index()) :
     k == 1 ? Int(B.xmax) - Int(B.xmin) + 1 :
     k == 2 ? Int(B.ymax) - Int(B.ymin) + 1 :
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

# Union of bounding-boxes:
Base.:(∪)(A::BoundingBox, B::BoundingBox) =
    BoundingBox(min(A.xmin, B.xmin), max(A.xmax, B.xmax),
                min(A.ymin, B.ymin), max(A.ymax, B.ymax))

# Intersection of bounding-boxes:
Base.:(∩)(A::BoundingBox, B::BoundingBox) =
    BoundingBox(max(A.xmin, B.xmin), min(A.xmax, B.xmax),
                max(A.ymin, B.ymin), min(A.ymax, B.ymax))

# Use bounding-boxes to extract a sub-array or a view.
@propagate_inbounds Base.getindex(A::AbstractMatrix, B::BoundingBox{<:Integer}) =
    A[B.xmin:B.xmax, B.ymin:B.ymax]

Base.view(A::AbstractMatrix, B::BoundingBox{<:Integer}) =
    view(A, B.xmin:B.xmax, B.ymin:B.ymax)

"""

```julia
interior([T,] B)
```

yields the largest bounding -box with integer valued bounds which is contained
by the bounding-box `B`.  Optional argument `T` is to specify the type of the
result or of the coordinates of the result which is the same as `B` by default.

See also: [`exterior`](@ref), [`round`](@ref).

"""
interior(obj::BoundingBox{T}) where {T} = interior(T, obj)
interior(::Type{BoundingBox{T}}, obj::BoundingBox) where {T} = interior(T, obj)
interior(::Type{T}, obj::BoundingBox{T}) where {T<:Integer} = obj
interior(::Type{T}, obj::BoundingBox{T}) where {T<:Real} =
    BoundingBox(ceil(obj.xmin), floor(obj.xmax),
                ceil(obj.ymin), floor(obj.ymax))
interior(::Type{T}, obj::BoundingBox{<:Integer}) where {T<:Integer} =
    BoundingBox{T}(obj)
interior(::Type{T}, obj::BoundingBox{<:Real}) where {T<:Integer} =
    BoundingBox(ceil(T, obj.xmin), floor(T, obj.xmax),
                ceil(T, obj.ymin), floor(T, obj.ymax))
interior(::Type{T}, obj::BoundingBox{<:Integer}) where {T<:Real} =
    BoundingBox{T}(obj)
interior(::Type{T}, obj::BoundingBox{U}) where {T<:Real,U<:Real} =
    BoundingBox{T}(interior(U, obj))

"""

```julia
exterior([T,] B)
```

yields the smallest bounding-box with integer valued bounds which contains the
bounding-box `B`.  Optional argument `T` is to specify the type of the result
or of the coordinates of the result which is the same as `B` by default.

See also: [`interior`](@ref), [`round`](@ref).

"""
exterior(obj::BoundingBox{T}) where {T} = exterior(T, obj)
exterior(::Type{BoundingBox{T}}, obj::BoundingBox) where {T} = exterior(T, obj)
exterior(::Type{T}, obj::BoundingBox{T}) where {T<:Integer} = obj
exterior(::Type{T}, obj::BoundingBox{T}) where {T<:Real} =
    BoundingBox(floor(obj.xmin), ceil(obj.xmax),
                floor(obj.ymin), ceil(obj.ymax))
exterior(::Type{T}, obj::BoundingBox{<:Integer}) where {T<:Integer} =
    BoundingBox{T}(obj)
exterior(::Type{T}, obj::BoundingBox{<:Real}) where {T<:Integer} =
    BoundingBox(floor(T, obj.xmin), ceil(T, obj.xmax),
                floor(T, obj.ymin), ceil(T, obj.ymax))
exterior(::Type{T}, obj::BoundingBox{<:Integer}) where {T<:Real} =
    BoundingBox{T}(obj)
exterior(::Type{T}, obj::BoundingBox{U}) where {T<:Real,U<:Real} =
    BoundingBox{T}(exterior(U, obj))

"""

```julia
center(B::BoundingBox) -> c::Point
```

yields the central point of the bounding-box `B`.

"""
center(B::BoundingBox{<:Integer}) =
    Point(0.5*(B.xmin + B.xmax), 0.5*(B.ymin + B.ymax))
center(B::BoundingBox{T}) where {T<:AbstractFloat} =
    Point(half(T)*(B.xmin + B.xmax), half(T)*(B.ymin + B.ymax))

half(::Type{T}) where {T<:AbstractFloat} = one(T)/convert(T, 2)

"""

`area(B)` yields the area of the bounding-box `B`.

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

# Scaling of bounding-box bounds (e.g. to change units).
Base.:(*)(B::BoundingBox, α::Real) = α*B
Base.:(*)(α::Real, B::BoundingBox) = BoundingBox(α*B.xmin, α*B.xmax,
                                                 α*B.ymin, α*B.ymax)
Base.:(/)(B::BoundingBox{T}, α::Real) where {T} = (one(T)/α)*B
Base.:(\)(α::Real, B::BoundingBox) = B/α

# Add or remove a margin δ to a bounding-box B.
Base.:(+)(B::BoundingBox, δ::Real) =
    BoundingBox(B.xmin - δ, B.xmax + δ,
                B.ymin - δ, B.ymax + δ)
Base.:(-)(B::BoundingBox, δ::Real) =
    BoundingBox(B.xmin + δ, B.xmax - δ,
                B.ymin + δ, B.ymax - δ)

# Translate a bounding-box.
Base.:(+)(B::BoundingBox, P::Point) =
    BoundingBox(B.xmin + P.x, B.xmax + P.x,
                B.ymin + P.y, B.ymax + P.y)
Base.:(-)(B::BoundingBox, P::Point) =
    BoundingBox(B.xmin - P.x, B.xmax - P.x,
                B.ymin - P.y, B.ymax - P.y)

#------------------------------------------------------------------------------

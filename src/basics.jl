#
# basics.jl --
#
# Basic types and method for 2-D points and bounding-boxes.
#
#-------------------------------------------------------------------------------
#
# This file is part of the TwoDimensional Julia package licensed under the MIT
# license (https://github.com/emmt/TwoDimensional.jl).
#
# Copyright (c) 2019-2024, Éric Thiébaut.
#

using Base: @propagate_inbounds

Base.Tuple(pnt::Point) = (pnt.x, pnt.y)
Base.getindex(pnt::Point, i::Integer) =
    i == 1 ? pnt.x :
    i == 2 ? pnt.y : throw(BoundsError(pnt, i))
Base.iterate(pnt::Point, i::Int=1) =
    i == 1 ? (pnt.x, 2) :
    i == 2 ? (pnt.y, 3) : nothing

Base.Tuple(pnt::WeightedPoint) = (pnt.w, pnt.x, pnt.y)
Base.getindex(pnt::WeightedPoint, i::Integer) =
    i == 1 ? pnt.w :
    i == 2 ? pnt.x :
    i == 3 ? pnt.y : throw(BoundsError(pnt, i))
Base.iterate(pnt::WeightedPoint, i::Int=1) =
    i == 1 ? (pnt.w, 2) :
    i == 2 ? (pnt.x, 3) :
    i == 3 ? (pnt.y, 4) : nothing

Base.Tuple(box::BoundingBox) = (box.xmin, box.xmax, box.ymin, box.ymax)
Base.getindex(box::BoundingBox, i::Integer) =
    i == 1 ? box.xmin :
    i == 2 ? box.xmax :
    i == 3 ? box.ymin :
    i == 4 ? box.ymax : throw(BoundsError(box, i))
Base.iterate(box::BoundingBox, i::Int=1) =
    i == 1 ? (box.xmin, 2) :
    i == 2 ? (box.xmax, 3) :
    i == 3 ? (box.ymin, 4) :
    i == 4 ? (box.ymax, 5) : nothing

for Class in (:Point, :BoundingBox)
    @eval Base.map(f, obj::$Class) = $Class(map(f, Tuple(obj)))
end

# Constructors of points and conversion to/from a Cartesian index.
Point(pnt::Point) = pnt
Point{T}(pnt::Point{T}) where {T} = pnt
Point(pnt::AbstractPoint) = Point(pnt.x, pnt.y)
Point{T}(pnt::AbstractPoint) where {T} = Point{T}(pnt.x, pnt.y)
Point(x::Tx, y::Ty) where {Tx<:Real,Ty<:Real} = Point{promote_type(Tx,Ty)}(x, y)
Point(; x::Real, y::Real) = Point(x, y)
Point{T}(; x::Real, y::Real) where {T} = Point{T}(x, y)
Point(I::CartesianIndex{2}) = Point(I[1], I[2])
Point{T}(I::CartesianIndex{2}) where {T} = Point{T}(I[1], I[2])
Point((x,y)::NTuple{2,Real}) = Point(x, y)
Point{T}((x,y)::NTuple{2,Real}) where {T} = Point{T}(x, y)

# Conversion to points (rely on constructors).
Base.convert(::Type{T}, arg::T) where {T<:Point} = arg
Base.convert(::Type{T}, arg::PointLike) where {T<:Point} = T(arg)

# Other basic methods.
Base.CartesianIndex(pnt::Point{<:Integer}) = CartesianIndex(pnt.x, pnt.y)
Base.eltype(::AbstractPoint{T}) where {T} = T
Broadcast.broadcasted(::Type{T}, pnt::Point) where {T<:Real} = Point{T}(pnt)
Base.promote_type(::Type{Point{T}}, ::Type{Point{U}}) where {T,U} = Point{promote_type(T,U)}
Base.promote_type(::Type{Point{T}}, ::Type{Point{T}}) where {T} = Point{T}

for func in (:one, :oneunit, :zero)
    @eval Base.$func(pnt::Point) = $func(typeof(pnt))
end
Base.zero(::Type{Point{T}}) where {T} = Point(zero(T),zero(T))
Base.one(::Type{Point{T}}) where {T} = Point(one(T),one(T))
Base.oneunit(::Type{Point{T}}) where {T} = Point(oneunit(T),oneunit(T))

# Constructors of weighted points.
WeightedPoint(pnt::WeightedPoint) = pnt
WeightedPoint{T}(pnt::WeightedPoint{T}) where {T<:AbstractFloat} = pnt
WeightedPoint{T}(pnt::WeightedPoint) where {T} = WeightedPoint{T}(pnt.w, pnt.x, pnt.y)
WeightedPoint(w::Tw, x::Tx, y::Ty) where {Tw<:Real,Tx<:Real,Ty<:Real} =
    WeightedPoint{float(promote_type(Tw,Tx,Ty))}(w, x, y)
WeightedPoint(; w::Real, x::Real, y::Real) = WeightedPoint(w, x, y)
WeightedPoint(pnt::Point{T}) where {T} = WeightedPoint(one(T), pnt.x, pnt.y)
WeightedPoint((w,x,y)::NTuple{3,Real}) = WeightedPoint(w, x, y)
WeightedPoint{T}((w,x,y)::NTuple{3,Real}) where {T} = WeightedPoint{T}(w, x, y)
Broadcast.broadcasted(::Type{T}, pnt::WeightedPoint) where {T<:AbstractFloat} =
    WeightedPoint{T}(pnt)
Base.promote_type(::Type{WeightedPoint{T}}, ::Type{WeightedPoint{U}}) where {T,U} =
    WeightedPoint{promote_type(T,U)}
Base.promote_type(::Type{WeightedPoint{T}}, ::Type{WeightedPoint{T}}) where {T} =
    WeightedPoint{T}

"""
    round([T,] obj::Union{Point,BoundingBox}, [r::RoundingMode])

yields the object that is the nearest to `obj` by rounding its coordinates to
nearest integral values. Argument `T` can be the type of the result (a point or
a bounding-box) or the type of the coordinates of the result. Rounding mode may
be specified by optional argument `r`, the default being the same as the
`round` method for a scalar value.

For points, see also: [`floor(::Point)`](@ref), [`ceil(::Point)`](@ref).

For bounding-boxes, see also: [`interior`](@ref), [`exterior`](@ref).

""" Base.round

# NOTE: Using the following functor is a bit faster than map with an anonymous
#       function.
struct Round{T,R} <: Function
    r::R
    Round{T}(r::R) where {T,R} = new{T,R}(r)
end
Round{T}() where {T} = Round{T}(nothing)
(f::Round{Nothing,Nothing})(x) = round(x)
(f::Round{Nothing,<:RoundingMode})(x) = round(x, f.r)
(f::Round{T,Nothing})(x) where {T<:Number} = round(T, x)
(f::Round{T,<:RoundingMode})(x) where {T<:Number} = round(T, x, f.r)

for Class in (:Point, :BoundingBox)
    @eval begin
        Base.round(obj::$Class) = map(round, obj)
        Base.round(obj::$Class, r::RoundingMode) = map(Round{Nothing}(r), obj)
        Base.round(::Type{T}, obj::$Class) where {T<:Number} = map(Round{T}(), obj)
        Base.round(::Type{T}, obj::$Class, r::RoundingMode) where {T<:Number} = map(Round{T}(r), obj)

        Base.round(::Type{$Class}, obj::$Class) = round(obj)
        Base.round(::Type{$Class}, obj::$Class, r::RoundingMode) = round(obj, r)
        Base.round(::Type{$Class{T}}, obj::$Class) where {T} = round(T, obj)
        Base.round(::Type{$Class{T}}, obj::$Class, r::RoundingMode) where {T} = round(T, obj, r)
    end
end

"""
    floor([T,] P::Point)

yields the point with the largest integer coordinates smaller or equal those of
the point `P`.  Argument `T` can be the type of the result or the type of the
coordinates of the result.

See also: [`round`](@ref), [`ceil`](@ref).

""" Base.floor

"""
    ceil([T,] P::Point)

yields the point with the smallest integer coordinates larger or equal those of
the point `P`.  Argument `T` can be the type of the result or the type of the
coordinates of the result.

See also: [`round`](@ref), [`floor`](@ref).

""" Base.ceil

for Class in (:Point, :BoundingBox), func in (:floor, :ceil)
    @eval begin
        Base.$func(obj::$Class) = map($func, obj)
        Base.$func(::Type{T}, obj::$Class) where {T<:Number} = map(x -> $func(T, x), obj)

        Base.$func(::Type{$Class}, obj::$Class) = $func(obj)
        Base.$func(::Type{$Class{T}}, obj::$Class) where {T} = $func(T, obj)
    end
end

function Base.clamp(pnt::Point{T1}, box::BoundingBox{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return Point{T}(clamp(T(pnt.x), T(box.xmin), T(box.xmax)),
                    clamp(T(pnt.y), T(box.ymin), T(box.ymax)))
end

# Methods hypot() and atan() yield the polar coordinates of a point.
Base.hypot(pnt::Point) = hypot(pnt.x, pnt.y)
Base.atan(pnt::Point) = atan(pnt.y, pnt.x)

"""
    distance(A, B)

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
BoundingBox{T}(; xmin::Real, xmax::Real, ymin::Real, ymax::Real) where {T} =
    BoundingBox{T}(xmin, xmax, ymin, ymax)

BoundingBox(box::BoundingBox) = box
BoundingBox{T}(box::BoundingBox{T}) where {T<:Real} = box
BoundingBox{T}(box::BoundingBox) where {T<:Real} =
    BoundingBox{T}(box.xmin, box.xmax, box.ymin, box.ymax)

BoundingBox(arg::NTuple{2,CartesianIndex{2}}) = BoundingBox(arg...)
BoundingBox{T}(arg::NTuple{2,CartesianIndex{2}}) where {T} =
    BoundingBox{T}(arg...)

BoundingBox(arg::NTuple{2,AbstractPoint}) = BoundingBox(arg...)
BoundingBox{T}(arg::NTuple{2,AbstractPoint}) where {T} = BoundingBox{T}(arg...)

BoundingBox(I0::CartesianIndex{2}, I1::CartesianIndex{2}) =
    BoundingBox{Int}(I0, I1)
BoundingBox{T}(I0::CartesianIndex{2}, I1::CartesianIndex{2}) where {T} =
    BoundingBox{T}(I0[1], I1[1], I0[2], I1[2])

BoundingBox(pnt0::AbstractPoint, pnt1::AbstractPoint) =
    BoundingBox(pnt0.x, pnt1.x, pnt0.y, pnt1.y)
BoundingBox{T}(pnt0::AbstractPoint, pnt1::AbstractPoint) where {T} =
    BoundingBox{T}(pnt0.x, pnt1.x, pnt0.y, pnt1.y)

BoundingBox((x0,y0)::NTuple{2,Real}, (x1,y1)::NTuple{2,Real}) =
    BoundingBox(x0,x1,y0,y1)
BoundingBox{T}((x0,y0)::NTuple{2,Real}, (x1,y1)::NTuple{2,Real}) where {T} =
    BoundingBox{T}(x0,x1,y0,y1)

BoundingBox((xmin,xmax,ymin,ymax)::NTuple{4,Real}) =
    BoundingBox(xmin,xmax,ymin,ymax)
BoundingBox{T}((xmin,xmax,ymin,ymax)::NTuple{4,Real}) where {T} =
    BoundingBox{T}(xmin,xmax,ymin,ymax)

BoundingBox(inds::NTuple{2,AbstractUnitRange{<:Integer}}) =
    BoundingBox(inds[1], inds[2])
BoundingBox{T}(inds::NTuple{2,AbstractUnitRange{<:Integer}}) where {T} =
    BoundingBox{T}(inds[1], inds[2])

BoundingBox(X::AbstractUnitRange{<:Integer}, Y::AbstractUnitRange{<:Integer}) =
    BoundingBox{Int}(X, Y)
BoundingBox{T}(X::AbstractUnitRange{<:Integer}, Y::AbstractUnitRange{<:Integer}) where {T} =
    BoundingBox{T}(first(X), last(X), first(Y), last(Y))

BoundingBox(R::CartesianIndices{2}) = BoundingBox(first(R), last(R))
BoundingBox{T}(R::CartesianIndices{2}) where {T} = BoundingBox{T}(first(R), last(R))

@deprecate BoundingBox(A::AbstractMatrix) BoundingBox(axes(A))
BoundingBox(A::AbstractMatrix{Bool}) = BoundingBox(identity, A)

# Conversion to bounding-boxes (rely on constructors).
Base.convert(::Type{T}, arg::T) where {T<:BoundingBox} = arg
Base.convert(::Type{T}, arg::BoundingBoxLike) where {T<:BoundingBox} = T(arg)

# Other basic methods.
Base.eltype(::BoundingBox{T}) where {T} = T
Broadcast.broadcasted(::Type{T}, obj::BoundingBox) where {T<:Real} =
    BoundingBox{T}(obj)
Base.promote_type(::Type{BoundingBox{T}}, ::Type{BoundingBox{U}}) where {T,U} =
    BoundingBox{promote_type(T,U)}
Base.promote_type(::Type{BoundingBox{T}}, ::Type{BoundingBox{T}}) where {T} =
    BoundingBox{T}

# See
# https://stackoverflow.com/questions/9852159/calculate-bounding-box-of-arbitrary-pixel-based-drawing
# for the basic ideas under the following algorithm.
function BoundingBox(f::Function, A::AbstractMatrix)
    I, J = axes(A)
    i0, i1 = getaxisbounds(I)
    j0, j1 = getaxisbounds(J)
    imin = jmin = typemax(Int)
    imax = jmax = typemin(Int)
    # Assuming column-major order, first start by scanning rows to narrow the
    # subsequent searches along columns.
    #
    # 1. Find bottom bound `jmin` by scanning rows from bottom to top.
    flag = false
    @inbounds for j in j0:j1, i in i0:i1
        if f(A[i,j])
            # This definitively set the value of `jmin` and gives limits for
            # the other bounds.
            imin = imax = i
            jmin = jmax = j
            flag = true
            break
        end
    end
    if flag
        # 2. Find top bound `jmax` by scanning rows from top to bottom.  No
        #    needs to go beyond `jmax+1`.
        @inbounds for j in j1:-1:jmax+1, i in i0:i1
            if f(A[i,j])
                jmax = j
                imin = min(imin, i)
                imax = max(imax, i)
                break
            end
        end
        # 3. Find leftmost bound `imin` by scanning columns from left to right.
        #    No needs to go beyond `imin-1`.
        @inbounds for i in i0:imin-1, j in jmin:jmax
            if f(A[i,j])
                imin = i
                imax = max(imax, i)
                break
            end
        end
        # 4. Find rightmost bound `imax` by scanning columns from right to
        #    left.  No needs to go beyond `imax+1`.
        @inbounds for i in i1:-1:imax+1, j in jmin:jmax
            if f(A[i,j])
                imax = i
                break
            end
        end
    end
    return BoundingBox(imin, imax, jmin, jmax)
end

"""
    getaxisbounds(I) = (i0,i1)

yields the bounds `i0` and `i1` of index range `I` as a 2-tuple of `Int`'s and
such that `i0:i1` represents the same indices as `I` (although not in the same
order if `step(I) < 0`).  If `step(I)` is not equal to ±1, an `ArgumentError`
exception is thrown.

"""
@inline function getaxisbounds(I::AbstractRange{<:Integer})
    i0, i1, s = Int(first(I)), Int(last(I)), step(I)
    return (s == +1 ? (i0,i1) :
            s == -1 ? (i1,i0) :
            throw(ArgumentError("expecting a range with a step equal to ±1")))
end

# Empty bounding and unlimited boxes.
BoundingBox{T}(::Nothing) where {T<:Real} = typemin(BoundingBox{T})
Base.typemin(::Type{BoundingBox{T}}) where {T<:Real} =
    BoundingBox(typemax(T), typemin(T), typemax(T), typemin(T))
Base.typemax(::Type{BoundingBox{T}}) where {T<:Real} =
    BoundingBox(typemin(T), typemax(T), typemin(T), typemax(T))

# Conversion of bounding-boxes to/from Cartesian indices.
Base.CartesianIndices(box::BoundingBox{<:Integer}) =
    CartesianIndices((Int(box.xmin):Int(box.xmax), Int(box.ymin):Int(box.ymax)))

# Lower left and upper right corners of a bounding-box.
Base.first(box::BoundingBox) = Point(box.xmin, box.ymin)
Base.last(box::BoundingBox) = Point(box.xmax, box.ymax)

Base.isempty(box::BoundingBox) = (box.xmin > box.xmax)|(box.ymin > box.ymax)

Base.size(box::BoundingBox{Int}) =
    isempty(box) ? (0, 0) : (box.xmax - box.xmin + 1, box.ymax - box.ymin + 1)
Base.size(box::BoundingBox{Int}, k::Integer) =
    isempty(box) ? (k ≥ 1 ? 0 : throw_bad_dimension_index()) :
    k == 1 ? box.xmax - box.xmin + 1 :
    k == 2 ? box.ymax - box.ymin + 1 :
    k > 2 ? 1 : throw_bad_dimension_index()

Base.size(box::BoundingBox{<:Integer}) =
    isempty(box) ? (0, 0) :
    (Int(box.xmax) - Int(box.xmin) + 1, Int(box.ymax) - Int(box.ymin) + 1)

Base.size(box::BoundingBox{<:Integer}, k::Integer) =
    isempty(box) ? (k ≥ 1 ? 0 : throw_bad_dimension_index()) :
    k == 1 ? Int(box.xmax) - Int(box.xmin) + 1 :
    k == 2 ? Int(box.ymax) - Int(box.ymin) + 1 :
    k > 2 ? 1 : throw_bad_dimension_index()

Base.axes(box::BoundingBox{Int}) = (box.xmin:box.xmax, box.ymin:box.ymax)
Base.axes(box::BoundingBox{Int}, k::Integer) =
    k == 1 ? (box.xmin:box.xmax) :
    k == 2 ? (box.ymin:box.ymax) :
    k > 2 ? (1:1) : throw_bad_dimension_index()

Base.axes(box::BoundingBox{<:Integer}) = (Int(box.xmin):Int(box.xmax),
                                          Int(box.ymin):Int(box.ymax))
Base.axes(box::BoundingBox{<:Integer}, k::Integer) =
    k == 1 ? (Int(box.xmin):Int(box.xmax)) :
    k == 2 ? (Int(box.ymin):Int(box.ymax)) :
    k > 2 ? (1:1) : throw_bad_dimension_index()

@noinline throw_bad_dimension_index() =
    error("invalid dimension index")

Base.in(pnt::AbstractPoint, box::BoundingBox) =
    (box.xmin ≤ pnt.x ≤ box.xmax) & (box.ymin ≤ pnt.y ≤ box.ymax)

# To deprecate `A ∈ B` in favor of `A ⊆ B` for bounding boxes `A` and `B`, the
# `Base.in` method must be imported.
import Base: in
@deprecate in(A::BoundingBox, B::BoundingBox) (A ⊆ B) false

# Extend ⊆ operator.
Base.issubset(A::BoundingBox, B::BoundingBox) =
    (isempty(A)|((A.xmin ≥ B.xmin)&(A.xmax ≤ B.xmax)&
                 (A.ymin ≥ B.ymin)&(A.ymax ≤ B.ymax)))
# Union of bounding-boxes:
Base.:(∪)(A::BoundingBox, B::BoundingBox) =
    BoundingBox(min(A.xmin, B.xmin), max(A.xmax, B.xmax),
                min(A.ymin, B.ymin), max(A.ymax, B.ymax))

# Intersection of bounding-boxes:
Base.:(∩)(A::BoundingBox, B::BoundingBox) =
    BoundingBox(max(A.xmin, B.xmin), min(A.xmax, B.xmax),
                max(A.ymin, B.ymin), min(A.ymax, B.ymax))

# Use bounding-boxes to extract a sub-array or a view.
@propagate_inbounds function Base.getindex(A::AbstractMatrix,
                                           B::BoundingBox{<:Integer})
    A[B.xmin:B.xmax, B.ymin:B.ymax]
end

Base.view(A::AbstractMatrix, B::BoundingBox{<:Integer}) =
    view(A, B.xmin:B.xmax, B.ymin:B.ymax)

"""
    interior([T,] box)

yields the largest bounding-box with integer valued bounds and which is
contained by the bounding-box `box`. Optional argument `T` is to specify the
type of the result or of the coordinates of the result which is the same as
`box` by default.

See also: [`exterior`](@ref), [`round`](@ref).

"""
interior(box::BoundingBox{T}) where {T} = interior(T, box)
interior(::Type{BoundingBox{T}}, box::BoundingBox) where {T} = interior(T, box)
interior(::Type{T}, box::BoundingBox{T}) where {T<:Integer} = box
interior(::Type{T}, box::BoundingBox{T}) where {T<:Real} =
    BoundingBox(ceil(box.xmin), floor(box.xmax),
                ceil(box.ymin), floor(box.ymax))
interior(::Type{T}, box::BoundingBox{<:Integer}) where {T<:Integer} =
    BoundingBox{T}(box)
interior(::Type{T}, box::BoundingBox{<:Real}) where {T<:Integer} =
    BoundingBox(ceil(T, box.xmin), floor(T, box.xmax),
                ceil(T, box.ymin), floor(T, box.ymax))
interior(::Type{T}, box::BoundingBox{<:Integer}) where {T<:Real} =
    BoundingBox{T}(box)
interior(::Type{T}, box::BoundingBox{U}) where {T<:Real,U<:Real} =
    BoundingBox{T}(interior(U, box))

"""
    exterior([T,] box)

yields the smallest bounding-box with integer valued bounds and which contains
the bounding-box `box`.  Optional argument `T` is to specify the type of the
result or of the coordinates of the result which is the same as `box` by default.

See also: [`interior`](@ref), [`round`](@ref).

"""
exterior(box::BoundingBox{T}) where {T} = exterior(T, box)
exterior(::Type{BoundingBox{T}}, box::BoundingBox) where {T} = exterior(T, box)
exterior(::Type{T}, box::BoundingBox{T}) where {T<:Integer} = box
exterior(::Type{T}, box::BoundingBox{T}) where {T<:Real} =
    BoundingBox(floor(box.xmin), ceil(box.xmax),
                floor(box.ymin), ceil(box.ymax))
exterior(::Type{T}, box::BoundingBox{<:Integer}) where {T<:Integer} =
    BoundingBox{T}(box)
exterior(::Type{T}, box::BoundingBox{<:Real}) where {T<:Integer} =
    BoundingBox(floor(T, box.xmin), ceil(T, box.xmax),
                floor(T, box.ymin), ceil(T, box.ymax))
exterior(::Type{T}, box::BoundingBox{<:Integer}) where {T<:Real} =
    BoundingBox{T}(box)
exterior(::Type{T}, box::BoundingBox{U}) where {T<:Real,U<:Real} =
    BoundingBox{T}(exterior(U, box))

"""
    center(box::BoundingBox) -> c::Point

yields the central point of the bounding-box `box`.

"""
center(box::BoundingBox{<:Integer}) =
    Point(0.5*(box.xmin + box.xmax), 0.5*(box.ymin + box.ymax))
center(box::BoundingBox{T}) where {T<:AbstractFloat} =
    Point(half(T)*(box.xmin + box.xmax), half(T)*(box.ymin + box.ymax))

half(::Type{T}) where {T<:AbstractFloat} = one(T)/convert(T, 2)

"""
    area(box)

yields the area of the bounding-box `box`.

"""
area(box::BoundingBox{T}) where {T<:Real} =
    max(box.xmax - box.xmin, zero(T))*max(box.ymax - box.ymin, zero(T))

# Scaling of a point coordinates.
Base.:(*)(pnt::Point, α::Real) = α*pnt
Base.:(*)(α::Real, pnt::Point) = Point(α*pnt.x, α*pnt.y)
Base.:(/)(pnt::Point, α::Real) = Point(pnt.x/α, pnt.y/α)
Base.:(\)(α::Real, pnt::Point) = pnt/α

# Unary plus does nothing.
Base.:(+)(pnt::Point) = pnt
Base.:(+)(box::BoundingBox) = box

# Unary minus applied to a point negate coordinates.
Base.:(-)(pnt::Point{<:Union{AbstractFloat,Signed,Irrational}}) =
    Point(-pnt.x, -pnt.y)

# Unary minus applied to a bounding-box.
Base.:(-)(box::BoundingBox{<:Union{AbstractFloat,Signed,Irrational}}) =
    BoundingBox(-box.xmax, -box.xmin, -box.ymax, -box.ymin)

# Addition and subtraction of point coordinates.
Base.:(+)(A::Point, B::Point) = Point(A.x + B.x, A.y + B.y)
Base.:(-)(A::Point, B::Point) = Point(A.x - B.x, A.y - B.y)

# Scaling of bounding-box bounds (e.g. to change units).
Base.:(*)(box::BoundingBox, α::Real) = α*box
Base.:(*)(α::Real, box::BoundingBox) = BoundingBox(α*box.xmin, α*box.xmax,
                                                   α*box.ymin, α*box.ymax)
Base.:(/)(box::BoundingBox{T}, α::Real) where {T} = (one(T)/α)*box
Base.:(\)(α::Real, box::BoundingBox) = box/α

# Add or remove a margin δ to a bounding-box box.
Base.:(+)(box::BoundingBox, δ::Real) =
    BoundingBox(box.xmin - δ, box.xmax + δ,
                box.ymin - δ, box.ymax + δ)
Base.:(-)(box::BoundingBox, δ::Real) =
    BoundingBox(box.xmin + δ, box.xmax - δ,
                box.ymin + δ, box.ymax - δ)

# Translate a bounding-box.
Base.:(+)(box::BoundingBox, pnt::Point) =
    BoundingBox(box.xmin + pnt.x, box.xmax + pnt.x,
                box.ymin + pnt.y, box.ymax + pnt.y)
Base.:(-)(box::BoundingBox, pnt::Point) =
    BoundingBox(box.xmin - pnt.x, box.xmax - pnt.x,
                box.ymin - pnt.y, box.ymax - pnt.y)

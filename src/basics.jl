#
# basics.jl --
#
# Basic methods points and bounding-boxes.
#
#-------------------------------------------------------------------------------
#
# This file is part of the TwoDimensional Julia package licensed under the MIT
# license (https://github.com/emmt/TwoDimensional.jl).
#
# Copyright (c) 2019-2024, Éric Thiébaut.
#

# FIXME: use TypeUtils
as(::Type{T}, x) where {T} = convert(T, x)::T

# Basic constructors (with tuple argument) for points and bounding boxes and
# API to make them indexable iterators.
for (type, len) in ((:Point,         2),
                    (:WeightedPoint, 3),
                    (:BoundingBox,   4),)
    @eval begin
        Base.Tuple(obj::$type) = getfield(obj, 1)
        Base.getindex(obj::$type, i::Integer) = getindex(Tuple(obj), i)
        Base.iterate(obj::$type, i::Int = 1) =
            1 ≤ i ≤ $len ? (obj[i], i + 1) : nothing
        Base.IteratorSize(obj::$type) = Base.IteratorSize(typeof(obj))
        Base.IteratorSize(::Type{<:$type}) = Base.HasLength()
        Base.length(obj::$type) = $len
        Base.show(io::IO, ::MIME"text/plain", obj::$type) = show(io, obj)
    end
    if type === :WeightedPoint
        @eval begin
            $type(vals::NTuple{$len,T}) where {T<:AbstractFloat} = $type{T}(vals)
            $type(vals::NTuple{$len,Real}) = $type(map(float, promote(vals...)))
        end
    else
        @eval begin
            $type(vals::NTuple{$len,T}) where {T} = $type{T}(vals)
            $type(vals::NTuple{$len,Any}) = $type(promote(vals...))
            Base.map(f, obj::$type) = $type(map(f, Tuple(obj)))
        end
    end
end

Base.propertynames(::Point) = (:x, :y)
Base.getproperty(pnt::Point, key::Symbol) =
    key === :x ? pnt[1] :
    key === :y ? pnt[2] : throw(KeyError(key))
function Base.show(io::IO, pnt::Point{T}) where {T}
    print(io, "Point{")
    show(io, T)
    print(io, "}(x = "); show(io, pnt.x)
    print(io, ", y = "); show(io, pnt.y)
    print(io, ")")
end

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

Base.propertynames(::BoundingBox) = (:xmin, :xmax, :ymin, :ymax)
Base.getproperty(box::BoundingBox, key::Symbol) =
    key === :xmin ? box[1] :
    key === :xmax ? box[2] :
    key === :ymin ? box[3] :
    key === :ymax ? box[4] : throw(KeyError(key))
function Base.show(io::IO, box::BoundingBox{T}) where {T}
    print(io, "BoundingBox{")
    show(io, T)
    print(io, "}(xmin = "); show(io, box.xmin)
    print(io, ", xmax = "); show(io, box.xmax)
    print(io, ", ymin = "); show(io, box.ymin)
    print(io, ", ymax = "); show(io, box.ymax)
    print(io, ")")
end

# Other constructors of points.
Point{T}(x, y) where {T} = Point{T}((x, y))
Point(x, y) = Point(promote(x, y))
Point(x::T, y::T) where {T} = Point{T}(x, y)

Point(; x, y) = Point(x, y)
Point{T}(; x, y) where {T} = Point{T}(x, y)

Point(pnt::Point) = pnt
Point{T}(pnt::Point{T}) where {T} = pnt
Point{T}(pnt::Point) where {T} = Point{T}(Tuple(pnt))

Point(pnt::PointLike) = Point(get_xy(pnt))
Point{T}(pnt::PointLike) where {T} = Point{T}(get_xy(pnt))

# Extend some basic methods for points.
for func in (:one, :oneunit, :zero)
    @eval Base.$func(pnt::Point) = $func(typeof(pnt))
end
Base.zero(::Type{Point{T}}) where {T} = Point(zero(T),zero(T))
Base.one(::Type{Point{T}}) where {T} = Point(one(T),one(T)) # FIXME:
Base.oneunit(::Type{Point{T}}) where {T} = Point(oneunit(T),oneunit(T)) # FIXME:

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

# Extend basic methods for abstract points.
Base.eltype(::AbstractPoint{T}) where {T} = T
Base.CartesianIndex(pnt::AbstractPoint{<:Integer}) = CartesianIndex(get_xy(pnt)...)

# Extend basic methods for points and bounding-boxes.
for (type, like, n) in ((:Point,         :PointLike,         2),
                        (:WeightedPoint, :WeightedPointLike, 3),
                        (:BoundingBox,   :BoundingBoxLike,   4),)
    @eval begin
        Base.convert(::Type{T}, obj::T) where {T<:$type} = obj
        Base.convert(::Type{T}, obj::$like) where {T<:$type} = T(obj)
        Base.convert(::Type{Tuple}, obj::$type) = Tuple(obj)
        Base.convert(::Type{NTuple{$n,T}}, obj::$type) where {T} =
            map(#= FIXME: as(T) =# Base.Fix1(as, T), Tuple(obj))
        Base.eltype(obj::$type) = eltype(typeof(obj))
        Base.eltype(::Type{<:$type{T}}) where {T} = T
        Broadcast.broadcasted(::Type{T}, obj::$type) where {T} = $type{T}(obj)
        Base.promote_type(::Type{$type{T}}, ::Type{$type{T}}) where {T} = $type{T}
        Base.promote_type(::Type{$type{T₁}}, ::Type{$type{T₂}}) where {T₁,T₂} =
            $type{promote_type(T₁,T₂)}
    end
end

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

function Base.clamp(pnt::Point, box::BoundingBox)
    T = promote_type(eltype(pnt), eltype(box))
    x = clamp(as(T, pnt.x), as(T, box.xmin), as(T, box.xmax))
    y = clamp(as(T, pnt.y), as(T, box.ymin), as(T, box.ymax))
    return Point{T}(x, y)
end

# Methods hypot() and atan() yield the polar coordinates of a point.
Base.hypot(pnt::Point) = hypot(pnt.x, pnt.y)
Base.atan(pnt::Point) = atan(pnt.y, pnt.x)

"""
    distance(A, B)

yields the Euclidean distance between the 2 points `A` and `B`.

"""
distance(A::Point, B::Point) = hypot(A - B)

# Box limits specified by 4 arguments.
BoundingBox{T}(xmin, xmax, ymin, ymax) where {T} = BoundingBox{T}((xmin, xmax, ymin, ymax))
BoundingBox(xmin, xmax, ymin, ymax) = BoundingBox(promote(xmin, xmax, ymin, ymax))
BoundingBox(xmin::T, xmax::T, ymin::T, ymax::T) where {T} =
    BoundingBox{T}(xmin, xmax, ymin, ymax)

# Box limits specified by keywords.
BoundingBox(; xmin, xmax, ymin, ymax) = BoundingBox(xmin, xmax, ymin, ymax)
BoundingBox{T}(; xmin, xmax, ymin, ymax) where {T} = BoundingBox{T}(xmin, xmax, ymin, ymax)

# Box specified by a ... box.
BoundingBox(box::BoundingBox) = box
BoundingBox{T}(box::BoundingBox{T}) where {T} = box
BoundingBox{T}(box::BoundingBox) where {T} = BoundingBox{T}(Tuple(box))

# Box limits specified as unit-ranges.
let type = :(AbstractUnitRange{<:Integer}),
    expr = (:(first(xrng)), :(last(xrng)), :(first(yrng)), :(last(yrng)))
    for args in ((:(xrng::$type), :(yrng::$type)),
                 (:((xrng,yrng)::Tuple{$type,$type}),))
        @eval begin
            BoundingBox($(args...)) = BoundingBox($(expr...))
            BoundingBox{T}($(args...)) where {T} = BoundingBox{T}($(expr...))
        end
    end
end

# Box limits specified by a pair of points, of Cartesian indices, or of 2-tuple.
let expr = (:(get_x(min)), :(get_x(max)), :(get_y(min)), :(get_y(max)))
    for type in (:Point, :(CartesianIndex{2}), :(Tuple{Any,Any}))
        for args in ((:(min::$type), :(max::$type)),
                     (:((min,max)::Tuple{$type,$type}),))
            @eval begin
                BoundingBox($(args...)) = BoundingBox($(expr...))
                BoundingBox{T}($(args...)) where {T} = BoundingBox{T}($(expr...))
            end
        end
    end
end

# Box limits specified by Cartesian indices.
BoundingBox(inds::CartesianIndices{2}) = BoundingBox(first(inds), last(inds))
BoundingBox{T}(inds::CartesianIndices{2}) where {T} = BoundingBox{T}(first(inds), last(inds))

# See
# https://stackoverflow.com/questions/9852159/calculate-bounding-box-of-arbitrary-pixel-based-drawing
# for the basic ideas under the following algorithm.
BoundingBox(A::AbstractMatrix{Bool}) = BoundingBox(identity, A)
function BoundingBox(f::Function, A::AbstractMatrix)
    I, J = axes(A)
    i0, i1 = get_axis_bounds(I)
    j0, j1 = get_axis_bounds(J)
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
    TwoDimensional.get_axis_bounds(I) = (i0,i1)

yields the bounds `i0` and `i1` of index range `I` as a 2-tuple of `Int`'s and
such that `i0:i1` represents the same indices as `I` (although not in the same
order if `step(I) < 0`). If `step(I)` is not equal to ±1, an `ArgumentError`
exception is thrown.

"""
@inline function get_axis_bounds(I::AbstractRange{<:Integer})
    i0, i1, s = Int(first(I)), Int(last(I)), step(I)
    return (s == +oneunit(s) ? (i0,i1) :
            s == -oneunit(s) ? (i1,i0) :
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

Base.size(box::BoundingBox{<:Integer}) = map(length, axes(box))
Base.size(box::BoundingBox{<:Integer}, d::Integer) = length(axes(box, d))

Base.axes(box::BoundingBox{<:Integer}) = (UnitRange{Int}(box.xmin, box.xmax),
                                          UnitRange{Int}(box.ymin, box.ymax))
Base.axes(box::BoundingBox{<:Integer}, d::Integer) =
    d == 1 ? UnitRange{Int}(box.xmin, box.xmax) :
    d == 2 ? UnitRange{Int}(box.ymin, box.ymax) :
    d > 2 ? (1:1) : throw_bad_dimension_index()

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

half(::Type{T}) where {T<:AbstractFloat} = inv(as(T, 2))

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

"""
    TwoDimensional.get_x(pnt::TwoDimensional.PointLike) -> x

yields the abscissa of point-like object `pnt`.

See also [`TwoDimensional.get_y`](@ref), [`TwoDimensional.get_xy`](@ref), and
[`TwoDimensional.PointLike`](@ref).

""" get_x

"""
    TwoDimensional.get_y(pnt::TwoDimensional.PointLike) -> y

yields the ordinate of point-like object `pnt`.

See also [`TwoDimensional.get_x`](@ref), [`TwoDimensional.get_xy`](@ref), and
[`TwoDimensional.PointLike`](@ref).

""" get_y

for (c, i) in ((:x, 1), (:y, 2))
    func = Symbol("get_",c)
    @eval begin
        $(func)(pnt::Point) = pnt[$(i)]
        $(func)(pnt::Union{NTuple{2},CartesianIndex{2}}) = pnt[$(i)]
        $(func)(pnt::AbstractPoint) = pnt.$(c)
        $(func)(pnt::PointLike) = get_xy(pnt)[$(i)]
    end
end

"""
    TwoDimensional.get_xy(pnt::TwoDimensional.PointLike) -> (x::T, y::T)

yields a 2-tuple with the abscissa `x` and ordinate `y` of point-like object
`pnt`. This is equivalent to, but more economical than,
`(get_x(pnt),get_y(pnt))`.

See also [`TwoDimensional.get_x`](@ref), [`TwoDimensional.get_y`](@ref) and [`TwoDimensional.PointLike`](@ref).

"""
get_xy(pnt::Point) = Tuple(pnt)
get_xy(pnt::CartesianIndex{2}) = Tuple(pnt)
get_xy(pnt::AbstractPoint) = (pnt.x, pnt.y)
get_xy(pnt::NTuple{2}) = pnt
get_xy(pnt::NTuple{2,Any}) = promote(pnt...)

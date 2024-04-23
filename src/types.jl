#
# types.jl --
#
# Type definitions and constants for TwoDimensinal package.
#
#-------------------------------------------------------------------------------
#
# This file is part of the TwoDimensional Julia package licensed under the MIT
# license (https://github.com/emmt/TwoDimensional.jl).
#
# Copyright (c) 2019-2024, Éric Thiébaut.
#

"""
    TwoDimensional.GeometricObject{T}

is the super-type of geometric objects whose coordinates are of type `T`. Most
computations assume that `T` is floating-point (possibly with units).
Coordinate type `T` of a geometric object `obj` can be retrieved with
[`TwoDimensional.coord_type(obj)`](@ref `TwoDimensional.coord_type).

"""
abstract type GeometricObject{T} end

"""
    TwoDimensional.GeometricElement{T} <: TwoDimensional.GeometricObject{T}

is the super-type of elementary geometric objects whose coordinates are of type
`T`. Concrete objects of this type are the building blocks or more complex
geometric objects.

"""
abstract type GeometricElement{T} <: GeometricObject{T} end

"""
    TwoDimensional.ShapeElement{T} <: TwoDimensional.GeometricElement{T}

is the super-type of elementary geometric objects whose coordinates are of type
`T` and which can be used to specify boundaries.

"""
abstract type ShapeElement{T} <: GeometricElement{T} end

"""
    TwoDimensional.MaskElement(shape::ShapeElement, opaque::Bool)

builds a simple mask whose boundaries are defined by `shape` and which is
opaque (i.e., an obscuration) if `opaque` is true and transparent (i.e., an
aperture) otherwise.

"""
struct MaskElement{T,S} <: ShapeElement{T}
    shape::S
    opaque::Bool
    MaskElement(obj::S, opaque::Bool) where {T,S<:ShapeElement{T}} = new{T,S}(obj, opaque)
end

"""
    AbstractPoint{T} <: TwoDimensional.GeometricElement{T}

is the abstract type of objects with at least 2 properties: `x` and `y`, their
respective abscissa and ordinate, both of type `T`.

See also [`Point`](@ref).

"""
abstract type AbstractPoint{T} <: GeometricElement{T} end

struct Point{T} <: AbstractPoint{T}
    vec::NTuple{2,T} # x, y
    # Inner constructor: Point{T}(x::T, y::T)
    Point{T}(xy::Vararg{T,2}) where {T} = new{T}(xy)
end

struct Circle{T} <: ShapeElement{T}
    center::Point{T}
    radius::T
    function Circle{T}(center::Point, radius) where {T}
        radius ≥ zero(radius) || throw(ArgumentError("circle radius must be non-negative"))
        return new{T}(center, radius)
    end
end

struct Rectangle{T} <: ShapeElement{T}
    vec::NTuple{2,Point{T}} # Point(x0, y0), Point(x1, y1)
    function Rectangle{T}(start::Point{T}, stop::Point{T}) where {T}
        # The coordinates are ordered.
        x0, x1 = minmax(start.x, stop.x)
        y0, y1 = minmax(start.y, stop.y)
        return new{T}((Point{T}(x0, y0), Point{T}(x1, y1)))
    end
end

struct Polygon{T,V<:AbstractVector{Point{T}}} <: ShapeElement{T}
    vertices::V
    function Polygon{T}(vertices::V) where {T,V<:AbstractVector{<:Point{T}}}
        len = length(vertices)
        len ≥ 3 || throw_insufficent_number_of_polygon_vertices(len)
        isconcretetype(T) || throw(ArgumentError("coordinate type must be concrete"))
        return new{T,V}(vertices)
    end
end

struct BoundingBox{T} <: GeometricElement{T}
    vec::NTuple{2,Point{T}} # Point(xmin, ymin), Point(xmax, ymax)
    BoundingBox{T}(startstop::Vararg{Point{T},2}) where {T} = new{T}(startstop)
end

const RectangularObject{T} = Union{BoundingBox{T},Rectangle{T}}

"""
    TwoDimensional.VertexBasedObject{T}

is the union of types of objects defined by their verices and with coordinate
type `T`.

See also [`TwoDimensional.apply`].

"""
const VertexBasedObject{T} = Union{Point{T},
                                   Rectangle{T},
                                   Polygon{T},
                                   BoundingBox{T}}

const VERTEX_BASED_TYPES = (:Point, :Rectangle, :Polygon, :BoundingBox)

const TupleOrVector{T} = Union{Tuple{Vararg{T}},AbstractVector{<:T}}

"""
    TwoDimensional.PointLike

is the union of types that can be converted into a
[`TwoDimensional.Point`](@ref).

The [`Point`](@ref) constructor can build an instance from any argument of
these types. Accessors [`TwoDimensional.get_x`](@ref) and
[`TwoDimensional.get_y`](@ref) may be used on objects of such type to retrieve
their abscissa and ordinate.

"""
const PointLike = Union{AbstractPoint,
                        NTuple{2,Number}, # NOTE: Purposely more restrictive
                                          #       than when calling explicitly
                                          #       the constructor to avoid
                                          #       ambiguities.
                        CartesianIndex{2}}

"""
    TwoDimensional.RectangleLike

is the union of types that may be used to specify a rectangle in
`TwoDimensional` package.

The [`Rectangle`](@ref) constructor can build an instance from any argument of
these types.

"""
const RectangleLike = Union{Rectangle,
                            NTuple{2,AbstractPoint},
                            NTuple{2,NTuple{2,Number}},
                            NTuple{2,CartesianIndex{2}}}

"""
    TwoDimensional.CircleLike

is the union of types that may be used to specify a circle in
`TwoDimensional` package.

The [`Circle`](@ref) constructor can build an instance from any argument of
these types.

"""
const CircleLike = Union{Circle,
                         Tuple{PointLike,Number}}

"""
    TwoDimensional.PolygonLike

is the union of types that can be converted into a
[`TwoDimensional.Polygon`](@ref). These are tuples or vectors of point-like
objects.

The [`Polygon`](@ref) constructor can build an instance from any argument of
these types. Accessors [`TwoDimensional.get_x`](@ref) and
[`TwoDimensional.get_y`](@ref) may be used on objects of such type to retrieve
their abscissa and ordinate.

"""
const PolygonLike = Union{TupleOrVector{<:AbstractPoint},
                          TupleOrVector{<:NTuple{2,Number}},
                          TupleOrVector{<:CartesianIndex{2}}}

"""
    TwoDimensional.BoundingBoxLike

is the union of types that may be used to specify a bounding-box in
`TwoDimensional` package.

The [`BoundingBox`](@ref) constructor can build an instance from any argument
of these types.

"""
const BoundingBoxLike = Union{BoundingBox,
                              NTuple{2,AbstractPoint},
                              NTuple{2,NTuple{2,Number}},
                              NTuple{2,CartesianIndex{2}},
                              NTuple{2,AbstractUnitRange{<:Integer}},
                              CartesianIndices{2}}

struct AffineTransform{T<:AbstractFloat,R,S} <: Function
    factors::NTuple{4,R}
    offsets::NTuple{2,S}
    function AffineTransform{T,R,S}(Axx, Axy, Ax, Ayx, Ayy, Ay) where {T<:AbstractFloat,R,S}
        isconcretetype(T) || throw(ArgumentError(
            "type parameter `T = $T` is not a concrete floating-point type"))
        bare_type(R) === T || throw(ArgumentError(
            "bare type of parameter `R = $R` is not `T = $T`, got `bare_type(R) = $(bare_type(R))`"))
        bare_type(S) === T || throw(ArgumentError(
            "bare type of parameter `S = $S` is not `T = $T`, got `bare_type(S) = $(bare_type(S))`"))
        return new{T,R,S}((Axx, Axy, Ayx, Ayy), (Ax, Ay))
    end
end

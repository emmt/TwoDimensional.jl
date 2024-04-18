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
    AbstractPoint{T} <: TwoDimensional.GeometricElement{T}

is the abstract type of objects with at least 2 properties: `x` and `y`, their
respective abscissa and ordinate, both of type `T`.

See also [`Point`](@ref).

"""
abstract type AbstractPoint{T} <: GeometricElement{T} end

struct Point{T} <: AbstractPoint{T}
    vals::NTuple{2,T} # x, y
    Point{T}(vals::NTuple{2,Any}) where {T} = new{T}(vals)
end

"""
    TwoDimensional.PointLike

is the union of types that may be used to specify a point in `TwoDimensional`
package.

[`Point`](@ref) constructor can build an instance from any argument of these
types. Accessors [`TwoDimensional.get_x`](@ref) and
[`TwoDimensional.get_y`](@ref) may be used on objects of such type to retrieve
their abscissa and ordinate.

"""
const PointLike = Union{AbstractPoint,Tuple{Any,Any},CartesianIndex{2}}

struct BoundingBox{T} <: GeometricElement{T}
    vals::NTuple{4,T} # xmin, xmax, ymin, ymax
    BoundingBox{T}(vals::NTuple{4,Any}) where {T} = new{T}(vals)
end

"""
    TwoDimensional.BoundingBoxLike

is the union of types that may be used to specify a bounding-box in
`TwoDimensional` package.

[`BoundingBox`](@ref) constructor can build an instance from any argument of
these types.

"""
const BoundingBoxLike = Union{BoundingBox,NTuple{2,AbstractPoint},
                              Tuple{Tuple{Any,Any},Tuple{Any,Any}},
                              Tuple{Any,Any,Any,Any},
                              NTuple{2,CartesianIndex{2}},
                              NTuple{2,AbstractUnitRange{<:Integer}},
                              CartesianIndices{2}}

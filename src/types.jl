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
    AbstractPoint{T}

is the abstract type of objects with at least 2 properties: `x` and `y`, their
respective abscissa and ordinate, both of type `T`.

See also [`Point`](@ref).

"""
abstract type AbstractPoint{T} <: GeometricElement{T} end

"""
    Point(x,y)
    Point((x,y))

yield an instance of a 2D point of coordinates `(x,y)`.

A point may be multiplied or divided by a scalar to scale its coordinates. The
addition (resp. subtraction) of two points adds (resp. subtracts) their
coordinates.

Coordinates can be specified by keywords:

    Point(x=xval, y=yval)

There are no default values for keywords `x` and `y` so both must be specified.

The coordinates of a `Point`, say `pnt`, can be retrieved as follows:

    pnt.x  or  pnt[1]  ->  x
    pnt.y  or  pnt[2]  ->  y

or:

    x, y = pnt

See also [`AbstractPoint`](@ref).

"""
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

"""
    BoundingBox(xmin,xmax,ymin,ymax)
    BoundingBox((xmin,ymin),(xmax,ymax))

yield an instance of a 2D rectangular bounding-box whose sides are aligned with
the coordinate axes and containing points of coordinates `(x,y)` such that
`xmin ≤ x ≤ xmax` and `ymin ≤ y ≤ ymax`. The box is *empty* if `xmin > xmax` or
`ymin > ymax`.

A bounding-box can be constructed from the first and last points (i.e. at the
lower-left and upper right opposite corners) of the box:

    BoundingBox(P0::Point, P1::Point)
    BoundingBox(I0::CartesianIndex{2}, I1::CartesianIndex{2})

Coordinates can be specified by keywords:

    BoundingBox(xmin=x0, ymin=y0, xmax=x1, ymax=y1)

There are no default values for keywords `xmin`, `xmax`, `ymin` and `ymax` so
all must be specified.

The coordinates of a `BoundingBox`, say `box`, can be retrieved as follows:

    box.xmin  or  box[1]  ->  xmin
    box.xmax  or  box[2]  ->  xmax
    box.ymin  or  box[3]  ->  ymin
    box.ymax  or  box[4]  ->  ymax

or:

    xmin, xmax, ymin, ymax = box

See also [`Point`](@ref), [`interior`](@ref), [`exterior`](@ref).

"""
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

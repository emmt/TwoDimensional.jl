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
    AbstractPoint{T}

is the abstract type of objects with at least 2 properties: `x` and `y`, their
respective abscissa and ordinate, both of type `T`.

See also: [`Point`](@ref), [`WeightedPoint`](@ref).

"""
abstract type AbstractPoint{T<:Real} end

"""
    Point(x,y)
    Point((x,y))

yield an instance of a 2D point of coordinates `(x,y)`.

A point may be multiplied or divided by a scalar to scale its coordinates.  The
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

See also: [`WeightedPoint`](@ref), [`AbstractPoint`](@ref).

"""
struct Point{T} <: AbstractPoint{T}
    x::T
    y::T
end

"""
    TwoDimensional.PointLike

is the union of types that may be used to specify a point in `TwoDimensional`
package.

"""
const PointLike = Union{AbstractPoint,NTuple{2,Real},CartesianIndex{2}}

"""
    WeightedPoint{T}(w,x,y)

yields a weighted point which has 3 fields: `w` its weight, `x` its abscissa
and `y` its ordinate, all of type `T`.  By convention `w ≥ 0` but this is not
checked for efficiency reasons.

See also: [`Point`](@ref), [`AbstractPoint`](@ref).

""" WeightedPoint

struct WeightedPoint{T<:AbstractFloat}  <: AbstractPoint{T}
    w::T # weight
    x::T # abscissa
    y::T # ordinate
end

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
struct BoundingBox{T<:Real}
    xmin::T
    xmax::T
    ymin::T
    ymax::T
end

"""
    TwoDimensional.BoundingBoxLike

is the union of types that may be used to specify a bounding-box in
`TwoDimensional` package.

"""
const BoundingBoxLike = Union{BoundingBox,NTuple{2,AbstractPoint},
                              NTuple{4,Real},NTuple{2,CartesianIndex{2}},
                              NTuple{2,AbstractUnitRange{<:Integer}},
                              CartesianIndices{2}}

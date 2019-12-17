#
# TwoDimensional.jl --
#
# Utilities for 2D geometry in Julia.
#
#------------------------------------------------------------------------------
#
# This file if part of the TwoDimensional Julia package licensed under the MIT
# license (https://github.com/emmt/TwoDimensional.jl).
#
# Copyright (C) 2016-2019, Éric Thiébaut.
#
module TwoDimensional

# WeightedPoint is not exported by default because packages may implement a
# different definition.
export
    AbstractPoint,
    AffineTransform,
    BoundingBox,
    Point,
    area,
    center,
    compose,
    distance,
    exterior,
    intercept,
    interior,
    jacobian,
    rotate,
    scale,
    translate

include("basics.jl")

include("xforms.jl")
import .AffineTransforms:
    AffineTransform,
    compose,
    intercept,
    jacobian,
    rotate,
    scale,
    translate

#
# Define some aliases with a 2D suffix to avoid conflicts with other packages.
#
const AbstractPoint2D{T} = AbstractPoint{T}
@doc @doc(AbstractPoint) AbstractPoint2D

const Point2D{T} = Point{T}
@doc @doc(Point) Point2D

const WeightedPoint2D{T} = WeightedPoint{T}
@doc @doc(WeightedPoint) WeightedPoint2D

const BoundingBox2D{T} = BoundingBox{T}
@doc @doc(BoundingBox) BoundingBox2D

const AffineTransform2D{T} = AffineTransform{T}
@doc @doc(AffineTransform) AffineTransform2D


# The TwoDimensional.Suffixed can be imported/used instead of TwoDimensional to
# have 2D suffix for types and avoid conflicts with other packages.
module Suffixed

export
    AbstractPoint2D,
    AffineTransform2D,
    BoundingBox2D,
    Point2D,
    area,
    center,
    compose,
    distance,
    exterior,
    intercept,
    interior,
    jacobian,
    rotate,
    scale,
    translate

import ..TwoDimensional:
    AffineTransform2D,
    BoundingBox2D,
    Point2D,
    area,
    center,
    compose,
    distance,
    exterior,
    intercept,
    interior,
    jacobian,
    rotate,
    scale,
    translate

end # module Suffixed

end # module TwoDimensional

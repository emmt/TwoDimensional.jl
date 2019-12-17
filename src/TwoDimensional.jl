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

"""

`AffineTransform{T}` is an alias for [`TwoDimensional.AbstractPoint{T}`](@ref
AbstractPoint).

""" AbstractPoint2D
const AbstractPoint2D{T} = AbstractPoint{T}

"""

`Point2D{T}` is an alias for [`TwoDimensional.Point{T}`](@ref Point).

""" Point2D
const Point2D{T} = Point{T}

"""

`WeightedPoint2D{T}` is an alias for [`TwoDimensional.WeightedPoint{T}`](@ref
WeightedPoint).

""" WeightedPoint2D
const WeightedPoint2D{T} = WeightedPoint{T}

"""

`BoundingBox2D{T}` is an alias for [`TwoDimensional.BoundingBox{T}`](@ref
BoundingBox).

""" BoundingBox2D
const BoundingBox2D{T} = BoundingBox{T}

"""

`AffineTransform2D{T}` is an alias for
[`TwoDimensional.AffineTransform{T}`](@ref AffineTransform).

""" AffineTransform2D
const AffineTransform2D{T} = AffineTransform{T}


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
    AbstractPoint2D,
    AffineTransform2D,
    BoundingBox2D,
    Point2D,
    WeightedPoint2D,
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

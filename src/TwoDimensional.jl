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
# Copyright (c) 2016-2024, Éric Thiébaut.
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
    # FIXME: compose,
    distance,
    exterior,
    intercept,
    interior,
    jacobian,
    rotate,
    scale,
    translate,
    # Re-export from LinearAlgebra
    det

#using TypeUtils
using Unitless
using Base: @propagate_inbounds
using LinearAlgebra

include("types.jl")
include("basics.jl")
include("xforms.jl")
import .AffineTransforms:
    AffineTransform,
    compose,
    factors_type,
    intercept,
    jacobian,
    offsets_type,
    rotate,
    scale,
    translate
include("suffixed.jl")

end # module TwoDimensional

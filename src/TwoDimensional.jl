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

export
    AbstractPoint,
    AffineTransform,
    BoundingBox,
    Point,
    area,
    center,
    coord_type,
    distance,
    exterior,
    intercept,
    interior,
    jacobian,
    rotate,
    scale,
    translate,
    # Re-exports from LinearAlgebra
    det,
    norm

#using TypeUtils
using Unitless
using Base: @propagate_inbounds
using LinearAlgebra

include("types.jl")
include("basics.jl")
include("boxes.jl")
include("points.jl")
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
include("math.jl")
include("aliases.jl")

end # module TwoDimensional

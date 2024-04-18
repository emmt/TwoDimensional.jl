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

# Imports for extension.
import Base: +, -, *, ∘, /, \, inv
import Base: Float16, Float32, Float64
import Base.MPFR: BigFloat
import LinearAlgebra: ⋅, det

include("types.jl")
include("aliases.jl")
include("basics.jl")
include("boxes.jl")
include("points.jl")
include("xforms.jl")
include("math.jl")

end # module TwoDimensional

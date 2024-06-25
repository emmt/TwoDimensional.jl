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
    Circle,
    Mask,
    Point,
    Polygon,
    Rectangle,
    aperture,
    area,
    center,
    circular_aperture,
    circular_obscuration,
    convert_coord_type,
    coord_type,
    diameter,
    distance,
    exterior,
    forge_mask, forge_mask!,
    interior,
    obscuration,
    polygonal_aperture,
    polygonal_obscuration,
    promote_coord_type,
    radius,
    rectangular_aperture,
    rectangular_obscuration,
    rotate,
    scale,
    translate,
    # Re-exports from LinearAlgebra
    cross,
    det,
    norm

using TypeUtils
using Base: @propagate_inbounds, Callable, Fix1, Fix2, tail
using Base: IteratorSize, SizeUnknown, HasLength, HasShape, IsInfinite
using Base: IteratorEltype, EltypeUnknown, HasEltype
using LinearAlgebra

# Imports for extension.
import Base: ==, +, -, *, ∘, /, \, inv
import LinearAlgebra: ⋅, det

include("types.jl")
include("aliases.jl")
include("common.jl")
include("points.jl")
include("rectangles.jl")
include("circles.jl")
include("polygons.jl")
include("boxes.jl")
include("xforms.jl")
include("winding.jl")
include("math.jl")
include("masks.jl")

end # module TwoDimensional

#
# basics.jl --
#
# Basic methods points and bounding-boxes.
#
#-------------------------------------------------------------------------------
#
# This file is part of the TwoDimensional Julia package licensed under the MIT
# license (https://github.com/emmt/TwoDimensional.jl).
#
# Copyright (c) 2019-2024, Éric Thiébaut.
#

# FIXME: use TypeUtils
as(::Type{T}, x) where {T} = convert(T, x)::T

# Basic constructors (with tuple argument) for points and bounding boxes and
# API to make them indexable iterators.
for (type, len) in ((:Point,         2),
                    (:WeightedPoint, 3),
                    (:BoundingBox,   4),)
    @eval begin
        Base.Tuple(obj::$type) = getfield(obj, 1)
        Base.getindex(obj::$type, i::Integer) = getindex(Tuple(obj), i)
        Base.iterate(obj::$type, i::Int = 1) =
            1 ≤ i ≤ $len ? (obj[i], i + 1) : nothing
        Base.IteratorSize(obj::$type) = Base.IteratorSize(typeof(obj))
        Base.IteratorSize(::Type{<:$type}) = Base.HasLength()
        Base.length(obj::$type) = $len
        Base.show(io::IO, ::MIME"text/plain", obj::$type) = show(io, obj)
    end
    if type === :WeightedPoint
        @eval begin
            $type(vals::NTuple{$len,T}) where {T<:AbstractFloat} = $type{T}(vals)
            $type(vals::NTuple{$len,Real}) = $type(map(float, promote(vals...)))
        end
    else
        @eval begin
            $type(vals::NTuple{$len,T}) where {T} = $type{T}(vals)
            $type(vals::NTuple{$len,Any}) = $type(promote(vals...))
        end
    end
end

# Extend basic methods for abstract points.
Base.eltype(::AbstractPoint{T}) where {T} = T
Base.CartesianIndex(pnt::AbstractPoint{<:Integer}) = CartesianIndex(get_xy(pnt)...)

# Extend basic methods for points and bounding-boxes.
for (type, like, n) in ((:Point,         :PointLike,         2),
                        (:WeightedPoint, :WeightedPointLike, 3),
                        (:BoundingBox,   :BoundingBoxLike,   4),)
    @eval begin
        Base.convert(::Type{T}, obj::T) where {T<:$type} = obj
        Base.convert(::Type{T}, obj::$like) where {T<:$type} = T(obj)
        Base.convert(::Type{Tuple}, obj::$type) = Tuple(obj)
        Base.convert(::Type{NTuple{$n,T}}, obj::$type) where {T} =
            map(#= FIXME: as(T) =# Base.Fix1(as, T), Tuple(obj))
        Base.eltype(obj::$type) = eltype(typeof(obj))
        Base.eltype(::Type{<:$type{T}}) where {T} = T
        Base.promote_type(::Type{$type{T}}, ::Type{$type{T}}) where {T} = $type{T}
        Base.promote_type(::Type{$type{T₁}}, ::Type{$type{T₂}}) where {T₁,T₂} =
            $type{promote_type(T₁,T₂)}
    end
end

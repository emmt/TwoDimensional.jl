"""
    Point{T}(x,y)
    Point{T}((x,y))
    Point{T}(; x=..., y=...)
    Point{T}(; r=..., θ=...)

construct a point given its Cartesian coordinates `(x,y)` (in the 3 first
examples) or its polar coordinates (in the last example) with `r` the distance
to the origin and `θ` the counterclockwise angle relative to `x`-axis.
Parameter `T` is the type used to store coordinates, it may be omitted.

Note that, in `TwoDimensional`, `x` and `y` respectively correspond to the 1st
and 2nd dimensions of 2-dimensional arrays.

A point may be multiplied or divided by a scalar to scale its coordinates. The
addition (resp. subtraction) of two points adds (resp. subtracts) their
coordinates.

The coordinates of a `Point`, say `pnt`, can be retrieved as follows:

    pnt.x  or  pnt[1]  ->  x
    pnt.y  or  pnt[2]  ->  y

or:

    x, y = pnt

the polar coordinates of the point can be retrieved by `hypot(pnt) -> r`,
`abs(pnt) -> r`, or `norm(pnt) -> r`, and by `atan(pnt) -> θ`.

See also [`AbstractPoint`](@ref).

"""
Point{T}(x, y) where {T} = Point{T}((x, y))
Point(x, y) = Point(promote(x, y))
Point(x::T, y::T) where {T} = Point{T}(x, y)

# Keyword-only constructor.
@inline Point(; kwds...) = build(Point; kwds...)
@inline Point{T}(; kwds...) where {T} = build(Point{T}; kwds...)
@inline build(::Type{T}; x=nothing, y=nothing, r=nothing, θ=nothing) where {T<:Point} =
    if (x !== nothing) & (y !== nothing) & (r === nothing) & (θ === nothing)
        return T(x, y)
    elseif (x === nothing) & (y === nothing) & (r !== nothing) & (θ !== nothing)
        sinθ, cosθ = sincos(θ)
        return T(r*cosθ, r*sinθ)
    else
        throw(ArgumentError(
            "exclusively keywords `x` and `y` or `r` and `θ` must be specified"))
    end

Point(pnt::Point) = pnt
Point{T}(pnt::Point{T}) where {T} = pnt
Point{T}(pnt::Point) where {T} = Point{T}(Tuple(pnt))

Point(pnt::PointLike) = Point(get_xy(pnt))
Point{T}(pnt::PointLike) where {T} = Point{T}(get_xy(pnt))

# Properties.
Base.propertynames(::Point) = (:x, :y)
Base.getproperty(pnt::Point, key::Symbol) =
    key === :x ? pnt[1] :
    key === :y ? pnt[2] : throw(KeyError(key))

function Base.show(io::IO, pnt::Point{T}) where {T}
    print(io, "Point{")
    show(io, T)
    print(io, "}(x = "); show(io, pnt.x)
    print(io, ", y = "); show(io, pnt.y)
    print(io, ")")
end

"""
    TwoDimensional.get_x(pnt::TwoDimensional.PointLike) -> x

yields the abscissa of point-like object `pnt`.

See also [`TwoDimensional.get_y`](@ref), [`TwoDimensional.get_xy`](@ref), and
[`TwoDimensional.PointLike`](@ref).

""" get_x

"""
    TwoDimensional.get_y(pnt::TwoDimensional.PointLike) -> y

yields the ordinate of point-like object `pnt`.

See also [`TwoDimensional.get_x`](@ref), [`TwoDimensional.get_xy`](@ref), and
[`TwoDimensional.PointLike`](@ref).

""" get_y

for (c, i) in ((:x, 1), (:y, 2))
    func = Symbol("get_",c)
    @eval begin
        $(func)(pnt::Point) = pnt[$(i)]
        $(func)(pnt::Union{NTuple{2},CartesianIndex{2}}) = pnt[$(i)]
        $(func)(pnt::AbstractPoint) = pnt.$(c)
        $(func)(pnt::PointLike) = get_xy(pnt)[$(i)]
    end
end

"""
    TwoDimensional.get_xy(pnt::TwoDimensional.PointLike) -> (x::T, y::T)

yields a 2-tuple with the abscissa `x` and ordinate `y` of point-like object
`pnt`. This is equivalent to, but more economical than,
`(get_x(pnt),get_y(pnt))`.

See also [`TwoDimensional.get_x`](@ref), [`TwoDimensional.get_y`](@ref) and [`TwoDimensional.PointLike`](@ref).

"""
get_xy(pnt::Point) = Tuple(pnt)
get_xy(pnt::CartesianIndex{2}) = Tuple(pnt)
get_xy(pnt::AbstractPoint) = (pnt.x, pnt.y)
get_xy(pnt::NTuple{2}) = pnt
get_xy(pnt::NTuple{2,Any}) = promote(pnt...)
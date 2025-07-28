"""
    pnt = Point{T}(x,y)
    pnt = Point{T}((x,y))
    pnt = Point{T}(; x=..., y=...)
    pnt = Point{T}(; r=..., θ=...)

construct a point given its Cartesian coordinates `(x,y)` (in the 3 first examples) or its
polar coordinates (in the last example) with `r` the distance to the origin and `θ` the
counterclockwise angle relative to `x`-axis. Parameter `T` is the type used to store
coordinates, it may be omitted.

Note that, in `TwoDimensional`, `x` and `y` respectively correspond to the 1st and 2nd
dimensions of 2-dimensional arrays.

A point may be multiplied or divided by a scalar to scale its coordinates. The addition
(resp. subtraction) of two points adds (resp. subtracts) their coordinates.

Points have the following properties reflecting the keywords accepted by their constructor:

    pnt.x  ->  x::T
    pnt.y  ->  y::T

Points are indexable iterators:

    length(pnt)     -> 2
    firstindex(pnt) -> 1
    lastindex(pnt)  -> 2
    pnt[1]          -> pnt.x
    pnt[2]          -> pnt.y
    first(pnt)      -> pnt.x
    last(pnt)       -> pnt.y
    Tuple(pnt)      -> (pnt.x, pnt.y)

The coordinates of `pnt` can thus be retrieved by:

    x, y = pnt

the polar coordinates of the point can be retrieved by `hypot(pnt) -> r`, `abs(pnt) -> r`,
or `norm(pnt) -> r`, and by `atan(pnt) -> θ`.

See also [`AbstractPoint`](@ref), [`Rectangle`](@ref), and [`BoundingBox`](@ref).

"""
Point(x, y) = Point(promote(x, y))
Point(x::T, y::T) where {T} = Point{T}(x, y)
Point{T}(x, y) where {T} = Point{T}(as(T, x), as(T, y))

# Point constructors for coordinates specified as a 2-tuple.
Point((x, y)::NTuple{2,Any}) = Point(x, y)
Point((x, y)::NTuple{2,T}) where {T} = Point{T}(x, y)
Point{T}((x, y)::NTuple{2,Any}) where {T} = Point{T}(x, y)

# Point constructors for coordinates specified as a Cartesian index.
Point(ind::CartesianIndex{2}) = Point(Tuple(ind))
Point{T}(ind::CartesianIndex{2}) where {T} = Point{T}(Tuple(ind))

# Keyword-only constructor.
@inline Point(; kwds...) = build(Point; kwds...)
@inline Point{T}(; kwds...) where {T} = build(Point{T}; kwds...)
@inline build(::Type{T}; x=nothing, y=nothing, r=nothing, θ=nothing) where {T<:Point} =
    if is_something(x, y) && is_nothing(r, θ)
        return T(x, y)
    elseif is_nothing(x, y) && is_something(r, θ)
        sinθ, cosθ = sincos(θ)
        return T(r*cosθ, r*sinθ)
    else
        throw(ArgumentError(
            "exclusively keywords `x` and `y` or `r` and `θ` must be specified"))
    end

# Convert/copy constructors.
Point(pnt::Point) = pnt
Point(pnt::AbstractPoint) = Point(pnt.x, pnt.y)
Point{T}(pnt::Point{T}) where {T} = pnt
Point{T}(pnt::AbstractPoint) where {T} = Point{T}(pnt.x, pnt.y)
Base.convert(::Type{T}, pnt::T) where {T<:Point} = pnt
Base.convert(::Type{T}, obj::PointLike) where {T<:Point} = T(obj)

# Extend basic methods for abstract points.
Base.eltype(::AbstractPoint{T}) where {T} = T
Base.eltype(::Type{<:AbstractPoint{T}}) where {T} = T
Base.CartesianIndex(pnt::AbstractPoint{<:Integer}) = CartesianIndex(get_xy(pnt)...)

# Properties of points.
Base.propertynames(::Point) = (:x, :y)
Base.getproperty(pnt::Point, key::Symbol) =
    key === :x ? pnt[1] :
    key === :y ? pnt[2] : throw(KeyError(key))

Base.show(io::IO, pnt::Point{T}) where {T} =
    print(io, "Point{", T,
          "}(x = ", pnt.x,
          ", y = ", pnt.y, ")")

"""
    TwoDimensional.point_type(arg)

yields the equivalent point type of `arg`.

If `arg` is a point-like object (or the type of such an object) the result is the point type
of `arg` when converted to a `Point`.

If `arg` is a tuple of point-like objects, the result is the promoted type of the conversion
of each `arg` to a `Point`.

""" point_type
@public point_type
point_type(::Type{<:AbstractPoint{T}}) where {T} = Point{T}
point_type(::Type{CartesianIndex{2}}) = Point{Int}
point_type(::Type{Tuple{X,Y}}) where {X<:Number,Y<:Number} = Point{promote_type(X, Y)}

for type in (:AbstractPoint, :(CartesianIndex{2}), :(NTuple{2,Number}))
    @eval begin
        point_type(pnt::$type) = point_type(typeof(pnt))
        point_type(arr::AbstractVector{<:$type}) = point_type(typeof(arr))
        point_type(::Type{<:AbstractVector{T}}) where {T<:$type} = point_type(T)
        point_type(tup::Tuple{$type,Vararg{$type}}) = _point_type(point_type(first(tup)), tail(tup))
    end
end
_point_type(::Type{T}, ::Tuple{}) where {T<:Point} = T
_point_type(::Type{T}, tup::Tuple) where {T<:Point} =
    _point_type(promote_type(T, point_type(first(tup))), tail(tup))

"""
    TwoDimensional.get_x(pnt::PointLike) -> x

yields the abscissa of point-like object `pnt`.

See also [`TwoDimensional.get_y`](@ref), [`TwoDimensional.get_xy`](@ref), and
[`TwoDimensional.PointLike`](@ref).

""" get_x
@public get_x

"""
    TwoDimensional.get_y(pnt::PointLike) -> y

yields the ordinate of point-like object `pnt`.

See also [`TwoDimensional.get_x`](@ref), [`TwoDimensional.get_xy`](@ref), and
[`TwoDimensional.PointLike`](@ref).

""" get_y
@public get_y

for (c, i) in ((:x, 1), (:y, 2))
    func = Symbol("get_",c)
    @eval begin
        $(func)(pnt::Point) = pnt[$i]
        $(func)(ind::CartesianIndex{2}) = ind[$i]
        $(func)(pnt::AbstractPoint) = pnt.$c
        # NOTE For non-homogeneous 2-tuple, `get_xy` is called to promote
        #      `x` and `y` to the same type.
        $(func)(xy::NTuple{2,T}) where {T<:Number} = xy[$i]
        $(func)(obj::PointLike) = get_xy(obj)[$i]
    end
end

"""
    TwoDimensional.get_xy(pnt::PointLike) -> (x::T, y::T)

yields a 2-tuple with the abscissa `x` and ordinate `y` of point-like object `pnt`. This is
equivalent to, but more economical than, `(get_x(pnt),get_y(pnt))`.

See also [`TwoDimensional.get_x`](@ref), [`TwoDimensional.get_y`](@ref) and
[`TwoDimensional.PointLike`](@ref).

""" get_xy
@public get_xy
get_xy(pnt::Point) = Tuple(pnt)
get_xy(ind::CartesianIndex{2}) = Tuple(ind)
get_xy(pnt::AbstractPoint) = promote(pnt.x, pnt.y) # convert to an homogeneous 2-tuple
get_xy(xy::NTuple{2,T}) where {T<:Number} = xy # homogeneous 2-tuple
get_xy(xy::NTuple{2,Number}) = promote(xy...) # convert to an homogeneous 2-tuple

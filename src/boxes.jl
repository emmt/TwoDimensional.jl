"""
    box = BoundingBox{T}(start, stop)
    box = BoundingBox{T}((xmin,ymin), (xmax,ymax))
    box = BoundingBox{T}(; start=..., stop=...)
    box = BoundingBox{T}(; xmin=..., xmax=..., ymin=..., ymax=...)

construct a rectangular bounding-box with edges aligned with the Cartesian axes
and given the coordinates of 2 opposite corners, `start` and `stop`, whose
coordinates, `(xmin,ymin)` and `(xmax,ymax)`, may be specified as points, as
2-tuples, as 2-dimensional Cartesian indices, or by keywords. Parameter `T` is
the type used to store coordinates, it may be omitted.

Bounding-boxes have the following properties reflecting the keywords
accepted by their constructor:

    box.xmin  -> xmin::T
    box.xmax  -> xmax::T
    box.ymin  -> ymin::T
    box.ymax  -> ymax::T
    box.start -> start::Point{T}
    box.stop  -> stop::Point{T}

A bounding-box is assumed to contain all points of coordinates `(x,y)` such
that `xmin ≤ x ≤ xmax` and `ymin ≤ y ≤ ymax`. If `xmin > xmax` or `ymin >
ymax`, the bounding-box is considered as **empty**. This can be checked with
`isempty(box)`.

Boxes are used to represent grid cells and bounding-boxes of other geometric
shape. Use [`TwoDimensional.Rectangle`](@ref) if you want to define rectangular
masks.

Bounding-boxes are indexable iterators:

    length(box)     -> 2
    firstindex(box) -> 1
    lastindex(box)  -> 2
    box[1]          -> box.start
    box[2]          -> box.stop
    first(box)      -> box.start
    last(box)       -> box.stop
    Tuple(box)      -> (box.start, box.stop)

Hence, the parameters of a bounding-box can be retrieved by:

    start, stop = box
    (xmin, ymin), (xmax, ymax) = box

See also [`Point`](@ref), [`Rectangle`](@ref), [`interior`](@ref), and
[`exterior`](@ref).

"""
BoundingBox(start::Point{T}, stop::Point{T}) where {T} = BoundingBox{T}(start, stop)
BoundingBox(start::Point, stop::Point) = BoundingBox(promote(start, stop)...)
BoundingBox{T}(start::Point, stop::Point) where {T} =
    BoundingBox{T}(Point{T}(start), Point{T}(stop))

# Bounding-box constructors for `start` and `stop` specified as other
# point-like objects than points.
for type in (:AbstractPoint, :(NTuple{2,Any}), :(CartesianIndex{2}))
    @eval begin
        # `start` and `stop` provided as 2 arguments.
        BoundingBox(start::$type, stop::$type) = BoundingBox(Point(start), Point(stop))
        BoundingBox{T}(start::$type, stop::$type) where {T} =
            BoundingBox{T}(Point(start), Point(stop))

        # `start` and `stop` provided as a 2-tuple.
        BoundingBox((start, stop)::NTuple{2,$type}) = BoundingBox(start, stop)
        BoundingBox{T}((start, stop)::NTuple{2,$type}) where {T} =
            BoundingBox{T}(start, stop)
    end
end

# Bounding-box constructors for bounds specified as integer-valued unit-ranges.
function BoundingBox(xrng::AbstractUnitRange{<:Integer},
                     yrng::AbstractUnitRange{<:Integer})
    T = promote_type(eltype(xrng), eltype(yrng))
    return BoundingBox{T}(xrng, yrng)
end
function BoundingBox{T}(xrng::AbstractUnitRange{<:Integer},
                        yrng::AbstractUnitRange{<:Integer}) where {T}
    return BoundingBox{T}(xmin = as(T, first(xrng)), xmax = as(T, last(xrng)),
                          ymin = as(T, first(yrng)), ymax = as(T, last(yrng)))
end
BoundingBox((xrng, yrng)::NTuple{2,AbstractUnitRange{<:Integer}}) = BoundingBox(xrng, yrng)
BoundingBox{T}((xrng, yrng)::NTuple{2,AbstractUnitRange{<:Integer}}) where {T} =
    BoundingBox{T}(xrng, yrng)

# Bounding-box constructors for box limits specified by Cartesian indices.
BoundingBox(inds::CartesianIndices{2}) = BoundingBox(first(inds), last(inds))
BoundingBox{T}(inds::CartesianIndices{2}) where {T} = BoundingBox{T}(first(inds), last(inds))

# Keyword-only bounding-box constructors.
@inline BoundingBox(; kwds...) = build(BoundingBox; kwds...)
@inline BoundingBox{T}(; kwds...) where {T} = build(BoundingBox{T}; kwds...)
@inline function build(::Type{T};
                       xmin=nothing, ymin=nothing, xmax=nothing, ymax=nothing,
                       start=nothing, stop=nothing) where{T<:BoundingBox}
    if is_something(start, stop) && is_nothing(xmin, ymin, xmax, ymax)
        return T(start, stop)
    elseif is_nothing(start, stop) && is_something(xmin, ymin, xmax, ymax)
        return T((xmin, ymin), (xmax, ymax))
    else
        throw(ArgumentError(
            "exclusively keywords `start` and `stop` or keywords `xmin`, `ymin`, `xmax`, and `ymax` must be specified"))
    end
end

"""
    BoundingBox{T}()

yields an empty bounding-box for coordinate type `T`.

"""
BoundingBox{T}() where {T} = BoundingBox((zero(T), zero(T)), (-oneunit(T), -oneunit(T)))

# Convert/copy constructors.
BoundingBox(box::BoundingBox) = box
BoundingBox{T}(box::BoundingBox{T}) where {T} = box
BoundingBox{T}(box::BoundingBox) where {T} = BoundingBox{T}(Tuple(box)...)
Base.convert(::Type{T}, box::T) where {T<:BoundingBox} = box
Base.convert(::Type{T}, obj::BoundingBoxLike) where {T<:BoundingBox} = T(obj)

# Make bounding-boxes indexable iterators.
Base.length(        box::BoundingBox) = 2
Base.firstindex(    box::BoundingBox) = 1
Base.lastindex(     box::BoundingBox) = 2
Base.first(         box::BoundingBox) = box[1]
Base.last(          box::BoundingBox) = box[2]
Base.eltype(        box::BoundingBox) = eltype(typeof(box))
Base.IteratorSize(  box::BoundingBox) = IteratorSize(typeof(box))
Base.IteratorEltype(box::BoundingBox) = IteratorEltype(typeof(box))
Base.eltype(        ::Type{<:BoundingBox{T}}) where {T} = Point{T}
Base.IteratorSize(  ::Type{<:BoundingBox}) = HasLength()
Base.IteratorEltype(::Type{<:BoundingBox}) = HasEltype()
@inline Base.iterate(box::BoundingBox, i::Int = 1) =
    i == 1 ? (first(box), 2) :
    i == 2 ? (last( box), 3) : nothing

# Properties of bounding-boxes.
Base.propertynames(::BoundingBox) = (:start, :stop, :xmin, :xmax, :ymin, :ymax,)
Base.getproperty(box::BoundingBox, key::Symbol) =
    key === :start ? first(box)   :
    key === :stop  ? last( box)   :
    key === :xmin  ? first(box).x :
    key === :ymin  ? first(box).y :
    key === :xmax  ? last( box).x :
    key === :ymax  ? last( box).y : throw(KeyError(key))

Base.show(io::IO, box::BoundingBox{T}) where {T} =
    print(io, "BoundingBox{", T,
          "}(xmin = ", box.xmin,
          ", xmax = ", box.xmax,
          ", ymin = ", box.ymin,
          ", ymax = ", box.ymax, ")")

"""
    BoundingBox{T}(obj)

yields the bounding-box of the geometric object `obj`. If the coordinate type
`T` is not provided, [`T = coord_type(obj)`](@ref TwoDimensional.coord_type) is
assumed.

This can be used to compute the union or the intersection of the bounding-boxes
of objects:

    mapreduce(BoundingBox, ∪, (obj1, obj2, obj3, ...))
    mapreduce(BoundingBox, ∩, [obj1, obj2, obj3, ...])

"""
BoundingBox(obj::GeometricObject{T}) where {T} = BoundingBox{T}(obj)
BoundingBox{T}(obj::MaskElement) where {T} = BoundingBox{T}(shape(obj))
BoundingBox{T}(pnt::Point) where {T} = BoundingBox{T}(pnt, pnt)
BoundingBox{T}(rect::Rectangle) where {T} = BoundingBox{T}(first(rect), last(rect))
function BoundingBox{T}(circ::Circle) where {T}
    c = center(circ)
    r = radius(circ)
    s = Point(r, r)
    return BoundingBox{T}(c - s, c + s)
end
function BoundingBox{T}(poly::Polygon) where {T}
    xmin = ymin = typemax(coord_type(poly))
    xmax = ymax = typemin(coord_type(poly))
    @inbounds for i in eachindex(poly)
        x, y = poly[i]
        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
    end
    return BoundingBox{T}((xmin, ymin), (xmax, ymax))
end

"""
    BoundingBox(f, A::AbstractMatrix)

yields the minimal bounding-box of the entries of `A` such that `f(A[i,j])` is
true. This can be used to extract this rectangular region:

    A[BoundingBox(f, A)]

If `A` is a matrix of Booleans, `f` is assumed to be the identity if not
specified.

"""
BoundingBox(A::AbstractMatrix{Bool}) = BoundingBox(identity, A)

# See
# https://stackoverflow.com/questions/9852159/calculate-bounding-box-of-arbitrary-pixel-based-drawing
# for the basic ideas under the following algorithm.
function BoundingBox(f::Function, A::AbstractMatrix)
    I, J = axes(A)
    i0, i1 = get_axis_bounds(I)
    j0, j1 = get_axis_bounds(J)
    imin = jmin = typemax(Int)
    imax = jmax = typemin(Int)
    # Assuming column-major order, first start by scanning rows to narrow the
    # subsequent searches along columns.
    #
    # 1. Find bottom bound `jmin` by scanning rows from bottom to top.
    flag = false
    @inbounds for j in j0:j1, i in i0:i1
        if f(A[i,j])
            # This definitively set the value of `jmin` and gives limits for
            # the other bounds.
            imin = imax = i
            jmin = jmax = j
            flag = true
            break
        end
    end
    if flag
        # 2. Find top bound `jmax` by scanning rows from top to bottom.  No
        #    needs to go beyond `jmax+1`.
        @inbounds for j in j1:-1:jmax+1, i in i0:i1
            if f(A[i,j])
                jmax = j
                imin = min(imin, i)
                imax = max(imax, i)
                break
            end
        end
        # 3. Find leftmost bound `imin` by scanning columns from left to right.
        #    No needs to go beyond `imin-1`.
        @inbounds for i in i0:imin-1, j in jmin:jmax
            if f(A[i,j])
                imin = i
                imax = max(imax, i)
                break
            end
        end
        # 4. Find rightmost bound `imax` by scanning columns from right to
        #    left.  No needs to go beyond `imax+1`.
        @inbounds for i in i1:-1:imax+1, j in jmin:jmax
            if f(A[i,j])
                imax = i
                break
            end
        end
    end
    return BoundingBox(imin:imax, jmin:jmax)
end

"""
    TwoDimensional.get_axis_bounds(I) = (i0,i1)

yields the bounds `i0` and `i1` of index range `I` as a 2-tuple of `Int`'s and
such that `i0:i1` represents the same indices as `I` (although not in the same
order if `step(I) < 0`). If `step(I)` is not equal to ±1, an `ArgumentError`
exception is thrown.

"""
@inline function get_axis_bounds(I::AbstractRange{<:Integer})
    i0, i1, s = Int(first(I)), Int(last(I)), step(I)
    s == +one(s) && return (i0, i1)
    s == -one(s) && return (i1, i0)
    throw(ArgumentError("expecting a range with a step equal to ±1, got $s"))
end

# Empty bounding and unlimited boxes.
BoundingBox{T}(::Nothing) where {T<:Real} = typemin(BoundingBox{T})
Base.typemin(::Type{BoundingBox{T}}) where {T<:Real} =
    BoundingBox(xmin = typemax(T), xmax = typemin(T),
                ymin = typemax(T), ymax = typemin(T))
Base.typemax(::Type{BoundingBox{T}}) where {T<:Real} =
    BoundingBox(xmin = typemin(T), xmax = typemax(T),
                ymin = typemin(T), ymax = typemax(T))

# Conversion of bounding-boxes to Cartesian indices.
Base.CartesianIndices(box::BoundingBox) = CartesianIndices(axes(box))

Base.isempty(box::BoundingBox) = (box.xmin > box.xmax)|(box.ymin > box.ymax)

Base.size(box::BoundingBox) = map(length, axes(box))
Base.size(box::BoundingBox, d::Integer) = length(axes(box, d))

Base.axes(box::BoundingBox{<:Integer}) = (UnitRange{Int}(box.xmin, box.xmax),
                                          UnitRange{Int}(box.ymin, box.ymax))
Base.axes(box::BoundingBox{<:Integer}, d::Integer) =
    d == 1 ? UnitRange{Int}(box.xmin, box.xmax) :
    d == 2 ? UnitRange{Int}(box.ymin, box.ymax) :
    d > 2  ? UnitRange{Int}(1, 1) : error("invalid dimension index")

Base.axes(box::BoundingBox) = throw_axes_restricted_box_with_integer_coordinates()
Base.axes(box::BoundingBox, d::Integer) = throw_axes_restricted_box_with_integer_coordinates()
@noinline throw_axes_restricted_box_with_integer_coordinates() = throw(ArgumentError(
    "`axes(box)` is restricted to bounding-boxes with integer coordinates"))

# Use bounding-boxes to extract a sub-array or a view.
Base.view(arr::AbstractMatrix, box::BoundingBox) = view(arr, axes(box)...)
@inline @propagate_inbounds Base.getindex(arr::AbstractMatrix, box::BoundingBox) =
    arr[axes(box)...]

"""
    TwoDimensional.grow(box, dx, dy=dx)

yields a new bounding-box object corresponding to the input `box` object with 1st and
2nd dimensions respectively grown by `dx` and `dy`.

Note that the algebraic (not absolute) values of `dx` and `dy` are applied.
Hence, if `dx` and `dy` are both negative, the bounding-box is effectively
shrunk by `abs(dx)` and `abs(dy)`.

See also [`TwoDimensional.shrink`](@ref).

"""
grow(box::BoundingBox{T}, delta::Point{T}) where {T} = BoundingBox(first(box) - delta,
                                                                   last( box) + delta)

"""
    TwoDimensional.shrink(box, dx, dy=dx)

yields a new bounding-box object corresponding to the input `box` object with 1st and
2nd dimensions respectively shrunk by `dx` and `dy`.

Note that the algebraic (not absolute) values are applied. Hence, if `dx` and
`dy` are both negative, the bounding-box is effectively grown by `abs(dx)` and
`abs(dy)`.

See also [`TwoDimensional.grow`](@ref).

"""
shrink(box::BoundingBox{T}, delta::Point{T}) where {T} = BoundingBox(first(box) + delta,
                                                                     last( box) - delta)

# Encode methods having identical code for `grow` and `shrink`.
for func in (:grow, :shrink)
    @eval begin
        $func(box::BoundingBox, dx::Number, dy::Number=dx) = $func(box, Point(dx, dy))
        $func(box::BoundingBox, delta::Point) = $func(promote_coord_type(box, delta)...)
    end
end

"""
    interior([T,] box)

yields the largest bounding-box with integer valued bounds and which is
contained by the bounding-box `box`. Optional argument `T` is to specify the
type of the result or the type of the coordinates of the result which is the
same as `box` by default.

See also [`exterior`](@ref), [`ceil`](@ref ceil(::Point)) and [`floor`](@ref
floor(::Point)).

"""
interior(::Type{T}, box::BoundingBox) where {T} = BoundingBox(ceil(T, first(box)),
                                                              floor(T, last(box)))

"""
    exterior([T,] box)

yields the smallest bounding-box with integer valued bounds and which contains
the bounding-box `box`. Optional argument `T` is to specify the type of the
result or the type of the coordinates of the result which is the same as `box`
by default.

See also [`interior`](@ref), [`ceil`](@ref ceil(::Point)) and [`floor`](@ref
floor(::Point)).

"""
exterior(::Type{T}, box::BoundingBox) where {T} = BoundingBox(floor(T, first(box)),
                                                              ceil(T, last(box)))

# Encode methods having identical code for `interior` and `exterior`.
for func in (:interior, :exterior)
    @eval begin
        # Deal with `T` and box coordinates all integers.
        $func(::Type{T}, box::BoundingBox{T}) where {T<:Integer} = box
        $func(::Type{T}, box::BoundingBox{<:Integer}) where {T<:Integer} =
            BoundingBox{T}(box)

        # Provide coordinate type for `interior` and `exterior`.
        $func(box::BoundingBox{T}) where {T} = $func(T, box)
        $func(::Type{BoundingBox}, box::BoundingBox{T}) where {T} = $func(T, box)
        $func(::Type{BoundingBox{T}}, box::BoundingBox) where {T} = $func(T, box)
        $func(::Type{BoundingBox{T}}, box::BoundingBox{T}) where {T} = $func(T, box)
    end
end

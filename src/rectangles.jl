"""
    rect = Rectangle{T}(start, stop)
    rect = Rectangle{T}((x0,y0), (x1,y1))
    rect = Rectangle{T}(; start=..., stop=...)
    rect = Rectangle{T}(; x0=..., x1=..., y0=..., y1=...)

construct a rectangular rectangle with edges aligned with the Cartesian axes
and given the coordinates of 2 opposite corners, `start` and `stop`, whose
coordinates, `(x0,y0)` and `(x1,y1)`, may be specified as points, as 2-tuples,
as 2-dimensional Cartesian indices, or by keywords. Parameter `T` is the type
used to store coordinates, it may be omitted.

Rectangles have the following properties reflecting the keywords accepted by
their constructor:

    rect.x0  -> min(x0, x1)::T
    rect.x1  -> max(x0, x1)::T
    rect.y0  -> min(y0, y1)::T
    rect.y1  -> max(y0, y1)::T
    rect.start -> start::Point{T}
    rect.stop  -> stop::Point{T}

Note that the coordinates are sorted. A rectangle is never empty and contains
at least a single point.

Rectangles are indexable iterators:

    length(rect)     -> 2
    firstindex(rect) -> 1
    lastindex(rect)  -> 2
    rect[1]          -> rect.start
    rect[2]          -> rect.stop
    first(rect)      -> rect.start
    last(rect)       -> rect.stop
    Tuple(rect)      -> (rect.start, rect.stop)

Hence, the parameters of a rectangle can be retrieved by:

    start, stop = rect
    (x0, y0), (x1, y1) = rect

See also [`Point`](@ref), [`BoundingBox`](@ref), [`interior`](@ref), and
[`exterior`](@ref).

"""
Rectangle(start::Point{T}, stop::Point{T}) where {T} = Rectangle{T}(start, stop)
Rectangle(start::Point, stop::Point) = Rectangle(promote(start, stop)...)
Rectangle{T}(start::Point, stop::Point) where {T} =
    Rectangle{T}(Point{T}(start), Point{T}(stop))

# Rectangle constructors for `start` and `stop` specified as other
# point-like objects than points.
for type in (:AbstractPoint, :(NTuple{2,Any}), :(CartesianIndex{2}))
    @eval begin
        # `start` and `stop` provided as 2 arguments.
        Rectangle(start::$type, stop::$type) = Rectangle(Point(start), Point(stop))
        Rectangle{T}(start::$type, stop::$type) where {T} =
            Rectangle{T}(Point(start), Point(stop))

        # `start` and `stop` provided as a 2-tuple.
        Rectangle((start, stop)::NTuple{2,$type}) = Rectangle(start, stop)
        Rectangle{T}((start, stop)::NTuple{2,$type}) where {T} =
            Rectangle{T}(start, stop)
    end
end

# Keyword-only rectangle constructors.
@inline Rectangle(; kwds...) = build(Rectangle; kwds...)
@inline Rectangle{T}(; kwds...) where {T} = build(Rectangle{T}; kwds...)
@inline function build(::Type{T};
                       x0=nothing, y0=nothing, x1=nothing, y1=nothing,
                       start=nothing, stop=nothing) where{T<:Rectangle}
    if is_something(start, stop) && is_nothing(x0, y0, x1, y1)
        return T(start, stop)
    elseif is_nothing(start, stop) && is_something(x0, y0, x1, y1)
        return T((x0, y0), (x1, y1))
    else
        throw(ArgumentError(
            "exclusively keywords `start` and `stop` or keywords `x0`, `y0`, `x1`, and `y1` must be specified"))
    end
end

# Build rectangles from other geometric objects.
Rectangle(pnt::Point) = Rectangle(pnt, pnt)
Rectangle{T}(pnt::Point) where {T} = Rectangle{T}(pnt, pnt)
Rectangle(pnt::AbstractPoint) = Rectangle(Point(pnt))
Rectangle{T}(pnt::AbstractPoint) where {T} = Rectangle{T}(Point(pnt))
Rectangle(obj::RectangularObject) = Rectangle(first(obj), last(obj))
Rectangle{T}(obj::RectangularObject) where {T} = Rectangle{T}(first(obj), last(obj))

# Convert/copy constructors.
Rectangle(rect::Rectangle) = rect
Rectangle{T}(rect::Rectangle{T}) where {T} = rect
Rectangle{T}(rect::Rectangle) where {T} = Rectangle{T}(Tuple(rect)...)
Base.convert(::Type{T}, rect::T) where {T<:Rectangle} = rect
Base.convert(::Type{T}, obj::RectangleLike) where {T<:Rectangle} = T(obj)

# Make rectangles indexable iterators.
Base.length(        rect::Rectangle) = 2
Base.firstindex(    rect::Rectangle) = 1
Base.lastindex(     rect::Rectangle) = 2
Base.first(         rect::Rectangle) = rect[1]
Base.last(          rect::Rectangle) = rect[2]
Base.eltype(        rect::Rectangle) = eltype(typeof(rect))
Base.IteratorSize(  rect::Rectangle) = IteratorSize(typeof(rect))
Base.IteratorEltype(rect::Rectangle) = IteratorEltype(typeof(rect))
Base.eltype(        ::Type{<:Rectangle{T}}) where {T} = Point{T}
Base.IteratorSize(  ::Type{<:Rectangle}) = HasLength()
Base.IteratorEltype(::Type{<:Rectangle}) = HasEltype()
@inline Base.iterate(rect::Rectangle, i::Int = 1) =
    i == 1 ? (first(rect), 2) :
    i == 2 ? (last( rect), 3) : nothing

# Properties of rectangles.
Base.propertynames(::Rectangle) = (:start, :stop, :x0, :x1, :y0, :y1,)
Base.getproperty(rect::Rectangle, key::Symbol) =
    key === :start ? first(rect)   :
    key === :stop  ? last( rect)   :
    key === :x0  ? first(rect).x :
    key === :y0  ? first(rect).y :
    key === :x1  ? last( rect).x :
    key === :y1  ? last( rect).y : throw(KeyError(key))

Base.show(io::IO, rect::Rectangle{T}) where {T} =
    print(io, "Rectangle{", T,
          "}(x0 = ", rect.x0,
          ", y0 = ", rect.y0,
          ", x1 = ", rect.x1,
          ", y1 = ", rect.y1, ")")

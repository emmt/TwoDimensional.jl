# Extend `map` and `Broadcast.broadcasted` for points and bounding-boxes.
# NOTE: `Base.Callable = Union{Type,Function}`
for type in (:Point, :BoundingBox)
    @eval begin
        Broadcast.broadcasted(::Type{T}, obj::$type) where {T} = map(T, obj)
        Broadcast.broadcasted(f::Function, obj::$type) = map(f, obj)
    end
end
@inline Base.map(f, pnt::Point) = Point(map(f, Tuple(pnt)))
@inline Base.map(f, box::BoundingBox; swap::Bool = false) =
    BoundingBox(map(f, swap ? (box.xmax, box.xmin, box.ymax, box.ymin) : Tuple(box)))

# Check for equality.
Base.:(==)(a::Point, b::Point) =
    (a.x == b.x) &
    (a.y == b.y)
Base.:(==)(a::BoundingBox, b::BoundingBox) =
    (a.xmin == b.xmin) &
    (a.ymin == b.ymin) &
    (a.xmax == b.xmax) &
    (a.ymax == b.ymax)

# Check for approximate equality.
Base.isapprox(a::Point, b::Point; kwds...) =
    isapprox(a.x, b.x; kwds...) &&
    isapprox(a.y, b.y; kwds...)
Base.isapprox(a::BoundingBox, b::BoundingBox; kwds...) =
    isapprox(a.xmin, b.xmin; kwds...) &&
    isapprox(a.ymin, b.ymin; kwds...) &&
    isapprox(a.xmax, b.xmax; kwds...) &&
    isapprox(a.ymax, b.ymax; kwds...)

# Extend ∈ operator.
Base.in(pnt::AbstractPoint, box::BoundingBox) =
    (box.xmin ≤ pnt.x ≤ box.xmax) & (box.ymin ≤ pnt.y ≤ box.ymax)

# Extend ⊆ operator.
Base.issubset(A::BoundingBox, B::BoundingBox) =
    (isempty(A)|((A.xmin ≥ B.xmin)&(A.xmax ≤ B.xmax)&
                 (A.ymin ≥ B.ymin)&(A.ymax ≤ B.ymax)))

# Union of bounding-boxes:
Base.:(∪)(A::BoundingBox, B::BoundingBox) =
    BoundingBox(min(A.xmin, B.xmin), max(A.xmax, B.xmax),
                min(A.ymin, B.ymin), max(A.ymax, B.ymax))

# Intersection of bounding-boxes:
Base.:(∩)(A::BoundingBox, B::BoundingBox) =
    BoundingBox(max(A.xmin, B.xmin), min(A.xmax, B.xmax),
                max(A.ymin, B.ymin), min(A.ymax, B.ymax))

# Unary plus does nothing.
Base.:(+)(pnt::Point) = pnt
Base.:(+)(box::BoundingBox) = box

# Unary minus negates coordinates.
Base.:(-)(pnt::Point) = map(-, pnt)
Base.:(-)(box::BoundingBox) = map(-, box; swap=true)

# Extend "trait" methods for points and bounding-boxes.
Base.zero(obj::Union{Point,BoundingBox}) = zero(typeof(obj))
Base.one(obj::Union{Point,BoundingBox}) = one(typeof(obj))

# Scaling of points and corresponding multiplicative identity.
Base.one(::Type{Point{T}}) where {T} = one(T)
Base.:(*)(pnt::Point, α::Number) = α*pnt
Base.:(*)(α::Number, pnt::Point) = map(Base.Fix1(*,α), pnt)
Base.:(/)(pnt::Point, α::Number) = map(Base.Fix2(/,α), pnt)
Base.:(\)(α::Number, pnt::Point) = pnt/α

# Scaling of bounding-box bounds (e.g. to change units) and corresponding
# multiplicative identity.
Base.one(::Type{BoundingBox{T}}) where {T} = one(T)
Base.:(*)(box::BoundingBox, α::Number) = α*box
Base.:(*)(α::Number, box::BoundingBox) = map(Base.Fix1(*,α), box; swap = α < zero(α))
Base.:(\)(α::Number, box::BoundingBox) = box/α
Base.:(/)(box::BoundingBox, α::Number) = map(Base.Fix2(/,α), box; swap = α < zero(α))

# Addition and subtraction of points and corresponding addtive identity.
Base.zero(::Type{Point{T}}) where {T} = Point(zero(T),zero(T))
Base.:(+)(A::Point, B::Point) = Point(A.x + B.x, A.y + B.y)
Base.:(-)(A::Point, B::Point) = Point(A.x - B.x, A.y - B.y)

# Addition amd subtraction of bounding-boxes (following the rules of the
# addition and subtraction of sets) and corresponding addtive identity.
Base.zero(::Type{BoundingBox{T}}) where {T} = BoundingBox(zero(T),zero(T),zero(T),zero(T))
Base.:(+)(A::BoundingBox, B::BoundingBox) =
    isempty(A) || isempty(B) ? BoundingBox{promote_type(eltype(A), eltype(B))}(nothing) :
    BoundingBox(map(+, Tuple(A), Tuple(B)))
Base.:(-)(A::BoundingBox, B::BoundingBox) = A + (-B)

# Add or remove a margin δ to a bounding-box box.
Base.:(+)(box::BoundingBox, δ::Number) =
    BoundingBox(box.xmin - δ, box.xmax + δ,
                box.ymin - δ, box.ymax + δ)
Base.:(-)(box::BoundingBox, δ::Number) =
    BoundingBox(box.xmin + δ, box.xmax - δ,
                box.ymin + δ, box.ymax - δ)

# Translate a bounding-box.
Base.:(+)(box::BoundingBox, pnt::Point) =
    BoundingBox(box.xmin + pnt.x, box.xmax + pnt.x,
                box.ymin + pnt.y, box.ymax + pnt.y)
Base.:(-)(box::BoundingBox, pnt::Point) =
    BoundingBox(box.xmin - pnt.x, box.xmax - pnt.x,
                box.ymin - pnt.y, box.ymax - pnt.y)


"""
    round([T,] obj::Union{Point,BoundingBox}, [r::RoundingMode])

yields the object that is the nearest to `obj` by rounding its coordinates to
nearest integral values. Argument `T` can be the type of the result (a point or
a bounding-box) or the type of the coordinates of the result. Rounding mode may
be specified by optional argument `r`, the default being the same as the
`round` method for a scalar value.

For points, see also: [`floor(::Point)`](@ref), [`ceil(::Point)`](@ref).

For bounding-boxes, see also: [`interior`](@ref), [`exterior`](@ref).

""" Base.round

# NOTE: Using the following functor is a bit faster than map with an anonymous
#       function.
struct Round{T,R} <: Function
    r::R
    Round{T}(r::R) where {T,R} = new{T,R}(r)
end
Round{T}() where {T} = Round{T}(nothing)
(f::Round{Nothing,Nothing})(x) = round(x)
(f::Round{Nothing,<:RoundingMode})(x) = round(x, f.r)
(f::Round{T,Nothing})(x) where {T<:Number} = round(T, x)
(f::Round{T,<:RoundingMode})(x) where {T<:Number} = round(T, x, f.r)

for type in (:Point, :BoundingBox)
    @eval begin
        Base.round(obj::$type) = map(round, obj)
        Base.round(obj::$type, r::RoundingMode) = map(Round{Nothing}(r), obj)
        Base.round(::Type{T}, obj::$type) where {T<:Number} = map(Round{T}(), obj)
        Base.round(::Type{T}, obj::$type, r::RoundingMode) where {T<:Number} = map(Round{T}(r), obj)

        Base.round(::Type{$type}, obj::$type) = round(obj)
        Base.round(::Type{$type}, obj::$type, r::RoundingMode) = round(obj, r)
        Base.round(::Type{$type{T}}, obj::$type) where {T} = round(T, obj)
        Base.round(::Type{$type{T}}, obj::$type, r::RoundingMode) where {T} = round(T, obj, r)
    end
end

"""
    floor([T,] pnt::Point)
    floor([T,] box::BoundingBox)

yield the point with the largest integer coordinates smaller or equal those of
`pnt` or the bounding-box resulting from applying this function to the content
of `box`. Argument `T` can be the type of the result or the type of the
coordinates of the result.

See also: [`round`](@ref), [`ceil`](@ref).

""" Base.floor

"""
    ceil([T,] pnt::Point)
    ceil([T,] box::BoundingBox)

yields the point with the smallest integer coordinates larger or equal those of
`pnt` or the bounding-box resulting from applying this function to the content
of `box`. Argument `T` can be the type of the result or the type of the
coordinates of the result.

See also: [`round`](@ref), [`floor`](@ref).

""" Base.ceil

for type in (:Point, :BoundingBox), func in (:floor, :ceil)
    @eval begin
        Base.$func(obj::$type) = map($func, obj)
        Base.$func(::Type{T}, obj::$type) where {T<:Number} = map(x -> $func(T, x), obj)

        Base.$func(::Type{$type}, obj::$type) = $func(obj) # FIXME:
        Base.$func(::Type{$type{T}}, obj::$type) where {T} = $func(T, obj)
    end
end

function Base.clamp(pnt::Point, box::BoundingBox)
    T = promote_type(eltype(pnt), eltype(box))
    x = clamp(as(T, pnt.x), as(T, box.xmin), as(T, box.xmax))
    y = clamp(as(T, pnt.y), as(T, box.ymin), as(T, box.ymax))
    return Point{T}(x, y)
end

# Mathematic functions for points.
Base.hypot(pnt::Point) = hypot(pnt.x, pnt.y)
Base.abs(pnt::Point) = hypot(pnt)
Base.abs2(pnt::Point) = abs2(pnt.x) + abs2(pnt.y)
LinearAlgebra.norm(pnt::Point) = abs(pnt)
Base.Math.atan(pnt::Point) = atan(pnt.y, pnt.x)
inner(a::Point, b::Point) = a.x*b.x + a.y*b.y # scalar/inner product
outer(a::Point, b::Point) = a.x*b.y - a.y*b.x # outer product
Base.:(*)(a::Point, b::Point) = inner(a, b) # scalar product

"""
    distance(A, B)

yields the Euclidean distance between the 2 points `A` and `B`.

"""
distance(A::Point, B::Point) = hypot(A - B)

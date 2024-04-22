# Check for equality.
Base.:(==)(a::Point, b::Point) =
    (a.x == b.x) & (a.y == b.y)
Base.:(==)(a::Rectangle, b::Rectangle) =
    (first(a) == first(b)) & (last(a) == last(b))
Base.:(==)(a::BoundingBox, b::BoundingBox) =
    isempty(a) ? isempty(b) : ((first(a) == first(b)) & (last(a) == last(b)))

# Check for approximate equality.
Base.isapprox(a::Point, b::Point; kwds...) =
    isapprox(a.x, b.x; kwds...) && isapprox(a.y, b.y; kwds...)
Base.isapprox(a::Rectangle, b::Rectangle; kwds...) =
    isapprox(first(a), first(b)) && isapprox(last(a), last(b))
Base.isapprox(a::BoundingBox, b::BoundingBox; kwds...) =
    isapprox(first(a), first(b)) && isapprox(last(a), last(b))

# Extend ∈ operator.
Base.in(pnt::AbstractPoint, box::BoundingBox) =
    (box.xmin ≤ pnt.x ≤ box.xmax) & (box.ymin ≤ pnt.y ≤ box.ymax)
Base.in(pnt::AbstractPoint, rect::Rectangle) =
    (rect.x0 ≤ pnt.x ≤ rect.x1) & (rect.y0 ≤ pnt.y ≤ rect.y1)
Base.in(A::Rectangle, B::BoundingBox) =
    (first(A) ∈ B) && (last(A) ∈ B)

# Extend ⊆ operator.
Base.issubset(A::BoundingBox, B::BoundingBox) =
    (isempty(A)|((A.xmin ≥ B.xmin)&(A.xmax ≤ B.xmax)&
                 (A.ymin ≥ B.ymin)&(A.ymax ≤ B.ymax)))

# Operator ∪, union of bounding-boxes.
Base.union(A::BoundingBox, B::BoundingBox) =
    BoundingBox(xmin = min(A.xmin, B.xmin), xmax = max(A.xmax, B.xmax),
                ymin = min(A.ymin, B.ymin), ymax = max(A.ymax, B.ymax))

# Operator ∩, intersection of bounding-boxes.
Base.intersect(A::BoundingBox, B::BoundingBox) =
    BoundingBox(xmin = max(A.xmin, B.xmin), xmax = min(A.xmax, B.xmax),
                ymin = max(A.ymin, B.ymin), ymax = min(A.ymax, B.ymax))

# Just permute operands for some common operations on geometric objects.
*(a::GeometricObject, b::Number) = b*a
\(a::Number, b::GeometricObject) = b/a
+(a::Point{T}, b::GeometricObject{T}) where {T} = b + a

# Promote operands to a common coordinate type for some geometric operations.
+(a::GeometricObject, b::Point) = +(promote_coord_type(a, b)...)
+(a::Point, b::GeometricObject) = +(promote_coord_type(a, b)...)
+(a::Point, b::Point) = +(promote_coord_type(a, b)...)
-(a::GeometricObject, b::Point) = -(promote_coord_type(a, b)...)
-(a::Point, b::GeometricObject) = -(promote_coord_type(a, b)...)
-(a::Point, b::Point) = -(promote_coord_type(a, b)...)

# Unary plus does nothing.
+(a::GeometricObject) = a

# Unary minus negates coordinates and should as multiplying by -1.
-(pnt::Point) = apply(-, pnt)
-(rect::Rectangle) = apply(-, rect)
-(box::BoundingBox) = apply(-, box; swap = true)

# Scaling of geometric objects and corresponding multiplicative identity.
Base.one(obj::GeometricObject) = one(typeof(obj))

Base.one(::Type{Point{T}}) where {T} = one(T)
*(α::Number, pnt::Point) = apply(Fix1(*, α), pnt)
/(pnt::Point, α::Number) = apply(Fix2(/, α), pnt)

Base.one(::Type{Rectangle{T}}) where {T} = one(T)
*(α::Number, rect::Rectangle) =  apply(Fix1(*, α), rect)
/(rect::Rectangle, α::Number) =  apply(Fix2(/, α), rect)

Base.one(::Type{BoundingBox{T}}) where {T} = one(T)
*(α::Number, box::BoundingBox) = apply(Fix1(*, α), box; swap = α < zero(α))
/(box::BoundingBox, α::Number) = apply(Fix2(/, α), box; swap = α < zero(α))

# Translate a geometric object by adding or subtracting a point and
# corresponding addtive identity.
Base.zero(obj::GeometricObject) = zero(typeof(obj))

Base.zero(::Type{Point{T}}) where {T} = Point(zero(T), zero(T))
+(A::Point{T}, B::Point{T}) where {T} = Point(A.x + B.x, A.y + B.y)
-(A::Point{T}, B::Point{T}) where {T} = Point(A.x - B.x, A.y - B.y)

Base.zero(::Type{Rectangle{T}}) where {T} = zero(Point{T})
+(rect::Rectangle{T}, pnt::Point{T}) where {T} = apply(Fix2(+, pnt), rect)
-(rect::Rectangle{T}, pnt::Point{T}) where {T} = apply(Fix2(-, pnt), rect)
-(pnt::Point{T}, rect::Rectangle{T}) where {T} = apply(Fix1(-, pnt), rect)

Base.zero(::Type{BoundingBox{T}}) where {T} = zero(Point{T})
+(box::BoundingBox{T}, pnt::Point{T}) where {T} = apply(Fix2(+, pnt), box; swap = false)
-(box::BoundingBox{T}, pnt::Point{T}) where {T} = apply(Fix2(-, pnt), box)
-(pnt::Point{T}, box::BoundingBox{T}) where {T} = apply(Fix1(-, pnt), box; swap = true)

# Addition and subtraction of bounding-boxes (following the rules of the
# addition and subtraction of sets) and corresponding addtive identity.
function +(A::BoundingBox, B::BoundingBox)
    T = promote(coord_type(A), coord_type(B))
    (isempty(A) || isempty(B)) && return BoundingBox{T}(nothing)
    return BoundingBox{T}(first(A) + first(B), last(A) + last(B))
end
function -(A::BoundingBox, B::BoundingBox)
    T = promote(coord_type(A), coord_type(B))
    (isempty(A) || isempty(B)) && return BoundingBox{T}(nothing)
    return BoundingBox{T}(first(A) - last(B), last(A) - first(B))
end

# Add or remove a margin δ to a bounding-box box.
+(box::BoundingBox, δ::Number) = grow(box, δ)
-(box::BoundingBox, δ::Number) = shrink(box, δ)

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

"""
    round([T,] obj, [r::RoundingMode])

applies the `round` function to the coordinates of the vertices defining the
vertex-based geometric object `obj` and returns an object of the same kind
built with the resulting vertices. Optional argument `T` can be the type of the
returned object or the type of the coordinates for the returned object.
Rounding mode may be specified by optional argument `r`, the default being the
same as the `round` method for a scalar value.

See also: [`floor`](@ref floor(::TwoDimensional.Point)), [`ceil`](@ref
ceil(::TwoDimensional.Point)).

"""
Base.round(obj::VertexBasedObject) = apply(round, obj)
Base.round(obj::VertexBasedObject, r::RoundingMode) = apply(Round{Nothing}(r), obj)
Base.round(::Type{T}, obj::VertexBasedObject) where {T} = apply(Round{T}(), obj)
Base.round(::Type{T}, obj::VertexBasedObject, r::RoundingMode) where {T} = apply(Round{T}(r), obj)

for type in (:Point, :Rectangle, :BoundingBox)
    @eval begin
        Base.round(::Type{$type}, obj::$type) = round(obj)
        Base.round(::Type{$type}, obj::$type, r::RoundingMode) = round(obj, r)
        Base.round(::Type{$type{T}}, obj::$type) where {T} = round(T, obj)
        Base.round(::Type{$type{T}}, obj::$type, r::RoundingMode) where {T} = round(T, obj, r)
    end
end

"""
    floor([T,] obj)

applies the `floor` function to the coordinates of the vertices defining the
vertex-based geometric object `obj` and returns an object of the same kind
built with the resulting vertices. Optional argument `T` can be the type of the
returned object or the type of the coordinates for the returned object.

See also: [`round`](@ref round(::TwoDimensional.Point)), [`ceil`](@ref
ceil(::TwoDimensional.Point)).

For bounding-boxes, see also: [`interior`](@ref TwoDimensional.interior),
[`exterior`](@ref TwoDimensional.exterior).

"""
Base.floor(obj::VertexBasedObject) = apply(floor, obj)

"""
    ceil([T,] obj)

applies the `ceil` function to the coordinates of the vertices defining the
vertex-based geometric object `obj` and returns an object of the same kind
built with the resulting vertices. Optional argument `T` can be the type of the
returned object or the type of the coordinates for the returned object.

See also: [`round`](@ref round(::TwoDimensional.Point)), [`floor`](@ref
floor(::TwoDimensional.Point)).

For bounding-boxes, see also: [`interior`](@ref TwoDimensional.interior),
[`exterior`](@ref TwoDimensional.exterior).

"""
Base.ceil(obj::VertexBasedObject) = apply(ceil, obj)

for type in (:Point, :Rectangle, :BoundingBox), func in (:floor, :ceil)
    @eval begin
        Base.$func(::Type{$type}, obj::$type) = $func(obj)
        Base.$func(::Type{$type{T}}, obj::$type) where {T} = $func(T, obj)
        Base.$func(::Type{T}, obj::$type) where {T} = apply(Fix1($func, T), obj)
    end
end

# `min()`, `max()`, and `minmax()` for points work as for Cartesian indices.
Base.min(a::Point{T}, b::Point{T}) where {T} = Point{T}(min(a.x, b.x), min(a.y, b.y))
Base.max(a::Point{T}, b::Point{T}) where {T} = Point{T}(max(a.x, b.x), max(a.y, b.y))
Base.minmax(a::Point{T}, b::Point{T}) where {T} = parts(Rectangle(a, b))
for func in (:min, :max, :minmax)
    @eval begin
        Base.$func(a::Point, b::Point) = $func(promote_coord_type(a, b)...)
    end
end

Base.clamp(pnt::Point, box::BoundingBox) = clamp(promote_coord_type(pnt, box)...)
function Base.clamp(pnt::Point{T}, box::BoundingBox{T}) where {T}
    x = clamp(pnt.x, box.xmin, box.xmax)
    y = clamp(pnt.y, box.ymin, box.ymax)
    return Point{T}(x, y)
end

# Mathematic functions for points.
Base.hypot(pnt::Point) = hypot(pnt.x, pnt.y)
Base.abs(pnt::Point) = hypot(pnt)
Base.abs2(pnt::Point) = abs2(pnt.x) + abs2(pnt.y)
LinearAlgebra.norm(pnt::Point) = abs(pnt)
Base.Math.atan(pnt::Point) = atan(pnt.y, pnt.x)

"""
    dot(a::Point, b::Point)
    a ⋅ b

yields the dot product (a.k.a. scalar product or inner product) of the two
points `a` and `b` which is given by:

    a.x*b.x + a.y*b.y

"""
LinearAlgebra.dot(a::Point, b::Point) = a.x*b.x + a.y*b.y

"""
    cross(a::Point, b::Point)
    a * b

yields the cross product (a.k.a. vector product or directed area product) of
the two points `a` and `b` which is given by:

    a.x*b.y - a.y*b.x

"""
LinearAlgebra.cross(a::Point, b::Point) = a.x*b.y - a.y*b.x
*(a::Point, b::Point) = cross(a, b)

"""
    distance(A, B)

yields the Euclidean distance between the 2 points `A` and `B`.

"""
distance(A::Point, B::Point) = hypot(A - B)

"""
    area(rect)

yields the area of the rectangle `rect`.

"""
area(rect::Rectangle) = (rect.x1 -rect.x0)*(rect.y1 - rect.y0)

"""
    center(obj::GeometricObject) -> c::Point

yields the central point of the geometric object `obj`.

"""
center(pnt::AbstractPoint) = Point(pnt.x, pnt.y)
center(pnt::Point) = pnt
center(rect::Rectangle) = Point((rect.x0 + rect.y0)/2, (rect.y0 + rect.y1)/2)
center(box::BoundingBox) =
    isempty(box) ? throw(ArgumentError("cannot get center of empty box")) :
    Point((box.xmin + box.xmax)/2, (box.ymin + box.ymax)/2)

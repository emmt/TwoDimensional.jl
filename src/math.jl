# Check for equality.
==(a::Point, b::Point) =
    (a.x == b.x) & (a.y == b.y)
==(a::Rectangle, b::Rectangle) =
    (first(a) == first(b)) & (last(a) == last(b))
==(a::Circle, b::Circle) =
    radius(a) == radius(b) && center(a) == center(b)
==(a::BoundingBox, b::BoundingBox) =
    isempty(a) ? isempty(b) : ((first(a) == first(b)) & (last(a) == last(b)))
==(a::Polygon, b::Polygon) = begin
    a === b && return true
    eachindex(b) == eachindex(a) || return false
    @inbounds for i in eachindex(a)
        a[i] == b[i] || return false
    end
    return true
end

# Check for approximate equality.
Base.isapprox(a::Point, b::Point; kwds...) =
    isapprox(a.x, b.x; kwds...) && isapprox(a.y, b.y; kwds...)
Base.isapprox(a::Rectangle, b::Rectangle; kwds...) =
    isapprox(first(a), first(b); kwds...) && isapprox(last(a), last(b); kwds...)
Base.isapprox(a::Circle, b::Circle; kwds...) =
    isapprox(radius(a), radius(b); kwds...) && isapprox(center(a), center(b); kwds...)
Base.isapprox(a::Polygon, b::Polygon; kwds...) =
    elements(a) === elements(b) || isapprox(elements(a), elements(b); kwds...)
Base.isapprox(a::BoundingBox, b::BoundingBox; kwds...) =
    isapprox(first(a), first(b); kwds...) && isapprox(last(a), last(b); kwds...)

# Extend ∈ operator for points.
Base.in(pnt::AbstractPoint, msk::MaskElement) = pnt ∈ shape(msk)

Base.in(pnt::AbstractPoint, box::BoundingBox) =
    (box.xmin ≤ pnt.x ≤ box.xmax) & (box.ymin ≤ pnt.y ≤ box.ymax)

Base.in(pnt::AbstractPoint, rect::Rectangle) =
    (rect.x0 ≤ pnt.x ≤ rect.x1) & (rect.y0 ≤ pnt.y ≤ rect.y1)

Base.in(pnt::AbstractPoint, circ::Circle) =
    distance(pnt - center(circ)) ≤ radius(circ)

Base.in(pnt::Point, poly::Polygon) =
    # FIXME: Converting to the same coordinate type has some cost...
    in(promote_coord_type(pnt, poly)...)
Base.in(pnt::Point{T}, poly::Polygon{T}) where {T} =
    winding_number_test(pnt, elements(poly))

# Extend ⊆ operator.
Base.issubset(obj::GeometricObject, msk::MaskElement) = obj ⊆ shape(msk)

function Base.issubset(rect::Rectangle, circ::Circle)
    # Return whether the corner the most distant form the center is inside the
    # circle.
    xmin, ymin = first(rect) - center(circ)
    xmax, ymax = last(rect) - center(circ)
    return hypot(max(-xmin, xmax), max(-ymin, ymax)) ≤ radius(circ)
end

Base.issubset(rect::Rectangle, box::BoundingBox) =
    (first(rect) ∈ box) && (last(rect) ∈ box)

function Base.issubset(box::BoundingBox, circ::Circle)
    # Return whether the corner the most distant form the center is inside the
    # circle.
    isempty(box) && return true
    xmin, ymin = first(box) - center(circ)
    xmax, ymax = last(box) - center(circ)
    return hypot(max(-xmin, xmax), max(-ymin, ymax)) ≤ radius(circ)
end

Base.issubset(A::Union{Circle,Polygon}, B::Union{Rectangle,BoundingBox}) =
    BoundingBox(A) ⊆ B # FIXME optimize as BoundingBox(A) is never empty

function Base.issubset(poly::Polygon, circ::Circle)
    # Return whether al vertices are inside the circle.
    for i in eachindex(poly)
        poly[i] ∈ circ || return false
    end
    return true
end

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

# Unary plus does nothing.
+(a::GeometricObjectLike) = a

# Unary minus negates coordinates and should behave as multiplying by -1.
-(a::GeometricObjectLike) = apply(-, a)
-(a::BoundingBox) = apply(-, a; swap = true)
-(a::Circle) = Circle(-center(a), radius(a))

# Scaling of geometric objects and corresponding multiplicative identity.
Base.one(obj::GeometricObjectLike) = one(typeof(obj))
Base.one(::Type{<:GeometricObjectLike{T}}) where {T} = one(T)

# NOTE Even though `T` seems to be unnecessary, do not remove `{T}... where {T}` for
#      correct dispatching on Julia ≤ 1.4
*(a::GeometricObjectLike, b::Number) = b*a
*(a::Number, b::GeometricObjectLike{T}) where {T} = apply(Fix1(*, a), b)
*(a::Number, b::BoundingBox{T}        ) where {T} = apply(Fix1(*, a), b; swap = a < zero(a))
*(a::Number, b::Circle{T}             ) where {T} = Circle(a*center(b), abs(a)*radius(b))

# NOTE Even though `T` seems to be unnecessary, do not remove `{T}... where {T}` for
#      correct dispatching on Julia ≤ 1.4
/(a::GeometricObjectLike, b::Number) = b\a
\(a::Number, b::GeometricObjectLike{T}) where {T} = apply(Fix1(\, a), b)
\(a::Number, b::BoundingBox{T}        ) where {T} = apply(Fix1(\, a), b; swap = a < zero(a))
\(a::Number, b::Circle{T}             ) where {T} = Circle(a\center(b), abs(a)\radius(b))

# Move a geometric object by adding or subtracting a point and corresponding addtive
# identity.
Base.zero(obj::GeometricObjectLike) = zero(typeof(obj))
Base.zero(::Type{<:GeometricObjectLike{T}}) where {T} = Point(zero(T), zero(T))

+(a::Point, b::Point) = +(promote_coord_type(a, b)...)
+(a::Point{T}, b::Point{T}) where {T} = Point(a.x + b.x, a.y + b.y)

-(a::Point, b::Point) = -(promote_coord_type(a, b)...)
-(a::Point{T}, b::Point{T}) where {T} = Point(a.x - b.x, a.y - b.y)

+(a::Point, b::GeometricObjectLike) = b + a

+(a::GeometricObjectLike,    b::Point) = +(promote_coord_type(a, b)...)
+(a::GeometricObjectLike{T}, b::Point{T}) where {T} = apply(Fix2(+, b), a)
+(a::BoundingBox{T},         b::Point{T}) where {T} = apply(Fix2(+, b), a; swap = false)
+(a::Circle{T},              b::Point{T}) where {T} = Circle(center(a) + b, radius(a))

-(a::GeometricObjectLike,    b::Point) = -(promote_coord_type(a, b)...)
-(a::GeometricObjectLike{T}, b::Point{T}) where {T} = apply(Fix2(-, b), a)
-(a::BoundingBox{T},         b::Point{T}) where {T} = apply(Fix2(-, b), a; swap = false)
-(a::Circle{T},              b::Point{T}) where {T} = Circle(center(a) - b, radius(a))

-(a::Point,    b::GeometricObjectLike) = -(promote_coord_type(a, b)...)
-(a::Point{T}, b::GeometricObjectLike{T}) where {T} = apply(Fix1(-, a), b)
-(a::Point{T}, b::BoundingBox{T}        ) where {T} = apply(Fix1(-, a), b; swap = true)
-(a::Point{T}, b::Circle{T}             ) where {T} = Circle(a - center(b), radius(b))

# Performing a geometric operation on a mask amounts to performing the
# operation on the embedded shape. Inverting a mask toggles its opacity.
Base.inv(msk::MaskElement) = MaskElement(shape(msk); opaque = !is_opaque(msk))
Base.inv(msk::Mask) = apply(inv, msk)

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

# Apply affine transform.
(A::AffineTransform)(msk::MaskElement) = MaskElement(A(shape(mask)); opaque = is_opaque(msk))
(A::AffineTransform)(rect::Rectangle) = A(Polygon(rect))
(A::AffineTransform)(circ::Circle) = error("not yet implemented")
(A::AffineTransform)(poly::Polygon) = apply(A, poly)
*(A::AffineTransform, obj::GeometricObject) = A(obj)

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

for type in VERTEX_BASED_TYPES
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

for type in VERTEX_BASED_TYPES, func in (:floor, :ceil)
    @eval begin
        Base.$func(::Type{$type}, obj::$type) = $func(obj)
        Base.$func(::Type{$type{T}}, obj::$type) where {T} = $func(T, obj)
        Base.$func(::Type{T}, obj::$type) where {T} = apply(Fix1($func, T), obj)
    end
end

# `min()`, `max()`, and `minmax()` for points work as for Cartesian indices.
Base.min(a::Point{T}, b::Point{T}) where {T} = Point{T}(min(a.x, b.x), min(a.y, b.y))
Base.max(a::Point{T}, b::Point{T}) where {T} = Point{T}(max(a.x, b.x), max(a.y, b.y))
Base.minmax(a::Point{T}, b::Point{T}) where {T} = elements(Rectangle(a, b))
for func in (:min, :max, :minmax)
    @eval begin
        Base.$func(a::Point, b::Point) = $func(promote_coord_type(a, b)...)
    end
end

Base.clamp(pnt::Point, box::BoundingBox) = clamp(promote_coord_type(pnt, box)...)
function Base.clamp(pnt::Point{T}, box::BoundingBox{T}) where {T}
    isempty(box) && throw(ArgumentError("cannot clamp to empty bounding-box"))
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
    TwoDimensional.area(obj)

yields the area of the geometric object `obj`.

"""
area(msk::MaskElement) = area(shape(msk))
area(pnt::AbstractPoint) = (z = zero(coord_type(pnt)); return z*z) # NOTE coord. may have units
area(rect::Rectangle) = (rect.x1 - rect.x0)*(rect.y1 - rect.y0)
area(circ::Circle) = (r = radius(circ); return (π*r)*r) # NOTE take care of correct conversion

"""
    center(obj::GeometricObject) -> c::Point

yields the central point of the geometric object `obj`. For a polygon, the
center of gravity of the vertices is returned.

"""
center(msk::MaskElement) = center(shape(msk))
center(pnt::AbstractPoint) = Point(pnt.x, pnt.y)
center(pnt::Point) = pnt
center(rect::Rectangle) = Point((rect.x0 + rect.y0)/2, (rect.y0 + rect.y1)/2)
center(circ::Circle) = getfield(circ, 1)
center(box::BoundingBox) =
    isempty(box) ? throw(ArgumentError("cannot get center of empty box")) :
    Point((box.xmin + box.xmax)/2, (box.ymin + box.ymax)/2)
function center(poly::Polygon)
    i, i_last = firstindex(poly), lastindex(poly)
    pnt = float(poly[i])
    @inbounds while i < i_last
        i += 1
        pnt += float(poly[i])
    end
    return pnt/length(poly)
end

"""
    TwoDimensional.radius(obj::GeometricObject)

yields the radius of the geometric object `obj`. The result is the radius of
the [smallest circle enclosing the
object](https://en.wikipedia.org/wiki/Smallest-circle_problem).

For circle-like and point-like objects with integer coordinate type, the radius
is also integer. For all other geometric objects, the raius is floating-point.

"""
radius(msk::MaskElement) = radius(shape(msk))
radius(pnt::AbstractPoint) = zero(coord_type(pnt))
radius(rect::Rectangle) = half(diameter(rect))
radius(circ::Circle) = getfield(circ, 2)
radius(poly::Polygon) = half(diameter(poly))
radius(box::BoundingBox) = half(diameter(box))

"""
    TwoDimensional.diameter(obj::GeometricObject)

yields the diameter of the geometric object `obj`. The result is the diameter
of the [smallest circle enclosing the
object](https://en.wikipedia.org/wiki/Smallest-circle_problem).

For circle-like and point-like objects with integer coordinate type, the
diameter is also integer. For all other geometric objects, the raius is
floating-point.

"""
diameter(msk::MaskElement) = diameter(shape(msk))
diameter(pnt::AbstractPoint) = zero(coord_type(pnt))
diameter(rect::Rectangle) = distance(first(rect), last(rect))
diameter(circ::Circle) = twice(radius(circ))
diameter(box::BoundingBox) =
    isempty(box) ? zero(float(coord_type(box))) : distance(first(box), last(box))
diameter(poly::Polygon) = error("not yet implemented")

"""
    TwoDimensional.is_convex(obj) -> bool

yields whether the geometric object `obj` is convex.

"""
is_convex(msk::MaskElement) = is_convex(shape(msk))
is_convex(pnt::Point) = true
is_convex(circ::Circle) = true
is_convex(rect::Rectangle) = true
is_convex(box::BoundingBox) = !isempty(box)
is_convex(poly::Polygon) = geometric_properties(poly).convex

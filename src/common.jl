"""
    TwoDimensional.elements(obj::GeometricElement)

yields the individual elements of elementary geometric `obj` from which it can be re-built
without ambiguities. For example, for a point `pnt`:

    Point(elements(pnt)...) === pnt
    Point(elements(pnt)) === pnt

both hold.

Geometrical objects that have homogeneous elements (see
[`TwoDimensional.VertexBasedObject`](@ref)) extend the `Base.Tuple` method to return these
elements, the `Base.getindex` method to directly index among these elements, and the
[`TwoDimensional.apply`](@ref) method.

"""
elements(pnt::Point) = getfield(pnt, 1)
elements(rect::Rectangle) = getfield(rect, 1)
elements(poly::Polygon) = getfield(poly, 1)
elements(circ::Circle) = (center(circ), radius(circ))
elements(box::BoundingBox) = getfield(box, 1)
elements(elem::MaskElement) = elements(shape(elem))
elements(msk::Mask) = getfield(msk, 1)

# Extend methods `Base.Tuple` and `Base.getindex` for some geometric objects.
Base.Tuple(obj::Union{Point,Rectangle,Circle,BoundingBox}) = elements(obj)
@inline @propagate_inbounds Base.getindex(obj::Union{Point,Rectangle,BoundingBox}, i::Integer) =
    getindex(elements(obj), as(Int, i))

Base.Tuple(obj::Polygon) = Tuple(elements(obj))
Base.Tuple(obj::Mask) = Tuple(elements(obj))

"""
    TwoDimensional.apply(f, obj)

applies function `f` to each part of geometric object `obj` and rebuild an object of the
same kind with the result. Here `f` is supposed to be a function implementing an
elementary geometric operation such as moving, scaling, etc. the geometric object `obj`.

If `obj` is a bounding-box, keyword, `swap` (default `false`) specifies whether to swap
the first and last end-points of the box.

See also [`TwoDimensional.elements`](@ref) and [`TwoDimensional.VertexBasedObject`](@ref).

"""
apply(f, pnt::Point) = Point(f(pnt[1]), f(pnt[2]))
apply(f, rect::Rectangle) = Rectangle(f(rect[1]), f(rect[2]))
apply(f, poly::Polygon) = Polygon(map(f, elements(poly)))
apply(f, box::BoundingBox; swap::Bool = false) = BoundingBox(f(box[swap ? 2 : 1]),
                                                             f(box[swap ? 1 : 2]))
apply(f, elem::MaskElement) = MaskElement(f(shape(elem)); opaque = is_opaque(elem))
apply(f, msk::Mask) = Mask(map(f, elements(msk)))

# Swap two elements.
swap((x, y)::NTuple{2,Any}) = (y, x)
swap(x, y) = (y, x)

twice(x) = x + x
half(x) = x/twice(one(x))

"""
    TwoDimensional.nearest(T, x) -> y::T

yields the nearest value of type `T` to `x`.

"""
nearest(::Type{T}, x::T) where {T<:Number} = x
nearest(::Type{T}, x::Number) where {T<:Number} =
    ((real_type(T) <: Integer) & !(real_type(x) <: Integer)) ? round(T, x) : as(T, x)

"""
    TwoDimensional.shape(obj)

yields the elementary geometric object defining the shape of `obj`. The result is `obj`
itself if it is an elementary geometric object.

"""
shape(obj::MaskElement) = obj.shape
shape(obj::ShapeElement) = obj

"""
    TwoDimensional.vertices(obj)

yields the vertices defining the vertex-based graphical object `obj`. The result is a
tuple or a vector of points.

"""
vertices(msk::Union{RectangularMask,PolygonalMask}) = vertices(shape(msk))
vertices(poly::Polygon) = elements(poly)
vertices(pnt::Point) = (pnt,)
function vertices(rect::Rectangle)
    (x0, y0), (x1, y1) = rect
    return (Point(x0, y0), Point(x1, y0), Point(x1, y1), Point(x0, y1))
end
function vertices(box::BoundingBox)
    (xmin, ymin), (xmax, ymax) = box
    return (Point(xmin, ymin), Point(xmax, ymin), Point(xmax, ymax), Point(xmin, ymax))
end

# Fast versions of `min` and `max` which return their first argument if any
# argument is a NaN. The fast version of `minmax` cannot warrant this
# property.
fastmin(x::T, y::T) where {T} = (x > y ? y : x)
fastmax(x::T, y::T) where {T} = (x < y ? y : x)
fastminmax(x::T, y::T) where {T} = (x > y ? (y,x) : (x,y))
for func in (:fastmax, :fasmin, :fastminmax)
    @eval $func(x, y) = $func(promote(x, y)...)
end

# Promotion rules.
for type in (:Point, :Rectangle, :Circle, :BoundingBox)
    @eval begin
        Base.promote_rule(::Type{$type{T}}, ::Type{$type{T}}) where {T} = $type{T}
        Base.promote_rule(::Type{$type{T}}, ::Type{$type{S}}) where {T,S} =
            $type{promote_type(T,S)}
    end
end
Base.promote_rule(::Type{Polygon{T,V}}, ::Type{Polygon{T,V}}) where {T,V} = Polygon{T,V}
function Base.promote_rule(::Type{Polygon{T,Vector{Point{T}}}},
                           ::Type{Polygon{S,Vector{Point{S}}}}) where {T,S}
    R = promote_type(T, S)
    return Polygon{R,Vector{Point{R}}}
end

Base.show(io::IO, ::MIME"text/plain", obj::GeometricObject) = show(io, obj)

"""
    coord_type(A) -> T

yields the coordinate type of a geometrical object or of its type. With several arguments:

    coord_type(A...) -> T

or if single argument a tuple or a vector of geometrical objects, the promoted coordinate
type is returned.

"""
coord_type(::GeometricObjectLike{T}) where {T} = T
coord_type(::Type{<:GeometricObjectLike{T}}) where {T} = T
coord_type(::Type{<:GeometricObjectLike}) = Any

# Methods for zero or more than one argument.
coord_type() = Union{}
@inline coord_type(A, B...) = coord_type((A, B...))

# Method for tuples or vectors of geometric objects.
coord_type(A::Tuple) = to_same_concrete_type(map(coord_type, A)...)
coord_type(A::AbstractVector) = mapreduce(coord_type, to_same_concrete_type, A; init=coord_type())

# Optimizations for multiple objects/types.
coord_type(A::Vararg{GeometricObjectLike{T}}) where {T} = T
coord_type(A::Tuple{Vararg{GeometricObjectLike{T}}}) where {T} = T
coord_type(A::AbstractVector{<:GeometricObjectLike{T}}) where {T} = T

# Fallbacks for any objects and for errors.
coord_type(A::Any) = coord_type(typeof(A))
@noinline coord_type(T::DataType) =
    throw(ArgumentError("objects of type `$T` have no definite coordinate type"))

"""
    promote_coord_type(objs...)

converts all arguments `objs...` to a common coordinate type and return them as a tuple.

"""
promote_coord_type(obj::GeometricObjectLike) = obj
@inline promote_coord_type(objs::GeometricObjectLike...) =
    convert_coord_type(coord_type(objs), objs...)

"""
    convert_coord_type(T::Type, obj::GeometricObject)
    GeometricObject{T}(obj)

convert the coordinate type of a geometrical object `obj` to `T`. If the coordinate type
of `obj` is already `T`, `obj` itself is returned.

"""
convert_coord_type(::Type{T}, obj::GeometricObject{T}) where {T} = obj
for type in (:Point, :Rectangle, :Circle, :Polygon, :BoundingBox, :MaskElement)
    @eval begin
        convert_coord_type(::Type{T}, obj::$type{T}) where {T} = obj
        convert_coord_type(::Type{T}, obj::$type) where {T} = $type{T}(obj)
    end
end
for type in (:GeometricObject, :GeometricElement, :ShapeElement, :MaskElement,
             :AbstractPoint, :Point, :Rectangle, :Circle, :BoundingBox, :Mask)
    @eval begin
        convert_coord_type(::Type{T}, ::Type{<:$type}) where {T} = $type{T}
    end
end
convert_coord_type(::Type{T}, ::Type{<:Polygon}) where {T} = Polygon{T,Vector{Point{T}}}
convert_coord_type(::Type{T}, ::Type{Polygon{T,V}}) where {T,V} = Polygon{T,V}

convert_coord_type(::Type{T}, msk::Mask{T}) where {T} = msk
convert_coord_type(::Type{T}, msk::Mask) where {T} = Mask{T}(elements(msk))

GeometricObject(obj::GeometricObject) = obj
GeometricObject{T}(obj::GeometricObject{T}) where {T} = obj
GeometricObject{T}(obj::GeometricObject) where {T} = convert_coord_type(T, obj)


Base.float(obj::GeometricObject{T}) where {T} = convert_coord_type(float(T), obj)
Base.float(::Type{G}) where {T,G<:GeometricObject{T}} = convert_coord_type(float(T), G)

for (func, type)  in ((:convert_bare_type,           Real),
                      (:convert_real_type,           Real),
                      (:convert_floating_point_type, AbstractFloat),)
    @eval begin
        TypeUtils.$func(::Type{T}, obj::GeometricObject) where {T<:$type} =
            convert_coord_type($func(T, coord_type(obj)), obj)
        TypeUtils.$func(::Type{T}, ::Type{G}) where {T<:$type,G<:GeometricObject} =
            convert_coord_type($func(T, coord_type(G)), G)
    end
end

"""
    convert_coord_type(T::Type, objs::GeometricObject...) -> objs′

yields the tuple of geometric objects `objs...` converted to the coordinate
type `T`.

"""
convert_coord_type(::Type{T}, objs::GeometricObject{T}...) where {T} = objs
convert_coord_type(::Type{T}, objs::GeometricObject...) where {T} =
    map(Fix1(convert_coord_type, T), objs)

"""
    TwoDimensional.is_nothing(x)
    TwoDimensional.is_nothing(x, y...)

yield whether `x` is `nothing` and, if other arguments `y...` are specified,
that all other `y...` are `nothing`.

See also [`TwoDimensional.is_nothing`](@ref).

"""
is_nothing(x::Nothing) = true
is_nothing(x::Any) = false
@inline is_nothing(x, y...) = is_nothing(x) && is_nothing(y...)

"""
    TwoDimensional.is_something(x)
    TwoDimensional.is_something(x, y...)

yield whether `x` is not `nothing` and, if other arguments `y...` are
specified, that none of the other `y...` is `nothing`.

See also [`TwoDimensional.is_nothing`](@ref).

"""
is_something(x) = !is_nothing(x)
@inline is_something(x, y...) = is_something(x) && is_something(y...)

"""
    Unsupported(T::DataType...)

yields an union of types `T...` and of type `Unsupported`. Such an union can be used to
mark unsupported argument type(s) and yet define a method applicable with that type(s)
(presumably a method that throws an instructive error) and which can be extended later
with the same signature except that with `Unsupported(T...)` replaced by `T...`. This
trick avoids conflicts that prevent pre-compilation with package extensions.

For example, in the main package:

```julia
some_method(some_arg::Unsupported{SomeType}) =
    error("package `SomeOtherPackage` has not yet been loaded")
```

and in the extension (e.g. automatically loaded when using `SomeOtherPackage`):

```julia
some_method(some_arg::SomeType) = do_something_with(some_arg)
```

"""
struct Unsupported
    Unsupported() = error("it is not possible to instanciate this type")
end
Unsupported(T::DataType...) = Union{T...,Unsupported}

new_array(::Type{T}, inds::eltype(ArrayShape)...) where {T} = new_array(T, inds)
new_array(::Type{T}, inds::ArrayShape{N}) where {T,N} = new_array(T, as_array_shape(inds))
new_array(::Type{T}, rngs::NTuple{N,Base.OneTo}) where {T,N} = new_array(T, as_array_size(rngs))
new_array(::Type{T}, dims::Dims{N}) where {T,N} = Array{T,N}(undef, dims)
new_array(::Type{T}, rngs::Unsupported(ArrayAxes{N})) where {T,N} =
    error("package `OffsetArrays` must be loaded for such array index ranges")

"""
    TwoDimensional.parts(obj::GeometricElement)

yields the individual elements of elementary geometric `obj` from which it can
be re-built without ambiguities. For example, for a point `pnt`:

    Point(parts(pnt)...) === pnt
    Point(parts(pnt)) === pnt

both hold.

Geometrical objects that have homogeneous parts (see
['TwoDimensional.VertexBasedObject`](@ref)) extend the `Base.Tuple` method to
return these parts, the `Base.getindex` method to directly index among these
parts, and the ['TwoDimensional.apply`](@ref) method.

"""
parts(pnt::Point) = getfield(pnt, 1)
parts(rect::Rectangle) = getfield(rect, 1)
parts(poly::Polygon) = vec(poly)
parts(circ::Circle) = (center(circ), radius(circ))
parts(box::BoundingBox) = getfield(box, 1)

# Extend methods `Base.Tuple` and `Base.getindex` for some geometric objects.
Base.Tuple(obj::Union{Point,Rectangle,Circle,BoundingBox}) = parts(obj)
@inline @propagate_inbounds Base.getindex(obj::Union{Point,Rectangle,BoundingBox}, i::Integer) =
    getindex(parts(obj), as(Int, i))

"""
    TwoDimensional.apply(f, obj)

applies function `f` to each part of geometric object `obj` and rebuild an
object of the same king with the result.

If `obj` is a bounding-box, keyword, `swap` (default `false`) specifies whether
to swap the first and las end-points of the box.

See also ['TwoDimensional.parts`](@ref) and
['TwoDimensional.VertexBasedObject`](@ref).

"""
@inline apply(f, pnt::Point) = Point(f(pnt[1]), f(pnt[2]))
@inline apply(f, rect::Rectangle) = Rectangle(f(rect[1]), f(rect[2]))
@inline apply(f, poly::Polygon) = Polygon(map(f, vec(poly)))
@inline apply(f, box::BoundingBox; swap::Bool = false) =
    BoundingBox(f(box[swap ? 2 : 1]), f(box[swap ? 1 : 2]))

# Swap two elements.
swap((x, y)::NTuple{2,Any}) = (y, x)
swap(x, y) = (y, x)

twice(x) = x + x
half(x) = x/twice(one(x))

"""
    TwoDimensional.vertices(obj)

yields the vertices defining the vertex-based graphical object `obj`. The
result is a tuple or a vector of points.

"""
vertices(poly::Polygon) = vec(poly)
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

# FIXME: use TypeUtils
as(::Type{T}, x::T) where {T} = x
as(::Type{T}, x) where {T} = convert(T, x)::T

# Basic constructors (with tuple argument) for points and bounding boxes and
# API to make them indexable iterators.
for type in (:Point, :Rectangle, :BoundingBox)
    @eval begin
        Base.promote_type(::Type{$type{T}}, ::Type{$type{T}}) where {T} = $type{T}
        Base.promote_type(::Type{$type{T₁}}, ::Type{$type{T₂}}) where {T₁,T₂} =
            $type{promote_type(T₁,T₂)}
    end
end

Base.show(io::IO, ::MIME"text/plain", obj::GeometricObject) = show(io, obj)

"""
    coord_type(obj) -> T

yields the coordinate type of a geometrical object or of its type.

"""
coord_type(::GeometricObject{T}) where {T} = T
coord_type(::Type{<:GeometricObject{T}}) where {T} = T

"""
    promote_coord_type(types::Type{<:GeometricObject}...) -> T

yields the promoted coordinate type of geometric object types in `types...`

"""
promote_coord_type(::Type{<:GeometricObject{T}}) where {T} = T
@inline function promote_coord_type(::Type{<:GeometricObject{T}},
                                    types::Type{<:GeometricObject}...) where {T}
    return _promote_coord_type(T, types...)
end

@inline function _promote_coord_type(::Type{R}, ::Type{<:GeometricObject{T}},
                                     types::Type{<:GeometricObject}...) where {R,T}
    return _promote_coord_type(promote_type(R, T), types...)
end
_promote_coord_type(::Type{R}) where {R} = R # end the recursion

"""
    promote_coord_type(objs::GeometricObject...) -> objs′

yields the tuple of geometric objects `objs...` converted to the same
coordinate type.

"""
promote_coord_type(obj::GeometricObject) = obj
@inline promote_coord_type(objs::GeometricObject...) =
    convert_coord_type(promote_coord_type(map(typeof, objs)...), objs...)

"""
    convert_coord_type(T::Type, obj::GeometricObject)

converts the coordinate type of a geometrical object `obj` to `T`. If the
coordinate type of `obj` is already `T`, `obj` itself is returned.

"""
convert_coord_type(::Type{T}, obj::GeometricObject{T}) where {T} = obj
for type in (:Point, :Rectangle, :Circle, :Polygon, :BoundingBox)
    @eval begin
        convert_coord_type(::Type{T}, obj::$type{T}) where {T} = obj
        convert_coord_type(::Type{T}, obj::$type) where {T} = $type{T}(obj)
    end
end
for type in (:GeometricObject, :GeometricElement, :ShapeElement, :MaskElement,
             :AbstractPoint, :Point, :Rectangle, :Circle, :BoundingBox,)
    @eval begin
        convert_coord_type(::Type{T}, ::Type{<:$type}) where {T} = $type{T}
    end
end
convert_coord_type(::Type{T}, ::Type{<:Polygon}) where {T} = Polygon{T,Vector{Point{T}}}
convert_coord_type(::Type{T}, ::Type{Polygon{T,V}}) where {T,V} = Polygon{T,V}

Base.float(obj::GeometricObject{T}) where {T} = convert_coord_type(float(T), obj)
Base.float(::Type{G}) where {T,G<:GeometricObject{T}} = convert_coord_type(float(T), G)

for (func, type)  in ((:convert_bare_type,           Real),
                      (:convert_real_type,           Real),
                      (:convert_floating_point_type, AbstractFloat),)
    @eval begin
        Unitless.$func(::Type{T}, obj::GeometricObject) where {T<:$type} =
            convert_coord_type($func(T, coord_type(obj)), obj)
        Unitless.$func(::Type{T}, ::Type{G}) where {T<:$type,G<:GeometricObject} =
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

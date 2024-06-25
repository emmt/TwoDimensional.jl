"""
    poly = Polygon{T}(pnts...)
    poly = Polygon{T}(pnts)

construct a polygon with vertices given by point-like objects `pnts...` or vector or tuple
of point-like objects `pnts`. Parameter `T` is the type used to store coordinates, if
omitted, a common coordinate type is automatically inferred.

The vertices of a polygon may be stored as a tuple or as a vector. Call `values(poly)` to
get the object backing the storage of the vertices of the polygon `poly`. The method
`Tuple(poly)` yields a tuple of the vertices of `poly`. The method `vec(poly)` yields a
vector of the vertices of the polygon which may be shared with `poly`. Call
`collect(poly)` to make an independent copy of the vector of vertices.

Vertices are directly accessible by indexing the polygon object:

    length(poly)     # number of vertices
    firstindex(poly) # index of first vertex
    lastindex(poly)  # index of last vertex
    poly[i]          # i-th vertex (a point)
    first(poly)      # first vertex
    last(poly)       # last vertex
    values(poly)     # vertices stored by `poly`
    vec(poly)        # vector of vertices
    Vector(poly)     # vector of vertices
    Tuple(poly)      # tuple of vertices

See also [`Point`](@ref TwoDimensional.Point), [`TwoDimensional.vertices`](@ref
TwoDimensional.vertices), and [`BoundingBox`](@ref TwoDimensional.BoundingBox).

"""
Polygon(vertices::List{Point{T}}) where {T} = Polygon{T}(vertices)

Polygon{T}(vertices::List{<:Point}) where {T} =
    Polygon{T}(map(Fix1(convert_coord_type, T), vertices))

# Build a polygon from an homogeneous tuple/vector of point-like objects.
for type in (:AbstractPoint, :(CartesianIndex{2}), :(NTuple{2,Number}))
    @eval begin
        Polygon(pnts::$type...) = Polygon(pnts)
        Polygon{T}(pnts::$type...) where {T} = Polygon{T}(pnts)
        function Polygon(pnts::List{<:$type})
            T = coord_type(point_type(pnts))
            return Polygon{T}(pnts)
        end
        function Polygon{T}(pnts::AbstractVector{<:$type}) where {T}
            # Check length before conversion.
            length(pnts) ≥ 3 || throw_insufficent_number_of_polygon_vertices(length(pnts))
            v = similar(pnts, Point{T})
            for (i, xy) in enumerate(pnts)
                v[i] = xy
            end
            return Polygon{T}(v)
        end
        function Polygon{T}(pnts::NTuple{N,$type}) where {T,N}
            # Check length before conversion.
            N ≥ 3 || throw_insufficent_number_of_polygon_vertices(N)
            return Polygon{T}(map(Point{T}, pnts))
        end
    end
end

# Tuples with zero item.
Polygon(::Tuple{}) = throw_insufficent_number_of_polygon_vertices(0)
Polygon{T}(::Tuple{}) where {T} = throw_insufficent_number_of_polygon_vertices(0)

# Build a polygon from other geometric objects.
Polygon{T}(obj::Union{Rectangle,BoundingBox}) where {T} = Polygon(convert_coord_type(T, obj))
Polygon(rect::Rectangle) = Polygon(collect(vertices(box)))
Polygon(box::BoundingBox) =
    isempty(box) ? throw(ArgumentError("cannot build a polygon from an empty bounding-box")) :
    Polygon(collect(vertices(box)))

# Convert/copy constructors.
Polygon(poly::Polygon) = poly
Polygon{T}(poly::Polygon{T}) where {T} = poly
Polygon{T}(poly::Polygon) where {T} = Polygon{T}(map(Fix1(convert_coord_type, T), elements(poly)))
Base.convert(::Type{T}, poly::T) where {T<:Polygon} = poly
Base.convert(::Type{T}, obj::PolygonLike) where {T<:Polygon} = T(obj)

# Make rectangles indexable iterators.
Base.length(        poly::Polygon) = length(elements(poly))
Base.size(          poly::Polygon{<:Any,<:AbstractVector}) = size(elements(poly))
Base.axes(          poly::Polygon{<:Any,<:AbstractVector}) = axes(elements(poly))
Base.size(          poly::Polygon{<:Any,<:Tuple}) = (length(poly),)
Base.axes(          poly::Polygon{<:Any,<:Tuple}) = (Base.OneTo(length(poly)),)
Base.firstindex(    poly::Polygon) = firstindex(elements(poly))
Base.lastindex(     poly::Polygon) = lastindex(elements(poly))
Base.eachindex(     poly::Polygon) = keys(poly)
Base.keys(          poly::Polygon) = eachindex(elements(poly))
Base.values(        poly::Polygon) = elements(poly)
Base.first(         poly::Polygon) = poly[firstindex(poly)]
Base.last(          poly::Polygon) = poly[lastindex(poly)]
Base.eltype(        poly::Polygon) = eltype(typeof(poly))
Base.IteratorSize(  poly::Polygon) = IteratorSize(typeof(poly))
Base.IteratorEltype(poly::Polygon) = IteratorEltype(typeof(poly))
Base.eltype(        ::Type{<:Polygon{T}}) where {T} = Point{T}
Base.IteratorSize(  ::Type{<:Polygon}) = HasLength()
Base.IteratorEltype(::Type{<:Polygon}) = HasEltype()
@inline Base.getindex(A::Polygon{<:Any,<:Tuple}, i::Integer) = getindex(elements(A), i)
@inline function Base.getindex(A::Polygon{<:Any,<:AbstractVector}, i::Integer)
    i = as(Int, i)
    @boundscheck checkbounds(elements(A), i)
    return @inbounds getindex(elements(A), i)
end
@inline function Base.setindex!(A::Polygon{<:Any,<:AbstractVector}, x, i::Integer)
    i = as(Int, i)
    @boundscheck checkbounds(elements(A), i)
    @inbounds setindex!(elements(A), x, i)
    return A
end

# Equality of polygons.
function ==(A::Polygon, B::Polygon)
    A === B && return true
    I = eachindex(A)
    eachindex(B) == I || return false
    @inbounds for i in I
        A[i] == B[i] || return false
    end
    return true
end

# Properties of polygons.
Base.propertynames(::Polygon) = (:vertices,)
Base.getproperty(poly::Polygon, key::Symbol) =
    key === :vertices ? elements(poly) : throw(KeyError(key))

function Base.show(io::IO, poly::Polygon{T,V}) where {T,V}
    print(io, "Polygon{", T, ",", V, "}((")
    flag = false
    for i in eachindex(poly)
        i > firstindex(poly) && print(io, "), (")
        print(io, poly[i].x, ", ", poly[i].y)
    end
    print(io, "))")
    nothing
end

@noinline throw_insufficent_number_of_polygon_vertices(len::Integer) =
    throw(ArgumentError("polygons have at least 3 vertices, got $(len)"))

"""
    vec(poly::TwoDimensional.Polygon)

yields a vector of the vertices of the polygon `poly`.

Call `collect(poly)` to make an independent copy of the vector of vertices.

"""
Base.vec(poly::Polygon{<:Any,<:AbstractVector}) = elements(poly)
function Base.vec(poly::Polygon)
    src = elements(poly)
    dst = Vector{eltype(src)}(undef, length(src))
    @inbounds for (i, x) in enumerate(src)
        dst[i] = x
    end
    return dst
end

Base.Vector(poly::Polygon{<:Any,<:Union{Tuple,Vector}}) = vec(poly)

Base.collect(poly::Polygon{<:Any,<:AbstractVector}) = copy(elements(poly))
Base.collect(poly::Polygon) = Vector(poly)

"""
    TwoDimensional.geometric_properties(poly) -> prop

yields a named tuple with the geometric properties of the polygon `poly`:

    prop.singular # `true` if successive polygon vertices are not distinct
    prop.convex   # `true` if polygon is convex
    prop.direct   # `true` if vertices are ordered with direct trigonometric
                  # orientation (anti-clockwise)

"""
function geometric_properties(poly::Polygon)
    i_first = firstindex(poly)
    i_last = lastindex(poly)
    xmin = typemax(coord_type(poly))
    orient = 0
    direct = true # direct orientation?
    convex = true # polygon is convex?
    singular = false # polygon is singular?
    @inbounds for i in eachindex(poly)
        # Compute cross-product to determine whether the polygon is convex and
        # its global direction (to guess where is interior). Determine
        # bounding-box.
        i_prev = ifelse(i > i_first, i - 1, i_last)
        i_next = ifelse(i < i_last, i + 1, i_first)
        s = cross(poly[i] - poly[i_prev], poly[i_next] - poly[i])
        if iszero(s)
            # Polygon vertices are not distinct.
            singular = true
            direct = false
            convex = false
        end
        this_orient = (s > zero(s) ? +1 : -1)
        if orient == 0
            orient = this_orient
        elseif orient != this_orient
            convex = false
        end
        x = poly[i].x
        if x < xmin
            # Set whether polygon is in direct trigonometric (anti-clockwise)
            # order according to the sign of the cross product of the two edges
            # at this extremum vertex.
            xmin = x
            direct = this_orient > 0
        end
    end
    return (convex = convex, direct = direct, singular = singular)
end

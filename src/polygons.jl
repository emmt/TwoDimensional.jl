"""
    poly = Polygon{T}(pnts...)
    poly = Polygon{T}(pnts)

construct a polygon with vertices given by point-like objects `pnts...` or
vector or tuple of point-like objects `pnts`. Parameter `T` is the type used to
store coordinates, if omitted, a common coordinate type is automatically
inferred.

The method `vec(poly)` yields the vector backing the storage of the vertices of
the polygon. Call `collect(poly)` or `copy(vec(poly))` to make an independent
copy of the vector of vertices.

Vertices are directly accessible by indexing the polygon object:

    length(poly)     # number of vertices
    firstindex(poly) # index of first vertex
    lastindex(poly)  # index of last vertex
    poly[i]          # i-th vertex (a point)
    first(poly)      # first vertex
    last(poly)       # last vertex
    vec(poly)        # vector of vertices

See also [`Point`](@ref TwoDimensional.Point), [`TwoDimensional.vertices`](@ref
TwoDimensional.vertices), and [`BoundingBox`](@ref TwoDimensional.BoundingBox).

"""
Polygon(vertices::AbstractVector{<:Point{T}}) where {T} = Polygon{T}(vertices)
Polygon{T}(vertices::AbstractVector{<:Point}) where {T} =
    Polygon{T}(convert(AbstractVector{Point{T}}, vertices))

# Build a polygon from an homogeneous tuple/vector of point-like objects.
for type in (:AbstractPoint, :(CartesianIndex{2}), :(NTuple{2,Number}))
    @eval begin
        Polygon(pnts::$type...) = Polygon(pnts)
        Polygon{T}(pnts::$type...) where {T} = Polygon{T}(pnts)
        function Polygon(pnts::TupleOrVector{$type})
            T = coord_type(point_type(pnts))
            return Polygon{T}(pnts)
        end
        function Polygon{T}(pnts::TupleOrVector{$type}) where {T}
            # Check length before conversion.
            len = length(pnts)
            len ≥ 3 || throw_insufficent_number_of_polygon_vertices(len)
            v = Vector{Point{T}}(undef, len)
            for (i, xy) in enumerate(pnts)
                v[i] = xy
            end
            return Polygon{T}(v)
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
Polygon{T}(poly::Polygon) where {T} = Polygon{T}(map(Fix1(convert_coord_type, T), vec(poly)))
Base.convert(::Type{T}, poly::T) where {T<:Polygon} = poly
Base.convert(::Type{T}, obj::PolygonLike) where {T<:Polygon} = T(obj)

# Make rectangles indexable iterators.
Base.length(        poly::Polygon) = length(vec(poly))
Base.size(          poly::Polygon) = size(vec(poly))
Base.axes(          poly::Polygon) = axes(vec(poly))
Base.firstindex(    poly::Polygon) = firstindex(vec(poly))
Base.lastindex(     poly::Polygon) = lastindex(vec(poly))
Base.eachindex(     poly::Polygon) = keys(poly)
Base.keys(          poly::Polygon) = eachindex(vec(poly))
Base.values(        poly::Polygon) = vec(poly)
Base.first(         poly::Polygon) = poly[firstindex(poly)]
Base.last(          poly::Polygon) = poly[lastindex(poly)]
Base.eltype(        poly::Polygon) = eltype(typeof(poly))
Base.IteratorSize(  poly::Polygon) = IteratorSize(typeof(poly))
Base.IteratorEltype(poly::Polygon) = IteratorEltype(typeof(poly))
Base.eltype(        ::Type{<:Polygon{T}}) where {T} = Point{T}
Base.IteratorSize(  ::Type{<:Polygon}) = HasLength()
Base.IteratorEltype(::Type{<:Polygon}) = HasEltype()
@inline function Base.getindex(poly::Polygon, i::Integer)
    i = as(Int, i)
    v = vec(poly)
    @boundscheck checkbounds(v, i)
    r = @inbounds getindex(v, i)
    return r
end
@inline function Base.setindex!(poly::Polygon, x::PointLike, i::Integer)
    i = as(Int, i)
    v = vec(poly)
    @boundscheck checkbounds(v, i)
    @inbounds setindex!(v, x, i)
    return poly
end

# Properties of polygons.
Base.propertynames(::Polygon) = (:vertices,)
Base.getproperty(poly::Polygon, key::Symbol) =
    key === :vertices ? vec(poly) : throw(KeyError(key))

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

yields the vector backing the storage of the vertices of the polygon `poly`.

Call `collect(poly)` or `copy(vec(poly))` to make an independent copy of the
vector of vertices.

"""
Base.vec(poly::Polygon) = getfield(poly, 1)
Base.collect(poly::Polygon) = copy(vec(poly))

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

"""
    TwoDimensional.is_convex(poly) -> bool

yields whether the polygon `poly` is strictly convex.

"""
is_convex(poly::Polygon) = geometric_properties(poly).convex
